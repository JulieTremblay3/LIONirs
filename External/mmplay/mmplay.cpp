#include "mmFilter.h"

// unfortunately to detect memory leaks, we have to put code in each .cpp file
#if defined(_DEBUG) && defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif

#include "mex.h"
#include "mmutils.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	HRESULT hr;
	int hasVideo = 0;
	int fullscreen = 0;

	if (nrhs < 1) mexErrMsgTxt("mmplay must have at least one parameter");
	
	for (int i=0; i< nrhs; i++) 
	{
		if ((!mxIsStruct(prhs[i]) || mxGetFieldNumber(prhs[i],"times") == -1 ||
			!((mxGetFieldNumber(prhs[i],"frames") != -1 && mxGetFieldNumber(prhs[i],"width") != -1 && mxGetFieldNumber(prhs[i],"height") != -1) ||
			(mxGetFieldNumber(prhs[i],"data") != -1 && mxGetFieldNumber(prhs[i],"rate") != -1))) &&
			!mxIsChar(prhs[i]))
				mexErrMsgTxt("mmplay accepts only structs that were created with the format of mmread");
	}

	// Create the graph builder
	IGraphBuilder* pGraphBuilder = NULL;
	if (FAILED(hr = ::CoCreateInstance(CLSID_FilterGraph, NULL, CLSCTX_INPROC_SERVER, IID_IGraphBuilder, (void**)&pGraphBuilder))) mexErrMsgTxt(message(hr));

	// make sure everything is back to "normal"; cleanUp() also releases the current pGraphBuilder, so we need to make a new one
	cleanUp(&pGraphBuilder);
	if (FAILED(hr = ::CoCreateInstance(CLSID_FilterGraph, NULL, CLSCTX_INPROC_SERVER, IID_IGraphBuilder, (void**)&pGraphBuilder))) mexErrMsgTxt(message(hr));

	mmFilter* mmfilt = new mmFilter();
	pGraphBuilder->AddFilter(mmfilt,L"mmfilt");

	for (int i=0; i< nrhs; i++) 
	{
		if (mxIsStruct(prhs[i]))
		{
			if (mxGetFieldNumber(prhs[i],"frames") != -1 && mxGetFieldNumber(prhs[i],"width") != -1 && mxGetFieldNumber(prhs[i],"height") != -1)
			{ // video...
				int nr = mxGetNumberOfElements(prhs[i]);
				for (int n=0;n<nr;n++)
				{
					mxArray* framesStruct = mxGetField(prhs[i],n,"frames");
					
					if (!mxIsStruct(framesStruct)) mexErrMsgTxt("video struct, field frames is not a struct as it should be.");
					
					int nrFrames = mxGetNumberOfElements(framesStruct);

					if (nrFrames==0) continue;
					if (mxGetFieldNumber(prhs[i],"times") == -1) mexErrMsgTxt("video struct, field 'times' is missing.");

					char** frames = (char**)newadd(new char*[nrFrames]);
					double* times = mxGetPr(mxGetField(prhs[i],n,"times"));
					int width = mxGetScalar(mxGetField(prhs[i],n,"width"));
					int height = mxGetScalar(mxGetField(prhs[i],n,"height"));
					
					if (mxGetNumberOfElements(mxGetField(prhs[i],n,"times")) != nrFrames) mexErrMsgTxt("the times vector doesn't have the same number of entries as the frames struct.");
					
					for (int f=0;f<nrFrames;f++)
					{
						mxArray* cdata = mxGetField(framesStruct,f,"cdata");
						frames[f] = (char*)mxGetPr(cdata);
						if (mxGetNumberOfElements(cdata) != width*height*3) mexErrMsgTxt("error one of the video frames, isn't of the proper size.");
					}
					
					mmfilt->addVideo(frames, times, nrFrames, width, height);
					hasVideo++;
				}
			} else if (mxGetFieldNumber(prhs[i],"data") != -1 && mxGetFieldNumber(prhs[i],"rate") != -1) {
				// audio...
				int nr = mxGetNumberOfElements(prhs[i]);
				for (int n=0;n<nr;n++)
				{
					int rate = mxGetScalar(mxGetField(prhs[i],n,"rate"));
					mxArray* data = mxGetField(prhs[i],n,"data");
					int nrChannels = mxGetN(data);
					int nrSamples = mxGetM(data);
					int nrFrames = nrSamples/rate+1; if ((nrFrames-1)*rate == nrSamples) nrFrames--;

					if (nrFrames==0) continue;
					if (mxGetFieldNumber(prhs[i],"times") == -1) mexErrMsgTxt("audio struct, field 'times' is missing.");

					double** frames = (double**)newadd(new double*[nrFrames]);
					int* lens = (int*)newadd(new int[nrFrames]);
					double* datap = mxGetPr(data);
					double* times = (double*)newadd(new double[nrFrames]);
					double* ftimes = mxGetPr(mxGetField(prhs[i],n,"times"));
					
					for (int f=0;f<nrFrames;f++)
					{
						lens[f] = (f==nrFrames-1)?nrSamples-rate*f:rate;
						
						frames[f] = (double*)newadd(new double[nrChannels*lens[f]]);
						for (int j=0;j<nrChannels;j++)
							for (int k=0;k<lens[f];k++)
								frames[f][j+k*nrChannels] = datap[f*rate+k+j*nrSamples];
						times[f] = f+ftimes[0];
					}
					
					mmfilt->addAudio((char**)frames,times,nrFrames,lens,nrChannels,rate);
				}
			}
		} else if (mxIsChar(prhs[i])) {
			char str[100];
			mxGetString(prhs[i],str,100);
			if (strcmp(str,"fullscreen")==0) fullscreen = 1;
			else mexWarnMsgTxt("Unrecognized string input");
		}
	}
	
	IMediaControl* pMediaControl;
	IMediaEvent* pMediaEventEx;

	pGraphBuilder->QueryInterface(IID_IMediaControl, (void **)&pMediaControl);
	pGraphBuilder->QueryInterface(IID_IMediaEvent, (void **)&pMediaEventEx);

	IPin* pin = NULL;
	IEnumPins* pinList;
	ULONG tmp;

	//get the output pins so that we can render them...
	if (SUCCEEDED(mmfilt->EnumPins(&pinList)))
	{
		pinList->Reset();
		while (pinList->Next(1, &pin, &tmp) == S_OK)
		{
			if (FAILED(hr = pGraphBuilder->Render(pin))) mexErrMsgTxt(message(hr));
			pin->Release();
		}
		pinList->Release();
	}

	if (hasVideo)
	{
		IVideoWindow *pVidWin = NULL;
		pGraphBuilder->QueryInterface(IID_IVideoWindow, (void **)&pVidWin);

		if (fullscreen)
		{
			pVidWin->put_FullScreenMode(OATRUE);
		} else {
			mexEvalString("['Figure ' num2str(gcf)] drawnow;");
			char str[100];
			mxGetString(mexGetVariable("base","ans"),str,100);
	
			HWND hwnd = FindWindow(NULL,str);
			pVidWin->put_Owner((OAHWND)hwnd);
			pVidWin->put_WindowStyle(WS_CHILD | WS_CLIPSIBLINGS);
	
			RECT r;
			GetClientRect(hwnd, &r);
			pVidWin->SetWindowPosition(0, 0, r.right, r.bottom);
		}
		pVidWin->Release();
	}

	if (FAILED(hr = pMediaControl->Run())) mexErrMsgTxt(message(hr));

	long evCode = 0;
	while (pMediaEventEx->WaitForCompletion(1000000, &evCode) == E_ABORT);

	// make sure everything has really stopped before returning.
	if (FAILED(hr = pMediaControl->Stop())) mexErrMsgTxt(message(hr));

	pMediaEventEx->Release();
	pMediaControl->Release();
		
	cleanUp(&pGraphBuilder);
	
//	delete mmfilt; // this is automatically deleted???
}