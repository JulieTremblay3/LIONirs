function varargout = GUI_GLMidentification(varargin)
% GUI_GLMIDENTIFICATION M-file for GUI_GLMidentification.fig
%      GUI_GLMIDENTIFICATION, by itself, creates a new GUI_GLMIDENTIFICATION or raises the existing
%      singleton*.
%
%      H = GUI_GLMIDENTIFICATION returns the handle to a new GUI_GLMIDENTIFICATION or the handle to
%      the existing singleton*.
%
%      GUI_GLMIDENTIFICATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GLMIDENTIFICATION.M with the given input arguments.
%
%      GUI_GLMIDENTIFICATION('Property','Value',...) creates a new GUI_GLMIDENTIFICATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_GLMidentification_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_GLMidentification_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_GLMidentification

% Last Modified by GUIDE v2.5 10-Jan-2019 15:38:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_GLMidentification_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_GLMidentification_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_GLMidentification is made visible.
function GUI_GLMidentification_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_GLMidentification (see VARARGIN)

% Choose default command line output for GUI_GLMidentification
handles.output = hObject;
tstart = varargin{1};
tstop = varargin{2};
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
d = PMI{currentsub}.data(cf).HRF.AvgC;

indt = tstart(end):tstop(end);
intensnorm = d(indt,:);
 set(handles.text_tstart,'string', ['Time :' num2str(PMI{currentsub}.data(cf).HRF.tHRF(tstart)), ' to ',num2str(PMI{currentsub}.data(cf).HRF.tHRF(tstop))] ) 
%Detrent DATA segment for centrering 
X = 1:1:size(intensnorm,1);
sous = intensnorm(end,:) - intensnorm(1,:);
Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
Mb2 =  intensnorm(1,:)'; %offset    
A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
spar = intensnorm - A; %with detrending 
spar = intensnorm;     %without detrending
listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
MeasListActplotLst = zeros(size(A,2),1);
MeasListActplotLst(PMI{1}.plotLst,:)=1;
if min(PMI{1}.plotLst-size(A,2)/2) >0
    MeasListActplotLst(PMI{1}.plotLst-size(A,2)/2,:)=1;
end
MeasListActplotLst((size(A,2)/2+1):end,1)=0; %forcer la fin a zeros
listgood = find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));

isfield( PMI{currentsub}.data(cf),'AUX') 
PMI{currentsub}.data(cf).AUX

t = PMI{currentsub}.data(cf).HRF.tHRF;
% spar = cat(4,spar(:,:,1:end/2),spar(:,:,end/2+1:end));
% spartmp = spar(:,:,listgood,:);
PMI{1}.tmpGLM.selected =1  ;   %used to visualised
PMI{1}.tmpGLM.spar = spar;                    %Data initial
PMI{1}.tmpGLM.listgood = listgood ; 
PMI{1}.tmpGLM.indt = [tstart(end),tstop(end)];%Time indice
load(PMI{1}.NIRSmat{1});
iAUX = 1 ;%each aux file could have a different sample rate
fsNIRS=NIRS.Cf.dev.fs;
idfile=PMI{1}.idfile;
iRegressor =  2;
if isfield(PMI{currentsub}.tmpGLM,'AUX')
    rmfield(PMI{currentsub}.tmpGLM,'AUX');
end

for iAUX = 1:numel(NIRS.Dt.AUX)
    nameAUX = NIRS.Dt.AUX(iAUX).pp(end).p{idfile};
    if isfield(NIRS.Dt.AUX(iAUX).pp(end),'sync_timesec')
        tstartf = NIRS.Dt.AUX(iAUX).pp(end).sync_timesec{idfile};
    else 
        tstartf = 1;
    end
    tstopf = tstartf+PMI{currentsub}.data(cf).HRF.tHRF(end);
    [pathtmp,filetmp,exttmp]=fileparts(nameAUX);
    try
    [data,infoBV,label,ind_dur_ch] = fopen_EEG(nameAUX, tstartf, tstopf);
      for ich=1:numel(infoBV.name_ele)
            PMI{currentsub}.tmpGLM.AUX.label{iRegressor} =[filetmp,' ',infoBV.name_ele{ich}];           
            fsAUX =1/(infoBV.SamplingInterval/1000000); %Frequence echantillonage Hz
            tmp = data(:,ich);
            [p,q] = rat(fsNIRS/fsAUX,0.0001);
            %tmpr=resample( tmp , p, q); %ancien Julie - a cause
                    %d'un jump au début fin, LCD l'a modifie par
                    %downsample - environ le meme fonctionnement - voir la
                    %photo Guide DownsampleEEgdata dans GITHUB student
                     tmpr=downsample( tmp , q);
            % we resample aux to the data initial size
            
            %%addition to the script: normalize the AUX data LCD
            tmpr=(tmpr-mean(tmpr))/std(tmpr);
            
            if numel(tmpr)<numel(t)
                nplus = numel(t)-numel(tmpr);
                tmpr = [tmpr;tmpr(end-nplus:end) ];
            elseif numel(tmpr)>numel(t)
                tmpr = tmpr(1:numel(t));
            end
            %we cut to fit the spar selection                      
            PMI{currentsub}.tmpGLM.AUX.fs{iRegressor} = fsNIRS;
            PMI{currentsub}.tmpGLM.AUX.data{iRegressor} = tmpr(PMI{1}.tmpGLM.indt(1):PMI{1}.tmpGLM.indt(end));
            clear tmpr
            PMI{currentsub}.tmpGLM.AUX.view{iRegressor} = 1;
            iRegressor = iRegressor +1;
      end
    catch
        
    end
end

PMI{1}.tmpGLM.AUX.label{1} = 'Constant';
PMI{1}.tmpGLM.AUX.fs{1} = fsNIRS;
PMI{1}.tmpGLM.AUX.data{1} = ones(numel(PMI{1}.tmpGLM.indt(1):PMI{1}.tmpGLM.indt(end) ),1);
PMI{1}.tmpGLM.AUX.view{1} = 1;

        set(handles.listbox_auxilary,'string',PMI{currentsub}.tmpGLM.AUX.label');
        set(handles.listbox_auxilary,'enable','on');
%load and resample the auxiliary for use in the GLM 

% % handles.NIRS.Dt.AUX(iAUX).pp(end).p{idfile}
%     if isfield(handles.NIRS.Dt,'AUX') %LOAD AUX ICI AUXILIARY AND ADD IN PMI structure
%         iAUX = get(handles.popupauxnb,'value')
%         PMI{currentsub}.data(cf).AUX = [];
%         nameAUX = handles.NIRS.Dt.AUX(iAUX).pp(end).p{idfile};
%         [data,infoBV,label,ind_dur_ch] = fopen_EEG(nameAUX);
%         for ich=1:numel(infoBV.name_ele)
%             PMI{currentsub}.data(cf).AUX.label{ich} =[infoBV.name_ele{ich},'  on'];
%             PMI{currentsub}.data(cf).AUX.data{ich} = data(:,ich);
%             PMI{currentsub}.data(cf).AUX.fs{ich} =1/(infoBV.SamplingInterval/1000000); %Frequence echantillonage Hz
%             PMI{currentsub}.data(cf).AUX.view{ich} = 1;
%         end
%         set(handles.listbox_AUXCH,'string',PMI{currentsub}.data(cf).AUX.label')
%         set(handles.listbox_AUXCH,'enable','on')
%     end
%     
% 
% load

PMI{currentsub}.data(cf).AUX
set(guiHOMER,'UserData',PMI);
% btn_MultipleLinearRegression_Callback(hObject, eventdata, handles)
% nbcomp = str2num(get(handles.edit_nbComponent,'string'));
% label = [];
% for i=1:nbcomp
%     label = [label,{['F',num2str(i)]}];
% end
% set(handles.popupmenu_Factor,'string',label)



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_GLMidentification wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_GLMidentification_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_listAuxiliary_Callback(hObject, eventdata, handles)
% hObject    handle to edit_listAuxiliary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_listAuxiliary as text
%        str2double(get(hObject,'String')) returns contents of edit_listAuxiliary as a double


% --- Executes during object creation, after setting all properties.
function edit_listAuxiliary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_listAuxiliary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_MultipleLinearRegression.
function btn_MultipleLinearRegression_Callback(hObject, eventdata, handles)
% hObject    handle to btn_MultipleLinearRegression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1;
handlesoutput = hObject;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
y=PMI{1}.tmpGLM.spar; 
Regressorlist = get(handles.listbox_AUXselected,'string');
if ~iscell(Regressorlist)
    Regressorlist = {Regressorlist};
end
%find indice of regressor in the list
idreg = [];
for iReg = 1:numel(Regressorlist)
    for ilist = 1:numel(PMI{1}.tmpGLM.AUX.label)
        if strcmp(PMI{1}.tmpGLM.AUX.label{ilist},Regressorlist{iReg})
        idreg = [idreg,ilist];
        end
    end
    
end
PMI{1}.tmpGLM.idreg = idreg; %regression auxiliary indice used 
disp('idreg = ');
disp(idreg);
X = [];%ones(size(PMI{1}.tmpGLM.spar,1),1);
for ireg = 1:numel(idreg)
    X = [X,    PMI{1}.tmpGLM.AUX.data{idreg(ireg)}];
end
if isempty(X)
    errordlg('The variable X is empty. Use the GUI to correct it.');
    return
end
%figure;plot(X)
%figure;plot(PMI{1}.tmpGLM.AUX.data{idreg(ireg)})
beta = zeros(size(X,2),numel(PMI{1}.tmpGLM.listgood));
bstd = zeros(size(X,2),numel(PMI{1}.tmpGLM.listgood));
for ich = 1:numel(PMI{1}.tmpGLM.listgood)
    idch = PMI{1}.tmpGLM.listgood(ich);
    y = PMI{1}.tmpGLM.spar(:,idch);
   %figure;plot(PMI{1}.tmpGLM.spar)
    if sum(isnan(y)) 
        b = zeros(size(X,2),1);
        beta(:,ich) = 0;
        bstd(:,ich) = 0;

    else
        [b,bint,r,rint,stats]=  regress(y,X);
        beta(:,ich) = b;
        bstd(:,ich) = stats(4);
    end
end
PMI{1}.tmpGLM.beta = beta;
PMI{1}.tmpGLM.std = bstd;
% PMI{1}.tmpGLM.selected = 1;
% % PMI{1}.tmpGLM.Xm = 
% % PMI{1}.tmpGLM.data = 
% set(handles.listbox_AUXselected,'value',1) %PMI{1}.tmpGLM.selected );
set(guiHOMER,'UserData',PMI);

% --- Executes on selection change in popup_factor.
function popup_factor_Callback(hObject, eventdata, handles)
% hObject    handle to popup_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_factor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_factor


% --- Executes during object creation, after setting all properties.
function popup_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_auxilary.
function listbox_auxilary_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_auxilary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_auxilary contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_auxilary


% --- Executes during object creation, after setting all properties.
function listbox_auxilary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_auxilary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_get_regressor.
function btn_get_regressor_Callback(hObject, eventdata, handles)
% hObject    handle to btn_get_regressor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listchoose =get(handles.listbox_AUXselected,'string');
if isempty(listchoose)
    listchoose = [];
elseif  ~iscell(listchoose)
    listchoose = {listchoose};
end
listaux = get(handles.listbox_auxilary,'string');
valaux = get(handles.listbox_auxilary,'value'); 
listchoose = [listchoose;listaux{valaux}];
set(handles.listbox_AUXselected,'string',listchoose)

% --- Executes on button press in btn_removeregressor.
function btn_removeregressor_Callback(hObject, eventdata, handles)
% hObject    handle to btn_removeregressor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listchoose = get(handles.listbox_AUXselected,'string');
valchoose = get(handles.listbox_AUXselected,'value'); 
listchoose(valchoose) = [];
set(handles.listbox_AUXselected,'string',listchoose);
set(handles.listbox_AUXselected,'value',1); 

% --- Executes on selection change in listbox_AUXselected.
function listbox_AUXselected_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_AUXselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_AUXselected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_AUXselected
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
PMI{1}.tmpGLM.selected = get(handles.listbox_AUXselected,'value');
set(guiHOMER,'UserData',PMI);


% --- Executes during object creation, after setting all properties.
function listbox_AUXselected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_AUXselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function  varargout = resample(varargin)
%RESAMPLE  Resample uniform or nonuniform data to a new fixed rate.
%   Y = RESAMPLE(X,P,Q) resamples the values, X, of a uniformly sampled
%   signal at P/Q times the original sample rate using a polyphase
%   anti-aliasing filter. If X is a matrix, then RESAMPLE treats each
%   column as an independent channel.
%
%   In its filtering process, RESAMPLE assumes the samples at times before
%   and after the given samples in X are equal to zero. Thus large
%   deviations from zero at the end points of the sequence X can cause
%   inaccuracies in Y at its end points.
%
%   [Y,Ty] = RESAMPLE(X,Tx) resamples the values, X, of a signal sampled at
%   the instants specified in vector Tx. RESAMPLE interpolates X linearly
%   onto a vector of uniformly spaced instants, Ty, with the same endpoints
%   and number of samples as Tx.  Tx may either be a numeric vector
%   expressed in seconds or a datetime object.  NaNs and NaTs (for datetime
%   objects) are treated as missing data and are ignored.
%
%   [Y,Ty] = RESAMPLE(X,Tx,Fs) uses interpolation and an anti-aliasing
%   filter to resample the signal at a uniform sample rate, Fs, expressed
%   in hertz (Hz).
% 
%   [Y,Ty] = RESAMPLE(X,Tx,Fs,P,Q) interpolates X to an intermediate
%   uniform grid with sample rate equal Q*Fs/P and filters the result
%   using UPFIRDN to upsample by P and downsample by Q.  Specify P and Q
%   so that Q*Fs/P is least twice as large as the highest frequency in the
%   input signal.  
%
%   [Y,Ty] = RESAMPLE(X,Tx,...,METHOD) specifies the interpolation method.
%   The default is linear interpolation.  Available methods are:
%     'linear' - linear interpolation
%     'pchip'  - shape-preserving piecewise cubic interpolation
%     'spline' - piecewise cubic spline interpolation
%
%   Y = RESAMPLE(...,P,Q,N) uses a weighted sum of 2*N*max(1,Q/P) samples
%   of X to compute each sample of Y.  The length of the FIR filter
%   RESAMPLE applies is proportional to N; by increasing N you will get
%   better accuracy at the expense of a longer computation time.  If you
%   don't specify N, RESAMPLE uses N = 10 by default.  If you let N = 0,
%   RESAMPLE performs a nearest neighbor interpolation; that is, the output
%   Y(n) is X(round((n-1)*Q/P)+1) ( Y(n) = 0 if round((n-1)*Q/P)+1 >
%   length(X) ).
%
%   Y = RESAMPLE(...,P,Q,N,BTA) uses BTA as the BETA design parameter for
%   the Kaiser window used to design the filter.  RESAMPLE uses BTA = 5 if
%   you don't specify a value.
%
%   Y = RESAMPLE(...,P,Q,B) uses B to filter X (after upsampling) if B is a
%   vector of filter coefficients.  RESAMPLE assumes B has odd length and
%   linear phase when compensating for the filter's delay; for even length
%   filters, the delay is overcompensated by 1/2 sample.  For non-linear
%   phase filters consider using UPFIRDN.
%
%   [Y,B] = RESAMPLE(X,P,Q,...) returns in B the coefficients of the filter
%   applied to X during the resampling process (after upsampling).
%
%   [Y,Ty,B] = RESAMPLE(X,Tx,...) returns in B the coefficients of the
%   filter applied to X during the resampling process (after upsampling).
%
%   % Example 1:
%   %   Resample a sinusoid at 3/2 the original rate.
%
%   tx = 0:3:300-3;         % Time vector for original signal
%   x = sin(2*pi*tx/300);   % Define a sinusoid 
%   ty = 0:2:300-2;         % Time vector for resampled signal        
%   y = resample(x,3,2);    % Change sampling rate
%   plot(tx,x,'+-',ty,y,'o:')
%   legend('original','resampled');
%   xlabel('Time')
%
%   % Example 2:
%   %   Resample a non-uniformly sampled sinusoid to a uniform 50 Hz rate.
%   
%   Fs = 50;
%   tx = linspace(0,1,21) + .012*rand(1,21);
%   x = sin(2*pi*tx);
%   [y, ty] = resample(x, tx, Fs);
%   plot(tx,x,'+-',ty,y,'o:')
%   legend('original','resampled');
%   xlabel('Time')
%   axis tight
%
%   See also UPFIRDN, INTERP, INTERP1, DECIMATE, FIRLS, KAISER, INTFILT.

%   NOTE-1: digital anti-alias filter is designed via windowing

%   Author(s): James McClellan, 6-11-93
%              Modified to use upfirdn, T. Krauss, 2-27-96
%   Copyright 1988-2018 The MathWorks, Inc.
%     

narginchk(2,8);

[varargin{:}] = convertStringsToChars(varargin{:});

% fetch the users desired interpolation method (if any)
[method, varargin] = getInterpMethod(varargin);

if isscalar(varargin{2})
  % [...] = RESAMPLE(X,P,Q,...)
  if nargin < 3
      error(message('signal:resample:MissingQ'));
  end
  % interpolation is not performed when using uniformly sampled data
  if ~strcmp(method,'')
    error(message('signal:resample:UnexpectedInterpolation',method));
  end  
  [varargout{1:max(1,nargout)}] = uniformResample(varargin{:});
else
  % [...] = RESAMPLE(X,Tx,...)
  if strcmp(method,'')
    % use linear method by default
    method = 'linear';
  end
  [varargout{1:max(1,nargout)}] = nonUniformResample(method,varargin{:});
end


function [y, ty, h] = nonUniformResample(method,varargin)

% fetch non-NaN input samples in time-sorted order
[x, tx] = getSamples(varargin{1:2});
  
% obtain sample rate
if numel(varargin)>2
  fs = varargin{3};
  validateFs(fs);
elseif isdatetime(tx)
  % compute the average sample rate if unspecified
  fs = (numel(tx)-1) / seconds(tx(end)-tx(1));
else
  fs = (numel(tx)-1) / (tx(end)-tx(1));
end

% get rational sample ratio  
if numel(varargin)>4
  p = varargin{4};
  q = varargin{5};
  validateResampleRatio(p,q);
elseif numel(varargin)>3
  error(message('signal:resample:MissingQ'));
else
  % get a close rational approximation of the ratio between the desired
  % sampling rate and the average sample rate
  [p, q] = getResampleRatio(tx,fs);
end

% use cubic spline interpolation onto a uniformly sampled grid
% with a target sample rate of Q*Fs/P
tsGrid = p/(q*fs);
if isdatetime(tx)
  tGrid = tx(1):seconds(tsGrid):tx(end);
else
  tGrid = tx(1):tsGrid:tx(end);
end


if isreal(x)
  xGrid = matInterp1(tx,x,tGrid,method);
else
  % compute real and imaginary channels independently
  realGrid = matInterp1(tx,real(x),tGrid,method);
  imagGrid = matInterp1(tx,imag(x),tGrid,method);
  xGrid = complex(realGrid,imagGrid);
end

% recover the desired sampling rate by resampling the grid at the 
% specified ratio
[y, h] = uniformResample(xGrid,p,q,varargin{6:end});

% create output time vector
if isvector(y)
  ny = length(y);
else
  ny = size(y,1);
end

if isdatetime(tx)
  ty = tx(1) + seconds((0:ny-1)/fs);
else
  ty = tx(1) + (0:ny-1)/fs;
end  

% match dimensionality of output time vector to input
if iscolumn(tx)
  ty = ty(:);
end 

%-------------------------------------------------------------------------
function  [y, h] = uniformResample( x, p, q, N, bta )

if nargin < 5
    bta = 5;
end   %--- design parameter for Kaiser window LPF
if nargin < 4 
    N = 10;
end

validateResampleRatio(p, q);

[p,q] = rat( p/q, 1e-12 );  %--- reduce to lowest terms 
   % (usually exact, sometimes not; loses at most 1 second every 10^12 seconds)
if (p==1) && (q==1)
    y = x; 
    h = 1;
    return
end
pqmax = max(p,q);
if length(N)>1      % use input filter
   L = length(N);
   h = N;
else                % design filter
   if( N>0 )
      fc = 1/2/pqmax;
      L = 2*N*pqmax + 1;
      h = firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,bta)' ;
      h = p*h/sum(h);
   else
      L = p;
      h = ones(1,p);
   end
end

Lhalf = (L-1)/2;
isvect = any(size(x)==1);
if isvect
    Lx = length(x);
else
    Lx = size(x, 1);
end

% Need to delay output so that downsampling by q hits center tap of filter.
nz = floor(q-mod(Lhalf,q));
z = zeros(1,nz);
h = [z h(:).'];  % ensure that h is a row vector.
Lhalf = Lhalf + nz;

% Number of samples removed from beginning of output sequence 
% to compensate for delay of linear phase filter:
delay = floor(ceil(Lhalf)/q);

% Need to zero-pad so output length is exactly ceil(Lx*p/q).
nz1 = 0;
while ceil( ((Lx-1)*p+length(h)+nz1 )/q ) - delay < ceil(Lx*p/q)
    nz1 = nz1+1;
end
h = [h zeros(1,nz1)];

% ----  HERE'S THE CALL TO UPFIRDN  ----------------------------
y = upfirdn(x,h,p,q);

% Get rid of trailing and leading data so input and output signals line up
% temporally:
Ly = ceil(Lx*p/q);  % output length
% Ly = floor((Lx-1)*p/q+1);  <-- alternately, to prevent "running-off" the
%                                data (extrapolation)
if isvect
    y(1:delay) = [];
    y(Ly+1:end) = [];
else
    y(1:delay,:) = [];
    y(Ly+1:end,:) = [];
end

h([1:nz (end-nz1+1):end]) = [];  % get rid of leading and trailing zeros 
                                 % in case filter is output

%-------------------------------------------------------------------------
function [x, tx] = removeMissingTime(x,tx)
if isdatetime(tx)
  idx = isnat(tx);
else
  idx = isnan(tx);
end

if ~isempty(idx)
  tx(idx) = [];
  if isvector(x)
    x(idx) = [];
  else
    x(idx,:) = [];
  end
end

%-------------------------------------------------------------------------
function y = matInterp1(tin, xin, tout, method)

% by default specify tout as a column to obtain column matrix output
tout = tout(:);

if isvector(xin)
  if isrow(xin)
    % preserve orientation of input vector
    tout = tout.';
  end
  % interpolate, excluding NaN
  idx = find(~isnan(xin));  
  y = vecInterp1(tin(idx), xin(idx), tout, method);
else
  % initialize matrix output
  nRows = size(tout,1);
  nCols = size(xin,2);
  y = zeros(nRows,nCols);
  
  % loop through each column of input x
  for col=1:nCols
    % interpolate, excluding NaN
    idx = find(~isnan(xin(:,col)));  
    y(:,col) = vecInterp1(tin(idx), xin(idx,col), tout, method);
  end
end

%-------------------------------------------------------------------------
function y = vecInterp1(tin, xin, tout, method)

% check sample times for duplicate entries
iDup = find(diff(tin)==0);

% copy indices to point to the repeated locations in xin/tin
iRepeat = 1 + iDup;

while ~isempty(iDup)
  % find the number of successive equal sample times
  numEqual = find(diff(iDup)~=1,1,'first');
  if isempty(numEqual)
    numEqual = numel(iDup);
  end
  
  % replace leading x with mean value of all duplicates
  xSelect = xin(iDup(1) + (0:numEqual));
  xMean = mean(xSelect(~isnan(xSelect)));
  xin(iDup(1)) = xMean;

  % move to next block of conflicting sample times
  iDup = iDup(2+numEqual:end);
end

% remove duplicates
xin(iRepeat) = [];
tin(iRepeat) = [];

% call interp
y = interp1(tin, xin, tout, method, 'extrap');

%-------------------------------------------------------------------------
function [x, tx] = getSamples(x, tx)

validateattributes(x, {'numeric'},{'2d'}, ...
    'resample','X',1);
  
if isvector(x) && numel(x)~=numel(tx)
  error(message('signal:resample:TimeVectorMismatch'));
end

if ~isvector(x) && size(x,1)~=numel(tx)
  error(message('signal:resample:TimeRowMismatch'));
end

if isdatetime(tx)
  validateattributes(tx, {'datetime'},{'vector'}, ...
      'resample','Tx',2);
else
  validateattributes(tx, {'numeric'},{'real','vector'}, ...
      'resample','Tx',2);
end

[x, tx] = removeMissingTime(x, tx);

% for efficiency, place samples in time-sorted order
[tx, idx] = sort(tx);
if isvector(x)
  % handle row vector input
  x = x(idx);
else
  x = x(idx,:);
end

% check finite behavior after sorting and NaN removal
if ~isdatetime(tx)
  validateattributes(tx, {'numeric'},{'finite'}, ...
      'resample','Tx',2);  
end

%-------------------------------------------------------------------------
function validateFs(fs)
validateattributes(fs, {'numeric'},{'real','finite','scalar','positive'}, ...
    'resample','Fs',3);

%-------------------------------------------------------------------------
function validateResampleRatio(p, q)
validateattributes(p, {'numeric'},{'integer','positive','finite','scalar'}, ...
    'resample','P');
validateattributes(q, {'numeric'},{'integer','positive','finite','scalar'}, ...
    'resample','Q');

%-------------------------------------------------------------------------
function [p, q] = getResampleRatio(t, fs)

% compute the average sampling interval
tsAvg = (t(end)-t(1))/(numel(t)-1);
if isduration(tsAvg)
  tsAvg = seconds(tsAvg);
end

% get a rational approximation of the ratio of the desired to the average
% sample rate

[p, q] = rat(tsAvg*fs,.01);

if p < 2
  % sample rate too small for input
  p = 1;
  q = round(1/(tsAvg*fs));
end

%-------------------------------------------------------------------------
function [method, arglist] = getInterpMethod(arglist)

method = '';
supportedMethods = {'linear','pchip','spline'};

iFound = 0;

for i=1:numel(arglist)
  if ischar(arglist{i})
    method = validatestring(arglist{i},supportedMethods,'resample','METHOD');
    iFound = i;
  end
end

if iFound
  arglist(iFound) = [];
end