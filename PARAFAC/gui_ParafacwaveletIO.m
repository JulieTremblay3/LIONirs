function varargout = gui_ParafacwaveletIO(varargin)
% GUI_PARAFACIO M-file for gui_ParafacIO.fig
%      GUI_PARAFACIO, by itself, creates a new GUI_PARAFACIO or raises the existing
%      singleton*.
%
%      H = GUI_PARAFACIO returns the handle to a new GUI_PARAFACIO or the handle to
%      the existing singleton*.
%
%      GUI_PARAFACIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PARAFACIO.M with the given input arguments.
%
%      GUI_PARAFACIO('Property','Value',...) creates a new GUI_PARAFACIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_ParafacIO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_ParafacIO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_ParafacIO

% Last Modified by GUIDE v2.5 09-Jan-2019 12:08:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_ParafacwaveletIO_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_ParafacwaveletIO_OutputFcn, ...
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

 
% --- Executes just before gui_ParafacIO is made visible.
function gui_ParafacwaveletIO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_ParafacIO (see VARARGIN)

% Choose default command line output for gui_ParafacIO


%CALCUL INTERACTIF PARAFAC SELECTIONNNER POUR LA FENETRE DE TEMPS DANS LE
%GUI
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
  
%Detrent DATA segment for centrering 
X = 1:1:size(intensnorm,1);
Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
Mb2 =  intensnorm(1,:)'; %offset    
A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
spar = intensnorm - A;

width = 8
Fs = 1/(PMI{1}.data.HRF.tHRF(2)-PMI{1}.data.HRF.tHRF(1))
freqvec = 0:0.1:Fs/2 
spar = morletcwtEd(A,freqvec,Fs,width);
figure
ch = 40
imagesc(PMI{1}.data.HRF.tHRF,freqvec, abs(spar(:,:, ch)').^2)
axis xy
listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
MeasListActplotLst = zeros(size(A,2),1);
MeasListActplotLst(PMI{1}.plotLst,:)=1;
if min(PMI{1}.plotLst-size(A,2)/2) >0
    MeasListActplotLst(PMI{1}.plotLst-size(A,2)/2,:)=1;
end
MeasListActplotLst((size(A,2)/2+1):end,1)=0; %forcer la fin a zeros

listgood = find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));



t = PMI{currentsub}.data(cf).HRF.tHRF;
spar = cat(4,spar(:,:,1:end/2),spar(:,:,end/2+1:end));
spartmp = spar(:,:,listgood,:);
PMI{1}.tmpPARAFACwavelet.checksumd = sum(intensnorm(:)); 
PMI{1}.tmpPARAFACwavelet.selected =1     
PMI{currentsub}.tmpPARAFACwavelet.spar = spartmp;                    %Data initial
PMI{currentsub}.tmpPARAFACwavelet.listgood = listgood ; 
PMI{currentsub}.tmpPARAFACwavelet.indt = [tstart(end),tstop(end)];%Time indice
PMI{currentsub}.tmpPARAFACwavelet.freqvec = freqvec;
set(guiHOMER,'UserData',PMI);
btn_Parafacwavelet_Callback(hObject, eventdata, handles)

nbcomp = str2num(get(handles.edit_nbComponent,'string'));
label = [];
for i=1:nbcomp
    label = [label,{['F',num2str(i)]}];
end
set(handles.popupmenu_Factor,'string',label)


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_ParafacIO wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_ParafacwaveletIO_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in popup_parafactID.
function popup_parafactID_Callback(hObject, eventdata, handles)
% hObject    handle to popup_parafactID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_parafactID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_parafactID


% --- Executes during object creation, after setting all properties.
function popup_parafactID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_parafactID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Parafacwavelet.
function btn_Parafacwavelet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Parafacwavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
spar = PMI{currentsub}.tmpPARAFACwavelet.spar;
Nc = str2num(get(handles.edit_nbComponent,'string'));
opt(1) = 1e-6; 
opt(2) = 0;    %initialisation 10 all methods
opt(3) = 0;     %plotting
opt(5) = 0;     %how often to show fit.
opt(6) = 0;     %Max num iterations
const = [0 0 0]; % constraints 0- nothing; 1-orthog; 2- nonneg; 3- unimod
Oldload{2} = []; % rand(21,3);

if get(handles.radio_equallambda,'value'); 
    fixMode = [0 0 1];
else
    fixMode = [0 0 0];
end
weights = []; %Mean for example put one at good time point and zero at noisy one
% figure 
% plot(spar(:,:,1))
[Factors,it,err,corcondia] = Parafac(abs(spar(:,:,:,:)),Nc,opt,const,Oldload,fixMode,weights);
PMI{currentsub}.tmpPARAFACwavelet.Factors = Factors;

set(guiHOMER,'UserData',PMI);
guidata(hObject, handles);

set(handles.text_ConcordiaVal,'string', num2str(corcondia))
set(handles.popupmenu_Factor,'value',1)
resetview(handles,0)

function edit_nbComponent_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nbComponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nbComponent as text
%        str2double(get(hObject,'String')) returns contents of edit_nbComponent as a double
nbcomp = str2num(get(handles.edit_nbComponent,'string'));
label = [];
for i=1:nbcomp
    label = [label,{['F',num2str(i)]}];
end

set(handles.popupmenu_Factor,'string',label)
 %for the topo
btn_Parafacwavelet_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_nbComponent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nbComponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function resetview(handles,newfigure)
%All normal 

axes(handles.axes_timecourse);
cla
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');


    FacTime = PMI{currentsub}.tmpPARAFACwavelet.Factors{1};
    FacWavelengt = PMI{currentsub}.tmpPARAFACwavelet.Factors{2};
    FacSpatial = PMI{currentsub}.tmpPARAFACwavelet.Factors{3};
    FacHemo = PMI{currentsub}.tmpPARAFACwavelet.Factors{4};
    Nc = str2num(get(handles.edit_nbComponent,'string'));
    colorlist = lines;
%     [0 0 1,
%         1 0 0,
%         0 1 0,
%         0 0 0,
%         1 0.5 0];
    cf = PMI{currentsub}.currentFile;
    if newfigure==0
        axes(handles.axes_timecourse);
        cla
        hold on
    else
        figure
        subplot(2,2,1)
        hold on
    end
    t = PMI{currentsub}.data(cf).HRF.tHRF(PMI{currentsub}.tmpPARAFACwavelet.indt(1):PMI{currentsub}.tmpPARAFACwavelet.indt(2));
    freqvec =  PMI{1}.tmpPARAFACwavelet.freqvec;
    for i=1:Nc
         if PMI{1}.tmpPARAFACwavelet.selected==i
            plot(t, FacTime(:,i),'color',colorlist(i,:),'linewidth',4 )
         else
            plot(t, FacTime(:,i),'color',colorlist(i,:) )
         end
    end

    if newfigure==0
        axes(handles.axes_Topo);
        cla
        hold on
    else
        subplot(2,2,2)
        hold on
    end
    for i=1:Nc
        if PMI{1}.tmpPARAFACwavelet.selected==i
            plot(FacSpatial(:,i),'color',colorlist(i,:),'linewidth',4 )
        else
            plot(FacSpatial(:,i),'color',colorlist(i,:) )
        end
    end
    if newfigure==0
        axes(handles.axes_hemo);
        cla
        hold on
     else
        subplot(2,2,3)
        hold on
    end
    for i=1:Nc
        if PMI{1}.tmpPARAFACwavelet.selected==i
            plot(FacHemo(:,i),'color',colorlist(i,:),'linewidth',4 )
        else
            plot(FacHemo(:,i),'color',colorlist(i,:) )        
        end
    end
     if newfigure==0
        axes(handles.axes_Frequence);
        cla
        hold on
     else
        subplot(2,2,4)
        hold on
    end
    for i=1:Nc
        if PMI{1}.tmpPARAFACwavelet.selected==i
            plot(freqvec,FacWavelengt(:,i),'color',colorlist(i,:),'linewidth',4 )
        else
            plot(freqvec,FacWavelengt(:,i),'color',colorlist(i,:) )        
        end
    end
    


% --- Executes on selection change in popupmenu_Factor.
function popupmenu_Factor_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Factor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Factor
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
PMI{1}.tmpPARAFACwavelet.selected =get(handles.popupmenu_Factor,'value')
set(guiHOMER,'UserData',PMI);
resetview(handles,0)
% guiHelmet = getappdata(0,'guiHelmet');
% if ~isempty(guiHelmet)  
%         Helmethandles = guihandles(guiHelmet);
%         fhresetview = getappdata(guiHelmet,'fhresetview');
%         fhresetview(Helmethandles);
% end
% --- Executes during object creation, after setting all properties.
function popupmenu_Factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function CreateNewFigure_Callback(hObject, eventdata, handles)
% hObject    handle to CreateNewFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_3dhelmet.
function btn_3dhelmet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_3dhelmet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

option.conc = 1;
option.SPMgui = 1;
IO_HelmetMTG_Display(handles,option)


% --- Executes on button press in radio_equallambda.
function radio_equallambda_Callback(hObject, eventdata, handles)
% hObject    handle to radio_equallambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_equallambda


% --------------------------------------------------------------------
function context_Newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to context_Newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetview(handles,1)


% --- Executes on button press in btn_Parafacwavelet.
function btn_ParafacWavelet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Parafacwavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1
btn_PARAFACWAVELET_Callback
% --- Executes during object creation, after setting all properties.
function text_ConcordiaVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_ConcordiaVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
