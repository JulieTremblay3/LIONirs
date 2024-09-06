function varargout = gui_MorletTG(varargin)
% GUI_MORLETTG M-file for gui_MorletTG.fig
%      GUI_MORLETTG, by itself, creates a new GUI_MORLETTG or raises the existing
%      singleton*.
%
%      H = GUI_MORLETTG returns the handle to a new GUI_MORLETTG or the handle to
%      the existing singleton*.
%
%      GUI_MORLETTG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MORLETTG.M with the given input arguments.
%
%      GUI_MORLETTG('Property','Value',...) creates a new GUI_MORLETTG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_MorletTG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_MorletTG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_MorletTG

% Last Modified by GUIDE v2.5 16-Nov-2017 15:51:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_MorletTG_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_MorletTG_OutputFcn, ...
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


% --- Executes just before gui_MorletTG is made visible.
function gui_MorletTG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_MorletTG (see VARARGIN)

% Choose default command line output for gui_MorletTG
handles.output = hObject;

time_start = varargin{1};
time_stop = varargin{2};
if numel(varargin)>2
    EEG = varargin{3};
end
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf = PMI{1}.currentFile;

PMI{1}.tmpWavelet.EEG = EEG
PMI{1}.tmpWavelet.time_start = time_start;
PMI{1}.tmpWavelet.time_stop = time_stop;
set(guiHOMER,'UserData',PMI);
NIRSlabelall = [];
if ~isfield(PMI{1},'zone')
   msgbox('Please make zone before to facilate the data display')
end
for izone=1:numel(PMI{1}.zone.label)
      plotLst =  PMI{1}.zone.plotLst{izone}
    for ich = 1:numel(plotLst)
     
        NIRSlabelall = [NIRSlabelall,{[PMI{1}.zone.label{izone},',',num2str(plotLst(ich))]}]
    end
end
set(handles.popup_NIRS,'string',NIRSlabelall)
set(handles.popup_EEGCH,'string',EEG.ele)


%  
% %Detrent DATA segment for centrering 
% X = 1:1:size(intensnorm,1);
% Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
% Mb2 =  intensnorm(1,:)'; %offset    
% A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
% spar = intensnorm - A;
% 
% listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
% MeasListActplotLst = zeros(size(A,2),1);
% MeasListActplotLst(PMI{1}.plotLst,:)=1;
% if min(PMI{1}.plotLst-size(A,2)/2) >0
%     MeasListActplotLst(PMI{1}.plotLst-size(A,2)/2,:)=1;
% end
% MeasListActplotLst((size(A,2)/2+1):end,1)=0; %forcer la fin a zeros
% 
% listgood = find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));
% 
% 
% 
% t = PMI{currentsub}.data(cf).HRF.tHRF;
% spar = cat(3,spar(:,1:end/2),spar(:,end/2+1:end));
% spartmp = spar(:,listgood,:);
%  
% PMI{1}.tmpPARAFAC.selected =1     
% PMI{currentsub}.tmpPARAFAC.spar = spartmp;                    %Data initial
% PMI{currentsub}.tmpPARAFAC.listgood = listgood ; 
% PMI{currentsub}.tmpPARAFAC.indt = [tstart(end),tstop(end)];%Time indice
% set(guiHOMER,'UserData',PMI);
% btn_Parafac_Callback(hObject, eventdata, handles)
% 
% nbcomp = str2num(get(handles.edit_nbComponent,'string'));
% label = [];
% for i=1:nbcomp
%     label = [label,{['F',num2str(i)]}];
% end
% set(handles.popupmenu_Factor,'string',label)


% Update handles structure

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_MorletTG wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_MorletTG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popup_EEGCH.
function popup_EEGCH_Callback(hObject, eventdata, handles)
% hObject    handle to popup_EEGCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_EEGCH contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_EEGCH
resetview_EEG(handles)

% --- Executes during object creation, after setting all properties.
function popup_EEGCH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_EEGCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_NIRS.
function popup_NIRS_Callback(hObject, eventdata, handles)
% hObject    handle to popup_NIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_NIRS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_NIRS
resetview_NIRS(handles)

% --- Executes during object creation, after setting all properties.
function popup_NIRS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_NIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_freqNIRS_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freqNIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freqNIRS as text
%        str2double(get(hObject,'String')) returns contents of edit_freqNIRS as a double


% --- Executes during object creation, after setting all properties.
function edit_freqNIRS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freqNIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_widthNIRS_Callback(hObject, eventdata, handles)
% hObject    handle to edit_widthNIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_widthNIRS as text
%        str2double(get(hObject,'String')) returns contents of edit_widthNIRS as a double


% --- Executes during object creation, after setting all properties.
function edit_widthNIRS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_widthNIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_runwavelet.
function btn_runwavelet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_runwavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure w% 
% currentsub=1;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf=  1;
currentsub = 1;
d = PMI{1}.data(cf).HRF.AvgC;
time_start = PMI{1}.tmpWavelet.time_start;
time_stop = PMI{1}.tmpWavelet.time_stop;

%Look Start stop interval
tstart=find(PMI{currentsub}.data(cf).HRF.tHRF<time_start);
if isempty(tstart)
    sprintf('Please check the time start')
end
    
if isempty(tstart);tstart=1;end
tstop=find(PMI{currentsub}.data(cf).HRF.tHRF<time_stop);
indt = tstart(end):tstop(end);
dint = d(indt,:);


Fs = 1/(PMI{currentsub}.data.HRF.tHRF(2)- PMI{currentsub}.data.HRF.tHRF(1));
freqVec = str2num(get(handles.edit_freqNIRS,'string')); % Hz.  Frequencies where to compute the CWT
width = str2num(get(handles.edit_widthNIRS,'string'));
tic; PMI{1}.tmpWavelet.NIRS.TF = morletcwtEd(dint,freqVec,Fs,width); tori=toc; % Complex Wavelets Coefficients
PMI{1}.tmpWavelet.NIRS.t = 1/Fs:1/Fs:(numel(indt)*1/Fs);
PMI{1}.tmpWavelet.NIRS.width = width;
PMI{1}.tmpWavelet.NIRS.freqVec = freqVec;


FsEEG = 1/(PMI{1}.tmpWavelet.EEG.t(2)- PMI{1}.tmpWavelet.EEG.t(1));
freqVec = str2num(get(handles.edit_freqEEG,'string')); % Hz.  Frequencies where to compute the CWT
width = str2num(get(handles.edit_widthEEG,'string')); %
n = str2num(get(handles.edit_downsampleEEG,'string'));
Fs = FsEEG/n;
t = downsample(PMI{1}.tmpWavelet.EEG.t,n);
for i=1:size(PMI{1}.tmpWavelet.EEG.data,2)
    data(:,i) =decimate(PMI{1}.tmpWavelet.EEG.data(:,i),n);
end
set(guiHOMER,'UserData',PMI);
resetview_NIRS(handles)

tic; PMI{1}.tmpWavelet.EEG.TF = morletcwtEd(data,freqVec,Fs,width); tori=toc; % Complex Wavelets Coefficients
PMI{1}.tmpWavelet.EEG.t = 1/Fs:1/Fs:(numel(t)*1/Fs);
PMI{1}.tmpWavelet.EEG.freqVec = freqVec;
PMI{1}.tmpWavelet.EEG.width = width;
PMI{1}.tmpWavelet.EEG.downsample = n;
set(guiHOMER,'UserData',PMI);
resetview_EEG(handles)


 



function edit_downsampleEEG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_downsampleEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_downsampleEEG as text
%        str2double(get(hObject,'String')) returns contents of edit_downsampleEEG as a double


% --- Executes during object creation, after setting all properties.
function edit_downsampleEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_downsampleEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_freqEEG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freqEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freqEEG as text
%        str2double(get(hObject,'String')) returns contents of edit_freqEEG as a double


% --- Executes during object creation, after setting all properties.
function edit_freqEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freqEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_widthEEG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_widthEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_widthEEG as text
%        str2double(get(hObject,'String')) returns contents of edit_widthEEG as a double


% --- Executes during object creation, after setting all properties.
function edit_widthEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_widthEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function resetview_NIRS(handles)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
axes(handles.axes1);
TF = PMI{1}.tmpWavelet.NIRS.TF;
freqVec = PMI{1}.tmpWavelet.NIRS.freqVec;
t = PMI{1}.tmpWavelet.NIRS.t ;
SEEG =  abs(TF).^2;
id = get(handles.popup_NIRS,'value');
labelall = get(handles.popup_NIRS,'string');
[tok,rem]=strtok(labelall{id},',');
ch =str2num(rem(2:end));
if 1==get(handles.popup_optionDisplay,'value')
        NORM= SEEG(:,:,ch);
elseif 2==get(handles.popup_optionDisplay,'value') 
    NORM= SEEG(:,:,ch)./(ones(size(SEEG),1)*mean(SEEG(:,:,ch),1));
end
imagesc(t,freqVec,NORM');
axis xy;
% figure
% hist(NORM(:))
% new = zscore(NORM);
% figure
% hist(new(:))
% idoulier = find(new>3);



function resetview_EEG(handles)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
axes(handles.axes2)
TF = PMI{1}.tmpWavelet.EEG.TF;
freqVec = PMI{1}.tmpWavelet.EEG.freqVec
t = PMI{1}.tmpWavelet.EEG.t 
SEEG =  abs(TF).^2;
ch = get(handles.popup_EEGCH,'value')
if 1==get(handles.popup_optionDisplay,'value')
    NORM= SEEG(:,:,ch);
elseif 2==get(handles.popup_optionDisplay,'value')    
    NORM= SEEG(:,:,ch)./(ones(size(SEEG),1)*mean(SEEG(:,:,ch),1));
end
imagesc(t,freqVec, NORM')   
axis xy


% --- Executes on selection change in popup_optionDisplay.
function popup_optionDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to popup_optionDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_optionDisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_optionDisplay
resetview_EEG(handles)
resetview_NIRS(handles)

% --- Executes during object creation, after setting all properties.
function popup_optionDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_optionDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
