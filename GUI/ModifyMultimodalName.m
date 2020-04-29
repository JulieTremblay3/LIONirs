function varargout = ModifyMultimodalName(varargin)
% MODIFYMULTIMODALNAME M-file for ModifyMultimodalName.fig
%      MODIFYMULTIMODALNAME, by itself, creates a new MODIFYMULTIMODALNAME or raises the existing
%      singleton*.
%
%      H = MODIFYMULTIMODALNAME returns the handle to a new MODIFYMULTIMODALNAME or the handle to
%      the existing singleton*.
%
%      MODIFYMULTIMODALNAME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODIFYMULTIMODALNAME.M with the given input arguments.
%
%      MODIFYMULTIMODALNAME('Property','Value',...) creates a new MODIFYMULTIMODALNAME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModifyMultimodalName_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModifyMultimodalName_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ModifyMultimodalName

% Last Modified by GUIDE v2.5 17-Jan-2020 15:46:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModifyMultimodalName_OpeningFcn, ...
                   'gui_OutputFcn',  @ModifyMultimodalName_OutputFcn, ...
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


% --- Executes just before ModifyMultimodalName is made visible.
function ModifyMultimodalName_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModifyMultimodalName (see VARARGIN)

% Choose default command line output for ModifyMultimodalName
handles.output = hObject;
NIRSpath = varargin{1};
load(NIRSpath {1});
if isfield(NIRS.Dt,'EEG')
    set(handles.listbox_EEGLIST,'string', NIRS.Dt.EEG.pp(end).p )
end
if isfield(NIRS.Dt,'AUX')
    set(handles.listbox_AUXLIST,'string', NIRS.Dt.AUX.pp(end).p )
end
if isfield(NIRS.Dt,'Video')
    set(handles.listbox_Video,'string', NIRS.Dt.Video.pp(end).p )
end
if isfield(NIRS.Dt,'Audio')
    set(handles.listbox_Audio,'string',NIRS.Dt.Audio.pp(end).p )
end
if isfield(NIRS.Cf.H,'prj')
    set(handles.listbox_MTG,'string',NIRS.Cf.H.prj )
end  
   
handles.NIRSpath = NIRSpath; 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ModifyMultimodalName wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ModifyMultimodalName_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_EEGLIST.
function listbox_EEGLIST_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_EEGLIST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_EEGLIST contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_EEGLIST


% --- Executes during object creation, after setting all properties.
function listbox_EEGLIST_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_EEGLIST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_modify_EEG.
function btn_modify_EEG_Callback(hObject, eventdata, handles)
% hObject    handle to btn_modify_EEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = get(handles.listbox_EEGLIST,'string' );

[file,path] = uigetfile(temp{1});
load(handles.NIRSpath{1});
NIRS.Dt.EEG.pp(end).p{1} = fullfile(path,file);
save(handles.NIRSpath{1},'NIRS');
set(handles.listbox_EEGLIST,'string',  NIRS.Dt.EEG.pp(end).p);

% --- Executes on selection change in listbox_AUXLIST.
function listbox_AUXLIST_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_AUXLIST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_AUXLIST contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_AUXLIST


% --- Executes during object creation, after setting all properties.
function listbox_AUXLIST_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_AUXLIST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Modify_AUX.
function btn_Modify_AUX_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Modify_AUX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = get(handles.listbox_AUXLIST,'string' )
[file,path] = uigetfile(temp{1})
load(handles.NIRSpath{1});
NIRS.Dt.AUX.pp(end).p{1} = fullfile(path,file)
save(handles.NIRSpath{1},'NIRS')
set(handles.listbox_AUXLIST,'string',  NIRS.Dt.AUX.pp(end).p)

% --- Executes on selection change in listbox_Video.
function listbox_Video_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_Video contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_Video


% --- Executes during object creation, after setting all properties.
function listbox_Video_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_Video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Modify_Video.
function btn_Modify_Video_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Modify_Video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = get(handles.listbox_Video,'string' )
[file,path] = uigetfile(temp{1})
load(handles.NIRSpath{1});
NIRS.Dt.Video.pp(end).p{1} = fullfile(path,file)
save(handles.NIRSpath{1},'NIRS')
set(handles.listbox_Video,'string',  NIRS.Dt.Video.pp(end).p)


% --- Executes on selection change in listbox_Audio.
function listbox_Audio_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_Audio contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_Audio


% --- Executes during object creation, after setting all properties.
function listbox_Audio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_Audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = get(handles.listbox_Audio,'string' )
[file,path] = uigetfile(temp{1})
load(handles.NIRSpath{1});
NIRS.Dt.Audio.pp(end).p{1} = fullfile(path,file)
save(handles.NIRSpath{1},'NIRS')
set(handles.listbox_Audio,'string',  NIRS.Dt.Audio.pp(end).p)


% --- Executes on selection change in listbox_MTG.
function listbox_MTG_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_MTG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_MTG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_MTG


% --- Executes during object creation, after setting all properties.
function listbox_MTG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_MTG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
