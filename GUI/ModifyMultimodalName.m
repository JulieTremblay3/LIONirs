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

% Last Modified by GUIDE v2.5 13-Jan-2022 14:31:01

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
filepath = [];
for i = 1:numel(NIRS.Dt.fir.pp(end).p)
   filepath = [filepath NIRS.Dt.fir.pp(end).p(i)];
end

set(handles.popupmenu_NIRS, 'String', filepath); %Set the file pop-up menu
set(handles.popupmenu_NIRS, 'Value', 1);

if isfield(NIRS.Dt,'EEG')
    set(handles.listbox_EEGLIST,'string', NIRS.Dt.EEG.pp(end).p{1} );
    try
        set(handles.edit_offsetEEG,'string', num2str(NIRS.Dt.EEG.pp(end).sync_timesec{1}) );
    catch
        set(handles.edit_offsetEEG,'string', ' ');
    end
end
if isfield(NIRS.Dt,'AUX')
    AUXname = [];
    for i=1:numel(NIRS.Dt.AUX) 
        AUXname= [AUXname;  NIRS.Dt.AUX(i).pp(end).p(1)] ;
    end
        set(handles.listbox_AUXLIST,'string', AUXname);
        set(handles.listbox_AUXLIST,'value', 1);
        set(handles.edit_offsetAUX,'string', num2str(NIRS.Dt.AUX(1).pp(end).sync_timesec{1}));
end
if isfield(NIRS.Dt,'Video')
    set(handles.listbox_Video,'string', NIRS.Dt.Video.pp(end).p{1} );
    set(handles.edit_offsetVideo,'string', num2str(NIRS.Dt.Video.pp(end).sync_timesec{1}));
end
if isfield(NIRS.Dt,'Audio')
    set(handles.listbox_Audio,'string',NIRS.Dt.Audio.pp(end).p{1} );    
    set(handles.edit_offsetAudio,'string', num2str(NIRS.Dt.Audio.pp(end).sync_timesec{1}));
end
if isfield(NIRS.Cf.H,'prj')
    set(handles.listbox_MTG,'string',NIRS.Cf.H.prj);
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

load(handles.NIRSpath{1}); 
id = get(handles.popupmenu_NIRS,'value');
NIRS.Dt.Video.pp(end).p{id} = get(handles.listbox_EEGLIST,'string');
save(handles.NIRSpath{1},'NIRS'); 
disp(['Adjust EEG file: ', get(handles.listbox_EEGLIST,'string')]);
 
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
[file,path] = uigetfile(temp);
if ~file==0
    load(handles.NIRSpath{1});
    NIRS.Dt.EEG.pp(end).p{1} = fullfile(path,file);
    save(handles.NIRSpath{1},'NIRS');
    set(handles.listbox_EEGLIST,'string',  NIRS.Dt.EEG.pp(end).p);
    disp(['Adjust EEG file: ', fullfile(path,file)]);
end

% --- Executes on selection change in listbox_AUXLIST.
function listbox_AUXLIST_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_AUXLIST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_AUXLIST contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_AUXLIST

load(handles.NIRSpath{1});
val = get(handles.listbox_AUXLIST,'value');
set(handles.edit_offsetAUX,'string', num2str(NIRS.Dt.AUX(val).pp(end).sync_timesec{1}));
 
% NIRS.Dt.AUX.pp(end).p{1} = get(handles.listbox_AUXLIST,'string');
% save(handles.NIRSpath{1},'NIRS'); 
% disp(['Adjust AUX file: ', get(handles.listbox_AUXLIST,'string')]);

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
temp = get(handles.listbox_AUXLIST,'string' );
value = get(handles.listbox_AUXLIST,'value');
[file,path] = uigetfile(temp{value});
if ~file==0
load(handles.NIRSpath{1});
id = get(handles.popupmenu_NIRS,'value');
NIRS.Dt.AUX(value).pp(end).p{id} = fullfile(path,file);
save(handles.NIRSpath{1},'NIRS');
 AUXname = [];
 for i=1:numel(NIRS.Dt.AUX) 
      AUXname= [AUXname;  NIRS.Dt.AUX(i).pp(end).p(1)] ;
 end
set(handles.listbox_AUXLIST,'string', AUXname);

disp(['Adjust AUX file: ', fullfile(path,file)])
end
% --- Executes on selection change in listbox_Video.
function listbox_Video_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_Video contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_Video
load(handles.NIRSpath{1}); 
id = get(handles.popupmenu_NIRS,'value');
NIRS.Dt.Video.pp(end).p{id} = get(handles.listbox_Video,'string');
save(handles.NIRSpath{1},'NIRS'); 
disp(['Adjust Video file: ', get(handles.listbox_Video,'string')]);

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
temp = get(handles.listbox_Video,'string' );
[file,path] = uigetfile(temp);
if ~file==0
load(handles.NIRSpath{1});
id = get(handles.popupmenu_NIRS,'value');
NIRS.Dt.Video.pp(end).p{id} = fullfile(path,file);
save(handles.NIRSpath{1},'NIRS');
set(handles.listbox_Video,'string',  NIRS.Dt.Video.pp(end).p);
disp(['Adjust Video file: ' ,fullfile(path,file) ]);
end

% --- Executes on selection change in listbox_Audio.
function listbox_Audio_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_Audio contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_Audio
load(handles.NIRSpath{1}); 
id = get(handles.popupmenu_NIRS,'value');
NIRS.Dt.Audio.pp(end).p{id} = get(handles.listbox_Audio,'string');
save(handles.NIRSpath{1},'NIRS'); 
disp(['Adjust Audio file: ', get(handles.listbox_Audio,'string')]);

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
temp = get(handles.listbox_Audio,'string' );
[file,path] = uigetfile(temp);
if  ~file==0
    load(handles.NIRSpath{1});
    id = get(handles.popupmenu_NIRS,'value');
    NIRS.Dt.Audio.pp(end).p{id} = fullfile(path,file);
    save(handles.NIRSpath{1},'NIRS');
    set(handles.listbox_Audio,'string',  NIRS.Dt.Audio.pp(end).p);
    disp(['Adjust Audio file: ' fullfile(path,file)])
end

% --- Executes on selection change in listbox_MTG.
function listbox_MTG_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_MTG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_MTG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_MTG


temp = get(handles.listbox_MTG,'string' );
load(handles.NIRSpath{1});
NIRS.Cf.H.prj = temp;
save(handles.NIRSpath{1},'NIRS');
disp(['Adjust MTG file: ', temp]);



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
temp = get(handles.listbox_MTG,'string' );
[file,path] = uigetfile(temp);
if ~file==0
    load(handles.NIRSpath{1});
    NIRS.Cf.H.prj = fullfile(path,file);
    save(handles.NIRSpath{1},'NIRS');
    set(handles.listbox_MTG,'string',  fullfile(path,file));
    disp(['Adjust MTG file: ', fullfile(path,file)]);
end


function edit_offsetVideo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_offsetVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_offsetVideo as text
%        str2double(get(hObject,'String')) returns contents of edit_offsetVideo as a double
load(handles.NIRSpath{1});
id = get(handles.popupmenu_NIRS,'value');
NIRS.Dt.Video.pp(end).sync_timesec{id} = str2num( get( handles.edit_offsetVideo,'string'));
save(handles.NIRSpath{1},'NIRS'); 
disp(['Adjust offset Video: ', get( handles.edit_offsetVideo,'string'),' seconds']);

% --- Executes during object creation, after setting all properties.
function edit_offsetVideo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_offsetVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_offsetAudio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_offsetAudio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_offsetAudio as text
%        str2double(get(hObject,'String')) returns contents of edit_offsetAudio as a double
load(handles.NIRSpath{1});
id = get(handles.popupmenu_NIRS,'value');
NIRS.Dt.Audio.pp(end).sync_timesec{id} = str2num( get( handles.edit_offsetAudio,'string'));
save(handles.NIRSpath{1},'NIRS'); 
disp(['Adjust offset Audio: ', get( handles.edit_offsetAudio,'string'),' seconds']);

% --- Executes during object creation, after setting all properties.
function edit_offsetAudio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_offsetAudio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_offsetAUX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_offsetAUX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_offsetAUX as text
%        str2double(get(hObject,'String')) returns contents of edit_offsetAUX as a double
load(handles.NIRSpath{1});
id = get(handles.popupmenu_NIRS,'value');
val = get(handles.listbox_AUXLIST,'value');
NIRS.Dt.AUX(val).pp(end).sync_timesec{id} = str2num( get( handles.edit_offsetAUX,'string'));
save(handles.NIRSpath{1},'NIRS');
disp(['Adjust offset AUX: ', get( handles.edit_offsetAUX,'string'),' seconds']);

% --- Executes during object creation, after setting all properties.
function edit_offsetAUX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_offsetAUX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_offsetEEG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_offsetEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_offsetEEG as text
%        str2double(get(hObject,'String')) returns contents of edit_offsetEEG as a double

load(handles.NIRSpath{1});
id = get(handles.popupmenu_NIRS,'value');
NIRS.Dt.EEG.pp(end).sync_timesec{id} = str2num(get( handles.edit_offsetEEG,'string'));
save(handles.NIRSpath{1},'NIRS'); 
disp(['Adjust offset EEG: ', get( handles.edit_offsetEEG,'string'),' seconds']); 

% --- Executes during object creation, after setting all properties.
function edit_offsetEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_offsetEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_NIRS.
function popupmenu_NIRS_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_NIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_NIRS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_NIRS

id = get(handles.popupmenu_NIRS,'value');
load(handles.NIRSpath{1});
set(handles.listbox_EEGLIST,'string', NIRS.Dt.EEG.pp(end).p{id} );
set(handles.listbox_AUXLIST,'string', NIRS.Dt.AUX.pp(end).p{id} );
set(handles.listbox_Video,'string', NIRS.Dt.Video.pp(end).p{id} );
set(handles.listbox_Audio,'string',NIRS.Dt.Audio.pp(end).p{id} );   

set(handles.edit_offsetEEG,'string',num2str(NIRS.Dt.EEG.pp(end).sync_timesec{id}));
set(handles.edit_offsetAUX,'string',num2str(NIRS.Dt.AUX.pp(end).sync_timesec{id}));
set(handles.edit_offsetVideo,'string',num2str(NIRS.Dt.Video.pp(end).sync_timesec{id}));
set(handles.edit_offsetAudio,'string',num2str(NIRS.Dt.Audio.pp(end).sync_timesec{id}));


% --- Executes during object creation, after setting all properties.
function popupmenu_NIRS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_NIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
