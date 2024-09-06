function varargout = NIRS_GUI(varargin)
% NIRS_GUI M-file for NIRS_GUI.fig
%      NIRS_GUI, by itself, creates a new NIRS_GUI or raises the existing
%      singleton*.
%
%      H = NIRS_GUI returns the handle to a new NIRS_GUI or the handle to
%      the existing singleton*.
%
%      NIRS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NIRS_GUI.M with the given input arguments.
%
%      NIRS_GUI('Property','Value',...) creates a new NIRS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NIRS_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NIRS_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NIRS_GUI

% Last Modified by GUIDE v2.5 03-Apr-2020 17:56:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NIRS_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NIRS_GUI_OutputFcn, ...
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


% --- Executes just before NIRS_GUI is made visible.
function NIRS_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NIRS_GUI (see VARARGIN)

% Choose default command line output for NIRS_GUI

handles.output = hObject;

try
    if numel(varargin)>=1
    handles.NIRSpath = varargin{1};
    [pathstr, ~, ~] = fileparts(handles.NIRSpath.name);
    handles.PROJECTpath = [pathstr, filesep];
    handles = load_NIRSmat(handles);
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NIRS_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% set(handles.figure1,'WindowStyle','modal');


% --- Outputs from this function are returned to the command line.
function varargout = NIRS_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.NIRSpath;
%uiwait(hObject);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    [FileName,PathName] = uigetfile('*.mat','Select the NIRS.mat file',handles.PROJECTpath);
catch
    [FileName,PathName] = uigetfile('*.mat','Select the NIRS.mat file');
end
handles.NIRSpath = [PathName,FileName];
handles.PROJECTpath = PathName;

handles = load_NIRSmat(handles);

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = load_NIRSmat(handles)
%
%
set(handles.edit1, 'String', handles.NIRSpath.name);

load(handles.NIRSpath.name);
handles.NIRS = NIRS;

%Find modules names
module_list = [];
handles.module_list = [];
for i = 1:length(NIRS.Dt.fir.pp)
    handles.module_list{i} = NIRS.Dt.fir.pp(i).pre;
    handles.module_path{i} = NIRS.Dt.fir.pp(i).p{1};
end
set(handles.listbox1, 'String', handles.module_list);
set(handles.listbox1, 'Value', 1);
set(handles.popupmenu1, 'String', handles.module_path{i});
set(handles.popupmenu1,'value',1);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.NIRSpath = get(hObject,'String');
[pathstr, ~, ~, ~] = fileparts(handles.NIRSpath);
handles.PROJECTpath = [pathstr, filesep];

load(handles.NIRSpath);
handles.NIRS = NIRS;

handles.current_field = fieldnames(NIRS);
set(handles.listbox1, 'String', handles.current_field);

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
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

NIRS = handles.NIRS;

directory_name = uigetdir(handles.PROJECTpath,'Select new path for project');
if directory_name==0
    return
end
directory_name = [directory_name, filesep];
old_path = NIRS.Dt.s.p;
NIRS.Dt.s.p = directory_name;

%Preprocessing steps
lst = length(NIRS.Dt.fir.pp);
for i = 1:lst
    nb_files = length(NIRS.Dt.fir.pp(i).p);
    for j = 1:nb_files
        [pathstr,nir_name,ext] = fileparts(NIRS.Dt.fir.pp(i).p{j});
        if ~strcmp(pathstr, old_path);
            subfolder = [pathstr(length(old_path)+1:end), filesep];
            NIRS.Dt.fir.pp(i).p{j} = [directory_name, subfolder, nir_name, ext];
        else
            NIRS.Dt.fir.pp(i).p{j} = [directory_name, nir_name, ext];
        end
    end
end

%Other fields (actually unused?)
fields = cellstr(strvcat('NIRS.Dt.fir.stax.p','NIRS.Dt.fir.rons','NIRS.Dt.fir.topo.p','NIRS.Dt.fir.tomo.p','NIRS.Dt.ana.p',...
    'NIRS.Dt.aux.log.p','NIRS.Dt.pro.p'));
for i = 1:length(fields)
    try
        actual_field = eval(fields{i});
        [pathstr,nir_name,ext] = fileparts(actual_field);
        if ~strcmp(pathstr, old_path);
            subfolder = [pathstr(length(old_path)+1:end), filesep];
            actual_field = [directory_name, subfolder, nir_name, ext];
        else
            actual_field = [directory_name, nir_name, ext];
        end
        eval([fields{i} ' = actual_field']);
    end
end

% Other fields with subfields...???
% NIRS.Dt.fir.vR()
% NIRS.Dt.fir.fR()
% NIRS.Dt.fir.ht{f,1}

% Option for deleting the old path?
save([handles.PROJECTpath,'NIRS.mat'],'NIRS','-mat')
msgbox('NIRS.mat path fields were correctly updated.','NIRS.mat','help');


handles.NIRS = NIRS;

guidata(hObject, handles);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NIRS = handles.NIRS;

% contents = cellstr(get(hObject,'String'));
selected_field = get(hObject,'Value');

%Find modules path
nb_file = length(NIRS.Dt.fir.pp(selected_field).p);
handles.filepath = [];
for i = 1:nb_file
    handles.filepath = [handles.filepath NIRS.Dt.fir.pp(selected_field).p(i)];
end

set(handles.popupmenu1, 'String', handles.filepath);

% set(handles.edit2, 'String', handles.fieldpath);
% set(hObject, 'Value', 1);

guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_lastoperation.
function btn_lastoperation_Callback(hObject, eventdata, handles)
% hObject    handle to btn_lastoperation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
NIRS = handles.NIRS;

%Delete subsequent fields
selected_field = get(handles.listbox1,'Value');
if selected_field < length(NIRS.Dt.fir.pp)
    %deleted old referenced file
    modsdel = selected_field+1;
    modedel = numel(NIRS.Dt.fir.pp);
    for imod = modsdel: modedel
        for ifile = 1:numel(NIRS.Dt.fir.pp(imod).p)
           [path, fname, extension] =fileparts(NIRS.Dt.fir.pp(imod).p{ifile})
           try; delete(fullfile(path, [fname,extension]));
           catch; disp(['error can''t delete',fullfile(path, [fname,extension])]);  end
           try;delete(fullfile(path, [fname,'.vmrk']));
           catch; disp(['error can''t delete',fullfile(path, [fname,'.vmrk'])]);  end
           try;delete(fullfile(path, [fname,'.vhdr']));
           catch; disp(['error can''t delete',fullfile(path, [fname,'.vhdr'])]);  end
        end
    end
    NIRS.Dt.fir.pp(selected_field+1:end) = [];

end

%Copy .nir files with "mod"suffix
if get(handles.radio_saveoldNIRS,'value')
    suffix = 'mod';
    for i = 1:length(NIRS.Dt.fir.pp(selected_field).p)
        filepath = NIRS.Dt.fir.pp(selected_field).p{i};
        [pathstr, name, ext, ~] = fileparts(filepath);
        NIRS.Dt.fir.pp(selected_field).p{i} = [pathstr, filesep, name, suffix, ext]; %Change files' names
        %Make new files
        copyfile(filepath,NIRS.Dt.fir.pp(selected_field).p{i});
         %write new .vmrk file
         [dir1,fil1,~] = fileparts(filepath);
         [dir2,fil2,~] = fileparts(NIRS.Dt.fir.pp(selected_field).p{i});
         infilewmrk = fullfile(dir1,[fil1 '.vmrk']);
         outfilevmrk = fullfile(dir2,[fil2 '.vmrk']);
         copyfile(infilewmrk,outfilevmrk);
    end
        saving_path = [handles.PROJECTpath, date,'_NIRS.mat'];
        copyfile(handles.NIRSpath,saving_path); %Copy and rename old NIRS.mat
        delete(handles.NIRSpath);
end

[FileName,PathName] = uiputfile('.mat','Save modified NIRS.mat',handles.NIRSpath.name);
save([PathName FileName], 'NIRS'); %Save the new NIRS.mat (with the name "NIRS.mat")

%Load modified NIRS 
handles.PROJECTpath = PathName;
handles.NIRSpath.name = [PathName FileName];
[handles] = load_NIRSmat(handles);

guidata(hObject,handles);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1); 
uiresume(gcf);


% --- Executes on button press in radio_saveoldNIRS.
function radio_saveoldNIRS_Callback(hObject, eventdata, handles)
% hObject    handle to radio_saveoldNIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_saveoldNIRS





function edit_stringfind_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stringfind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stringfind as text
%        str2double(get(hObject,'String')) returns contents of edit_stringfind as a double


% --- Executes during object creation, after setting all properties.
function edit_stringfind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stringfind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_replace_Callback(hObject, eventdata, handles)
% hObject    handle to edit_replace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_replace as text
%        str2double(get(hObject,'String')) returns contents of edit_replace as a double


% --- Executes during object creation, after setting all properties.
function edit_replace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_replace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_replace.
function btn_replace_Callback(hObject, eventdata, handles)
% hObject    handle to btn_replace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
NIRS = handles.NIRS;
pattern = get(handles.edit_stringfind,'string')
replacepattern = get(handles.edit_replace,'string') 
findonce = 0;
for imodule=1:numel(NIRS.Dt.fir.pp)
    for ifile = 1:numel(NIRS.Dt.fir.pp(imodule).p)
        str = NIRS.Dt.fir.pp(imodule).p{ifile}
        k = strfind(str, pattern);
        if numel(k)==1           
            strnew = [str(1:k-1), replacepattern, str(k+numel(pattern):end)];
            findonce =  findonce+1;
            NIRS.Dt.fir.pp(imodule).p{ifile} = strnew;
        elseif numel(k)>2
            ['Ensure the string ', pattern,' could be find once in the file'];
        end
    end
end
if findonce==0
   msgbox(['Ensure the string ', pattern,' is present in the file']);
   return
end
handles.NIRS = NIRS;
save([handles.PROJECTpath,'NIRS.mat'],'NIRS','-mat');
handles = load_NIRSmat(handles);
guidata(hObject,handles);
msgbox(['NIRS.mat path fields were correctly updated. ',num2str( findonce),' name are changed.'],'NIRS.mat','help');


% --- Executes on button press in btn_DeleteFile.
function btn_DeleteFile_Callback(hObject, eventdata, handles)
% hObject    handle to btn_DeleteFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imodule = get(handles.listbox1,'value');
NIRS =  handles.NIRS;
warning('off')
for i=1:numel(NIRS.Dt.fir.pp(imodule).p)
    [pathdir,file,ext]=fileparts(NIRS.Dt.fir.pp(imodule).p{i});
    try     
        delete(fullfile(pathdir,[file,ext]))
     catch
        disp([fullfile(pathdir,[file,ext]), ' could not be delete'])
    end
    try     
        delete(fullfile(pathdir,[file,'.vhdr']))
    catch
        disp([fullfile(pathdir,[file,'.vhdr']), ' could not be delete'])
    end
    try
        delete(fullfile(pathdir,[file,'.vmrk']))
    catch
        disp([fullfile(pathdir,[file,'.vmrk']), ' could not be delete'])
    end   
end
 disp('DONE')
