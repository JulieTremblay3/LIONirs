function varargout = ExportFigureOption(varargin)
% EXPORTFIGUREOPTION M-file for ExportFigureOption.fig
%      EXPORTFIGUREOPTION, by itself, creates a new EXPORTFIGUREOPTION or raises the existing
%      singleton*.
%
%      H = EXPORTFIGUREOPTION returns the handle to a new EXPORTFIGUREOPTION or the handle to
%      the existing singleton*.
%
%      EXPORTFIGUREOPTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXPORTFIGUREOPTION.M with the given input arguments.
%
%      EXPORTFIGUREOPTION('Property','Value',...) creates a new EXPORTFIGUREOPTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ExportFigureOption_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ExportFigureOption_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ExportFigureOption

% Last Modified by GUIDE v2.5 10-Sep-2010 11:02:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ExportFigureOption_OpeningFcn, ...
                   'gui_OutputFcn',  @ExportFigureOption_OutputFcn, ...
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


% --- Executes just before ExportFigureOption is made visible.
function ExportFigureOption_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ExportFigureOption (see VARARGIN)

% Choose default command line output for ExportFigureOption
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ExportFigureOption wait for user response (see UIRESUME)
% uiwait(handles.ExportFigureOptions);


% --- Outputs from this function are returned to the command line.
function varargout = ExportFigureOption_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radio_title.
function radio_title_Callback(hObject, eventdata, handles)
% hObject    handle to radio_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_title


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radio_colorbar.
function radio_colorbar_Callback(hObject, eventdata, handles)
% hObject    handle to radio_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_colorbar



function edit_xzoom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xzoom as text
%        str2double(get(hObject,'String')) returns contents of edit_xzoom as a double


% --- Executes during object creation, after setting all properties.
function edit_xzoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yzoom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yzoom as text
%        str2double(get(hObject,'String')) returns contents of edit_yzoom as a double


% --- Executes during object creation, after setting all properties.
function edit_yzoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_zzoom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zzoom as text
%        str2double(get(hObject,'String')) returns contents of edit_zzoom as a double


% --- Executes during object creation, after setting all properties.
function edit_zzoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xpixel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xpixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xpixel as text
%        str2double(get(hObject,'String')) returns contents of edit_xpixel as a double


% --- Executes during object creation, after setting all properties.
function edit_xpixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xpixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ypixel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ypixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ypixel as text
%        str2double(get(hObject,'String')) returns contents of edit_ypixel as a double


% --- Executes during object creation, after setting all properties.
function edit_ypixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ypixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in btn_OK.
function btn_OK_Callback(hObject, eventdata, handles)
% hObject    handle to btn_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setappdata(0,'guiExportFigureOption',handles.ExportFigureOptions)



% --- Executes on button press in btn_browsehorz.
function btn_browsehorz_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browsehorz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[path]=  uigetdir()
set(handles.edit_pathhorz,'string',path)



% --- Executes on selection change in list_horz.
function list_horz_Callback(hObject, eventdata, handles)
% hObject    handle to list_horz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_horz contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_horz


% --- Executes during object creation, after setting all properties.
function list_horz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_horz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_okhorz.
function btn_okhorz_Callback(hObject, eventdata, handles)
% hObject    handle to btn_okhorz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = [get(handles.edit_pathhorz,'string'),'\'];
subjectname = get(handles.edit_subjectname,'string');
time = str2num(get(handles.edit_nbtime,'string'));
type = get(handles.edit_type,'string');
limitex =  str2num(get(handles.edit_xlim,'string')); 
limitey = str2num(get(handles.edit_ylim,'string')); 
file=[]
for i = 1:numel(time)
    if isempty(file)
        file = {[subjectname,sprintf('%1.1f',time(i)),type,'.tif']};
    else
        file = [file,{[subjectname,sprintf('%1.1f',time(i)),type,'.tif']}]
    end
end

[fileout1,pathout]= combinetiffileH(path,file,limitex,limitey,time)
fileout = get(handles.list_vert,'string')
set(handles.edit_pathvert,'string',pathout)
if isempty(fileout)
    set(handles.list_vert,'string',{fileout1})
else
    fileout = [fileout;{fileout1}]
    set(handles.list_vert,'string',fileout)
end

% --- Executes on button press in btn_browsevert.
function btn_browsevert_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browsevert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uigetfile('.tif','MultiSelect','on')
if path==0
    set(handles.edit_pathvert,'string','')
    set(handles.list_vert,'string','')
return
end
set(handles.edit_pathvert,'string',path)
set(handles.list_vert,'string',file)


% --- Executes on selection change in list_vert.
function list_vert_Callback(hObject, eventdata, handles)
% hObject    handle to list_vert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_vert contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_vert


% --- Executes during object creation, after setting all properties.
function list_vert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_vert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_okvert.
function btn_okvert_Callback(hObject, eventdata, handles)
% hObject    handle to btn_okvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pathout = get(handles.edit_pathvert,'string');
fileall = get(handles.list_vert,'string');
cmin = str2num(get(handles.edit_cmin,'string'));
cmax = str2num(get(handles.edit_cmax,'string'));
thresh = str2num(get(handles.edit_tresh,'string'));
[pathout,fileout]= combinetiffilev(pathout,fileall,cmin,cmax,thresh)

function edit_pathhorz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathhorz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathhorz as text
%        str2double(get(hObject,'String')) returns contents of edit_pathhorz as a double


% --- Executes during object creation, after setting all properties.
function edit_pathhorz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathhorz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_nbtime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nbtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nbtime as text
%        str2double(get(hObject,'String')) returns contents of edit_nbtime as a double


% --- Executes during object creation, after setting all properties.
function edit_nbtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nbtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xlim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xlim as text
%        str2double(get(hObject,'String')) returns contents of edit_xlim as a double


% --- Executes during object creation, after setting all properties.
function edit_xlim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ylim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ylim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ylim as text
%        str2double(get(hObject,'String')) returns contents of edit_ylim as a double


% --- Executes during object creation, after setting all properties.
function edit_ylim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ylim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_subjectname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subjectname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subjectname as text
%        str2double(get(hObject,'String')) returns contents of edit_subjectname as a double


% --- Executes during object creation, after setting all properties.
function edit_subjectname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subjectname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cmin as text
%        str2double(get(hObject,'String')) returns contents of edit_cmin as a double


% --- Executes during object creation, after setting all properties.
function edit_cmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cmax as text
%        str2double(get(hObject,'String')) returns contents of edit_cmax as a double


% --- Executes during object creation, after setting all properties.
function edit_cmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tresh as text
%        str2double(get(hObject,'String')) returns contents of edit_tresh as a double


% --- Executes during object creation, after setting all properties.
function edit_tresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_type_Callback(hObject, eventdata, handles)
% hObject    handle to edit_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_type as text
%        str2double(get(hObject,'String')) returns contents of edit_type as a double


% --- Executes during object creation, after setting all properties.
function edit_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_pathvert_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathvert as text
%        str2double(get(hObject,'String')) returns contents of edit_pathvert as a double


% --- Executes during object creation, after setting all properties.
function edit_pathvert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


