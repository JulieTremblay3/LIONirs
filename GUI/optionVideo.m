function varargout = optionVideo(varargin)
% OPTIONVIDEO M-file for optionVideo.fig
%      OPTIONVIDEO, by itself, creates a new OPTIONVIDEO or raises the existing
%      singleton*.
%
%      H = OPTIONVIDEO returns the handle to a new OPTIONVIDEO or the handle to
%      the existing singleton*.
%
%      OPTIONVIDEO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTIONVIDEO.M with the given input arguments.
%
%      OPTIONVIDEO('Property','Value',...) creates a new OPTIONVIDEO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before optionVideo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to optionVideo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help optionVideo

% Last Modified by GUIDE v2.5 07-Nov-2019 14:16:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @optionVideo_OpeningFcn, ...
                   'gui_OutputFcn',  @optionVideo_OutputFcn, ...
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


% --- Executes just before optionVideo is made visible.
function optionVideo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to optionVideo (see VARARGIN)

% Choose default command line output for optionVideo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes optionVideo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = optionVideo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_codecoption.
function popupmenu_codecoption_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_codecoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_codecoption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_codecoption

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
if get(handles.popupmenu_codecoption,'value')==1
    PMI{1}.videooption.codec = 'VideoReader'; 
else
    PMI{1}.videooption.codec = 'mmread'; 
end
set(guiHOMER,'UserData',PMI);

% --- Executes during object creation, after setting all properties.
function popupmenu_codecoption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_codecoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
temp = get(handles.slider1,'value');

PMI{1}.videooption.brighness = temp*255; 

set(guiHOMER,'UserData',PMI);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in radio_DisplayVideo.
function radio_DisplayVideo_Callback(hObject, eventdata, handles)
% hObject    handle to radio_DisplayVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_DisplayVideo
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
PMI{1}.videoview = get(handles.radio_DisplayVideo,'value'); 

set(guiHOMER,'UserData',PMI);
