function varargout = EnterANewComment(varargin)
% ENTERANEWCOMMENT M-file for EnterANewComment.fig
%      ENTERANEWCOMMENT, by itself, creates a new ENTERANEWCOMMENT or raises the existing
%      singleton*.
%
%      H = ENTERANEWCOMMENT returns the handle to a new ENTERANEWCOMMENT or the handle to
%      the existing singleton*.
%
%      ENTERANEWCOMMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENTERANEWCOMMENT.M with the given input arguments.
%
%      ENTERANEWCOMMENT('Property','Value',...) creates a new ENTERANEWCOMMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Zonename_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EnterANewComment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EnterANewComment

% Last Modified by GUIDE v2.5 23-Jul-2024 09:41:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EnterANewComment_OpeningFcn, ...
                   'gui_OutputFcn',  @EnterANewComment_OutputFcn, ...
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


% --- Executes just before EnterANewComment is made visible.
function EnterANewComment_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EnterANewComment (see VARARGIN)

% Choose default command line output for EnterANewComment
handles.output = hObject;
if ~isempty(varargin)
    set(handles.edit1,'string',varargin{1})
end
guidata(hObject, handles);
uiwait(handles.figure1);
% Update handles structure


% UIWAIT makes EnterANewComment wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EnterANewComment_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


varargout{1} = handles.output ;
varargout{2} = handles.output2;
close(handles.figure1)


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(handles.edit1,'String');
handles.output2 = get(handles.edit_duration,'String');
% Update handles structure
guidata(handles.figure1, handles);
% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);




function edit_duration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_duration as text
%        str2double(get(hObject,'String')) returns contents of edit_duration as a double


% --- Executes during object creation, after setting all properties.
function edit_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
