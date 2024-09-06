function varargout = GUI_ViewNetwork(varargin)
% GUI_VIEWNETWORK M-file for GUI_ViewNetwork.fig
%      GUI_VIEWNETWORK, by itself, creates a new GUI_VIEWNETWORK or raises the existing
%      singleton*.
%
%      H = GUI_VIEWNETWORK returns the handle to a new GUI_VIEWNETWORK or the handle to
%      the existing singleton*.
%
%      GUI_VIEWNETWORK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_VIEWNETWORK.M with the given input arguments.
%
%      GUI_VIEWNETWORK('Property','Value',...) creates a new GUI_VIEWNETWORK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ViewNetwork_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ViewNetwork_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ViewNetwork

% Last Modified by GUIDE v2.5 14-Mar-2017 14:23:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ViewNetwork_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ViewNetwork_OutputFcn, ...
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


% --- Executes just before GUI_ViewNetwork is made visible.
function GUI_ViewNetwork_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ViewNetwork (see VARARGIN)

% Choose default command line output for GUI_ViewNetwork
handles.output = hObject;
cameratoolbar(handles.figure1);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_ViewNetwork wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ViewNetwork_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_subnetwork.
function popupmenu_subnetwork_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_subnetwork (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_subnetwork contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_subnetwork


% --- Executes during object creation, after setting all properties.
function popupmenu_subnetwork_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_subnetwork (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
