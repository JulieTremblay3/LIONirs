function varargout = CombineFile(varargin)
% COMBINEFILE M-file for CombineFile.fig
%      COMBINEFILE, by itself, creates a new COMBINEFILE or raises the existing
%      singleton*.
%
%      H = COMBINEFILE returns the handle to a new COMBINEFILE or the handle to
%      the existing singleton*.
%
%      COMBINEFILE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMBINEFILE.M with the given input arguments.
%
%      COMBINEFILE('Property','Value',...) creates a new COMBINEFILE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CombineFile_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CombineFile_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CombineFile

% Last Modified by GUIDE v2.5 17-Nov-2011 17:22:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CombineFile_OpeningFcn, ...
                   'gui_OutputFcn',  @CombineFile_OutputFcn, ...
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


% --- Executes just before CombineFile is made visible.
function CombineFile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CombineFile (see VARARGIN)

% Choose default command line output for CombineFile
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CombineFile wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CombineFile_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_Horizontal.
function btn_Horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('.tif','multiselect','on')
Total = [];
for i = 1:numel(file)
    [X] = imread([path,file{i}],'tif'); 
    Total = [Total,X];
end
imwrite(Total,[path,'all', file{1},'.tif'],'tif');


% --- Executes on button press in btn_Vertical.
function btn_Vertical_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile('.tif','multiselect','on')
Total = [];
for i = 1:numel(file)
    [X] = imread([path,file{i}],'tif'); 
    Total = [Total;X];
end
imwrite(Total,[path,'all', file{1},'.tif'],'tif');
