function varargout = optionEEG(varargin)
% OPTIONEEG M-file for optionEEG.fig
%      OPTIONEEG, by itself, creates a new OPTIONEEG or raises the existing
%      singleton*.
%
%      H = OPTIONEEG returns the handle to a new OPTIONEEG or the handle to
%      the existing singleton*.
%
%      OPTIONEEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTIONEEG.M with the given input arguments.
%
%      OPTIONEEG('Property','Value',...) creates a new OPTIONEEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before optionEEG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to optionEEG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help optionEEG

% Last Modified by GUIDE v2.5 13-Jan-2022 09:03:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @optionEEG_OpeningFcn, ...
                   'gui_OutputFcn',  @optionEEG_OutputFcn, ...
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


% --- Executes just before optionEEG is made visible.
function optionEEG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to optionEEG (see VARARGIN)

% Choose default command line output for optionEEG
handles.output = hObject;

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
if isfield(PMI{1},'EEGoption')
   if isfield(PMI{1}.EEGoption,'filter');
       if PMI{1}.EEGoption.filter
            set(handles.text_LPF ,'visible','on')
            set(handles.text_HPF,'visible','on')
            set(handles.edit_LPF,'visible','on')
            set(handles.edit_HPF,'visible','on')
       else
            set(handles.text_LPF ,'visible','off')
            set(handles.text_HPF,'visible','off')
            set(handles.edit_LPF,'visible','off')
            set(handles.edit_HPF,'visible','off')
       end
   end
end
PMI = get(guiHOMER,'UserData');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes optionEEG wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = optionEEG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_ELE.
function listbox_ELE_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_ELE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_ELE contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_ELE


% --- Executes during object creation, after setting all properties.
function listbox_ELE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_ELE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LPF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LPF as text
%        str2double(get(hObject,'String')) returns contents of edit_LPF as a double
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    PMI = get(guiHOMER,'UserData');
    PMI{1}.EEGoption.filterLPF = str2num(get(handles.edit_LPF,'string'));
    set(guiHOMER,'UserData',PMI);

% --- Executes during object creation, after setting all properties.
function edit_LPF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_HPF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_HPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_HPF as text
%        str2double(get(hObject,'String')) returns contents of edit_HPF as a double
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    PMI = get(guiHOMER,'UserData');
    PMI{1}.EEGoption.filterHPF = str2num(get(handles.edit_HPF,'string'));
    set(guiHOMER,'UserData',PMI);

% --- Executes during object creation, after setting all properties.
function edit_HPF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_HPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_BDPF.
function radiobutton_BDPF_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_BDPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_BDPF
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
if get(handles.radiobutton_BDPF,'value')
    set(handles.text_LPF ,'visible','on')
    set(handles.text_HPF,'visible','on')
    set(handles.edit_LPF,'visible','on')
    set(handles.edit_HPF,'visible','on')
    PMI{1}.EEGoption.filter = 1;
    PMI{1}.EEGoption.filterLPF = str2num(get(handles.edit_LPF,'string'));
    PMI{1}.EEGoption.filterHPF = str2num(get(handles.edit_HPF,'string'));
else
    set(handles.text_LPF ,'visible','off')
    set(handles.text_HPF,'visible','off')
    set(handles.edit_LPF,'visible','off')
    set(handles.edit_HPF,'visible','off')
    PMI{1}.EEGoption.filter = 0;
end
set(guiHOMER,'UserData',PMI);
