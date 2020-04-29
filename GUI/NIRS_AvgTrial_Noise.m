function varargout = NIRS_AvgTrial_Noise(varargin)
% NIRS_AVGTRIAL_NOISE M-file for NIRS_AvgTrial_Noise.fig
%      NIRS_AVGTRIAL_NOISE, by itself, creates a new NIRS_AVGTRIAL_NOISE or raises the existing
%      singleton*.
%
%      H = NIRS_AVGTRIAL_NOISE returns the handle to a new NIRS_AVGTRIAL_NOISE or the handle to
%      the existing singleton*.
%
%      NIRS_AVGTRIAL_NOISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NIRS_AVGTRIAL_NOISE.M with the given input arguments.
%
%      NIRS_AVGTRIAL_NOISE('Property','Value',...) creates a new NIRS_AVGTRIAL_NOISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MultiStim_View_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NIRS_AvgTrial_Noise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NIRS_AvgTrial_Noise

% Last Modified by GUIDE v2.5 04-Oct-2011 11:17:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NIRS_AvgTrial_Noise_OpeningFcn, ...
                   'gui_OutputFcn',  @NIRS_AvgTrial_Noise_OutputFcn, ...
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


% --- Executes just before NIRS_AvgTrial_Noise is made visible.
function NIRS_AvgTrial_Noise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NIRS_AvgTrial_Noise (see VARARGIN)

% Choose default command line output for NIRS_AvgTrial_Noise
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NIRS_AvgTrial_Noise wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NIRS_AvgTrial_Noise_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Noise_Callback(hObject, eventdata, handles)
% hObject    handle to Noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
name = get(gco,'displayname');
set(handles.Context_NAME,'Label',name);

% --------------------------------------------------------------------
function Context_Remove_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    %HOMERhandles = guihandles(guiHOMER); 
currentsub=1;
PMI = get(guiHOMER,'UserData'); 
cf = PMI{currentsub}.currentFile; 
name = get(gco,'displayname'); 
[tok,name]=strtok(name,'ch');
filenb = str2num(tok(5:end));
[tok,rem]=strtok(name,'ind');
set(gco,'color','y','marker','+','markersize',10);
plotLst = str2num(tok(4:end));
NChalf = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
if isempty(find(plotLst - NChalf < 0 ))
    plotLst = [plotLst;plotLst-NChalf];
else
    plotLst = [plotLst;plotLst+NChalf];
end
ind = str2num(rem(5:end));

noise = PMI{currentsub}.filenoise{filenb}.noise;
noise(ind,plotLst)= 1;
PMI{currentsub}.filenoise{filenb}.noise = noise;
set(guiHOMER,'UserData',PMI);

% --------------------------------------------------------------------
function Context_Restore_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData'); 
cf = PMI{currentsub}.currentFile; 
name = get(gco,'displayname');
[tok,name]=strtok(name,'ch');
filenb = str2num(tok(5:end));
[tok,rem]=strtok(name,'ind');
plotLst = str2num(tok(4:end));
idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList,plotLst, numel(PMI{currentsub}.color)/3);
set(gco,'color',PMI{currentsub}.color(idxc,:));
set(gco,'marker','none');
NChalf = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
if isempty(find(plotLst - NChalf < 0 ))
    plotLst = [plotLst;plotLst-NChalf];
else
    plotLst = [plotLst;plotLst+NChalf];
end
ind = str2num(rem(5:end));
noise = PMI{currentsub}.filenoise.noise{filenb};
noise(ind,plotLst)= 0;
PMI{currentsub}.filenoise{filenb}.noise = noise;
set(guiHOMER,'UserData',PMI);
% --------------------------------------------------------------------
function Context_NAME_Callback(hObject, eventdata, handles)
% hObject    handle to Context_NAME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in radio_viewnoise.
function radio_viewnoise_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewnoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_viewnoise

hch = get(gca,'children');   
for i = 1:numel(hch)
    mark = get(hch(i),'marker');
    if strcmp(mark,'+')
        if get(handles.radio_viewnoise,'value')
            set(hch(i),'visible','on');
        else
            set(hch(i),'visible','off');
        end
    end
end
