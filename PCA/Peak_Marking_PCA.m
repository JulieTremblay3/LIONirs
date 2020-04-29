function varargout = Peak_Marking_PCA(varargin)
% PEAK_MARKING_PCA M-file for Peak_Marking_PCA.fig
%      PEAK_MARKING_PCA, by itself, creates a new PEAK_MARKING_PCA or raises the existing
%      singleton*.
%
%      H = PEAK_MARKING_PCA returns the handle to a new PEAK_MARKING_PCA or the handle to
%      the existing singleton*.
%
%      PEAK_MARKING_PCA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEAK_MARKING_PCA.M with the given input arguments.
%
%      PEAK_MARKING_PCA('Property','Value',...) creates a new PEAK_MARKING_PCA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Peak_Marking_PCA_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Peak_Marking_PCA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Peak_Marking_PCA

% Last Modified by GUIDE v2.5 16-Feb-2012 14:50:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Peak_Marking_PCA_OpeningFcn, ...
                   'gui_OutputFcn',  @Peak_Marking_PCA_OutputFcn, ...
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


% --- Executes just before Peak_Marking_PCA is made visible.
function Peak_Marking_PCA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Peak_Marking_PCA (see VARARGIN)

% Choose default command line output for Peak_Marking_PCA
handles.output = hObject;
handles.data = varargin{1};
handles.zone = varargin{2};
handles.t = varargin{3};
handles.listremove = varargin{4};
handles.dataPCA = zeros(size(handles.data));

resetview(handles)

set(handles.popupmenu_zone,'string',handles.zone.label)
% Update handles structure

guidata(hObject, handles);

% UIWAIT makes Peak_Marking_PCA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Peak_Marking_PCA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
    varargout{1} = handles.output;
    varargout{2} = handles.dataPCA;
    varargout{3} = handles.listremove;
    uiwait(hObject);
    
    %If 'hObject' is not a handle anymore, user has exited Dlg with 'X'
    if( ishandle(hObject) )
        handles = guidata( hObject );
        varargout{1} = handles.output;
        varargout{2} = handles.dataPCA;
        varargout{3} = handles.listremove;
        close(hObject);
    end
% --------------------------------------------------------------------
function Peakselection_Callback(hObject, eventdata, handles)
% hObject    handle to Peakselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function Remove_Callback(hObject, eventdata, handles)
% hObject    handle to Remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(handles.axes1,'CurrentPoint');
yli = get(handles.axes1,'YLim');
pos = get(handles.axes1,'CurrentPoint');
ind = find(pos(1,1)> handles.PeakMarked -1 &  pos(1,1)< handles.PeakMarked + 1);
handles.PeakMarked(ind) = [];
delete(handles.hPeakMarked(ind));
handles.hPeakMarked(ind) =[];
ylim([0.8,1.2]);
guidata(hObject, handles);


% --------------------------------------------------------------------
function Mark_Callback(hObject, eventdata, handles)
% hObject    handle to Mark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


yli = get(handles.axes1,'YLim');
pos = get(handles.axes1,'CurrentPoint');
ind = find(pos(1,1)>handles.t);
posx = ind(end);
axes(handles.axes1)
h = plot([handles.t(posx),handles.t(posx)],[yli(1), yli(2)]);
button = questdlg([num2str(pos(1,1)),'        ', num2str(pos(1,2))],'Confirm','yes','no','yes');
if strcmp(button,'yes')
    handles.PeakMarked = [handles.PeakMarked, posx];
    handles.hPeakMarked = [handles.hPeakMarked , h];
else
    delete(h)
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = []; %handles.PeakMarked;
guidata(gcf, handles);
uiresume(gcf)


% --- Executes on button press in btn_stopmarking.
function btn_stopmarking_Callback(hObject, eventdata, handles)
% hObject    handle to btn_stopmarking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = 0;
guidata(gcf, handles);
uiresume(gcf)


% --- Executes on button press in btn_marknoise.
function btn_marknoise_Callback(hObject, eventdata, handles)
% hObject    handle to btn_marknoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gca
plotbrowser('on')


% --------------------------------------------------------------------
function RemoveNoise_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% nbch = get(gco,'displayname')
% set(handles.RemoveNoisech,'label',['RemoveChannel : ',nbch])
% cmax = max(max(handles.data(:,find(handles.listremove))));
% cmin = min(min(handles.data(:,find(handles.listremove))));
% amp = max(abs([cmin,cmax]))
% ylim([-amp,amp])

% --------------------------------------------------------------------
function RemoveNoisech_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveNoisech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nbtoremove = get(gco,'displayname')
set(gco,'color','y','linestyle','-.')
handles.listremove(str2num(nbtoremove))= 0;
guidata(hObject, handles);

% --- Executes on button press in btn_mask_removech.
function btn_mask_removech_Callback(hObject, eventdata, handles)
% hObject    handle to btn_mask_removech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetview(handles)



% --------------------------------------------------------------------
function menu_Restore_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nbtoremove = get(gco,'displayname')
handles.listremove(str2num(nbtoremove))= 1;
guidata(hObject, handles);
resetview(handles);

% --- Executes on selection change in popupmenu_zone.
function popupmenu_zone_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_zone contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_zone
resetview(handles)

      

% --- Executes during object creation, after setting all properties.
function popupmenu_zone_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_PCA.
function btn_PCA_Callback(hObject, eventdata, handles)
% hObject    handle to btn_PCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

izone = get(handles.popupmenu_zone,'value');
label = get(handles.popupmenu_zone,'string');
label{izone}= [label{izone},'_PCA'];
set(handles.popupmenu_zone,'string',label);

plotLst=handles.zone.plotLst{izone};
d1 = handles.data(:,plotLst);
% %figure;plot(d1)
% A = handles.data(124,:)-1
% save('TESTTOPO.mat','A','-mat')
c = d1'*d1;
[v,s,foo]=svd(c);
svs = diag(s);
u = d1*v*inv(s);
lstSV =1; %composante 1
temp = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
Norm_Cov = d1- u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)'+1;
handles.dataPCA(:,plotLst) = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)'-1;%Norm_Cov;

guidata(hObject, handles);
resetview(handles);

% --- Executes on button press in btn_Artefact.
function btn_Artefact_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Artefact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

izone = get(handles.popupmenu_zone,'value');
plotLst=handles.zone.plotLst{izone};
label = get(handles.popupmenu_zone,'string');
label{izone}= [label{izone},'_bad']; 
set(handles.popupmenu_zone,'string',label);

handles.listremove(plotLst)=0;
guidata(hObject, handles);
resetview(handles)

function resetview(handles)
%All normal 
axes(handles.axes1);
cla
hold on 
plotLst = [];
colorchoice = [];
    for izone=numel(handles.zone.plotLst):-1:1
        plotLst =[plotLst;handles.zone.plotLst{izone}];
        for izoneplotLst=1:numel(handles.zone.plotLst{izone});
            colorchoice = [colorchoice;handles.zone.color(izone,:)];
        end
    end
for i = 1:numel(plotLst) 
       ch = plotLst(i);     
        h = plot(handles.t,handles.data(:,ch));
        set(h,'uicontext',handles.RemoveNoise);
        set(h,'displayname',num2str(ch));
         if handles.listremove(ch)
             set(h,'color',colorchoice(i,:));
         else
             if get(handles.btn_mask_removech,'value')
                set(h,'visible','off');
             else
                 set(h,'color','y');
             end
         end
end
xlim([handles.t(1),handles.t(end)]);
% ylim([0.8,1.2]);

%Region
axes(handles.axes2);
cla
izone = get(handles.popupmenu_zone,'value');
hold on
plotLstizone=handles.zone.plotLst{izone};
d1 = handles.data(:,plotLstizone);
for i=1:numel(plotLstizone)
    ch = plotLstizone(i);
    h = plot(handles.t,handles.data(:,ch));
    set(h,'displayname',num2str(ch),'color',handles.zone.color(izone,:));
end
c = d1'*d1;
[v,s,foo]=svd(c);
svs = diag(s);
u = d1*v*inv(s);
lstSV =1; %composante 1
temp = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
Norm_Cov = d1- u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)'+1;
h = plot(handles.t,Norm_Cov ,'k');
xlim([handles.t(1),handles.t(end)]);
% ylim([0.80,1.20]);

%Principal components display
axes(handles.axes_PCA);
cla
hold on
for i = 1:numel(plotLst) 
       ch = plotLst(i);
        h = plot(handles.t,handles.dataPCA(:,ch));
        set(h,'uicontext',handles.RemoveNoise);
        set(h,'displayname',num2str(ch));
         if handles.listremove(ch)
             set(h,'color',colorchoice(i,:));
         else
             if get(handles.btn_mask_removech,'value')
                set(h,'visible','off');
             else
                set(h,'color','y','visible','off');
             end
         end
end

xlim([handles.t(1),handles.t(end)]);

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
