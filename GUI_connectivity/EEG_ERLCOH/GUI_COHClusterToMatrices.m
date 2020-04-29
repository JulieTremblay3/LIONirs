function varargout = GUI_COHClusterToMatrices(varargin)
% GUI_COHCLUSTERTOMATRICES M-file for GUI_COHClusterToMatrices.fig
%      GUI_COHCLUSTERTOMATRICES, by itself, creates a new GUI_COHCLUSTERTOMATRICES or raises the existing
%      singleton*.
%
%      H = GUI_COHCLUSTERTOMATRICES returns the handle to a new GUI_COHCLUSTERTOMATRICES or the handle to
%      the existing singleton*.
%
%      GUI_COHCLUSTERTOMATRICES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_COHCLUSTERTOMATRICES.M with the given input arguments.
%
%      GUI_COHCLUSTERTOMATRICES('Property','Value',...) creates a new GUI_COHCLUSTERTOMATRICES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_COHClusterToMatrices_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_COHClusterToMatrices_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_COHClusterToMatrices

% Last Modified by GUIDE v2.5 07-Feb-2019 10:32:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_COHClusterToMatrices_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_COHClusterToMatrices_OutputFcn, ...
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


% --- Executes just before GUI_COHClusterToMatrices is made visible.
function GUI_COHClusterToMatrices_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_COHClusterToMatrices (see VARARGIN)

% Choose default command line output for GUI_COHClusterToMatrices
handles.output = hObject;
%setappdata(0,'GUI_COHClusterToMatrices',handles.GUI_COHClusterToMatrices)
setappdata(0,'gui_COHCluster',handles.figure1);
guiCOHc= getappdata(0,'gui_COHCluster');
data = []; %where gui COHERENCE data will be store
set(guiCOHc,'UserData',data);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_COHClusterToMatrices wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_COHClusterToMatrices_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_COHERENCE_CLUSTER_Callback(hObject, eventdata, handles)
% hObject    handle to edit_COHERENCE_CLUSTER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_COHERENCE_CLUSTER as text
%        str2double(get(hObject,'String')) returns contents of edit_COHERENCE_CLUSTER as a double


% --- Executes during object creation, after setting all properties.
function edit_COHERENCE_CLUSTER_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_COHERENCE_CLUSTER (see GCBO)
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
[file,path]=uigetfile('.mat')

guiCOHc= getappdata(0,'gui_COHCluster');
data = load([path,file])
set(handles.edit_COHERENCE_CLUSTER,'string',[path,file])
set(guiCOHc,'UserData',data);
labelcluster = [];
for i=1:numel(data.C)
    labelcluster =  [labelcluster;{['cluster',num2str(i)]}];  
end
  set(handles.popup_Cluster,'string',labelcluster)
 val=get(handles.popup_Cluster,'value')
 labelclusterid = [];
for i=1:numel(data.C{val})
    labelclusterid =  [labelclusterid;{['IDx',num2str(i)]}];  
end
set(handles.popup_clusterid,'string', labelclusterid)

 
updateview(handles)

function edit_link_Callback(hObject, eventdata, handles)
% hObject    handle to edit_link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_link as text
%        str2double(get(hObject,'String')) returns contents of edit_link as a double
updateview(handles)

% --- Executes during object creation, after setting all properties.
function edit_link_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_Cluster.
function popup_Cluster_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Cluster contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Cluster
set(handles.popup_clusterid,'value',1)
updateview(handles)

% --- Executes during object creation, after setting all properties.
function popup_Cluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function updateview(handles)
guiCOHc= getappdata(0,'gui_COHCluster');
data = get(guiCOHc,'UserData');

%affiche un lien de coherence
link = str2num(get(handles.edit_link,'string'))
axes(handles.axes_COHERENCE)
fs = 1/data.Args.dt;
t = 1/fs:1/fs:(1/fs*size(data.Args.coix,2));
if isfield(data.Args,'time')
    t = data.Args.time
end
period = data.Args.period
try
    if isfield(data.Args,'layer')
        imagesc(t ,data.Args.layer,abs(data.MATCORR(:,:,link)))
         axis xy
         ylabel('Frequency Hz')
         xlabel('Time (s)')
         
    else
        imagesc( t,log2(data.Args.period),data.MATCORR(:,:,link))
        Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period)))); % period ticks
        Fticks = 1./(Yticks.*data.Args.tRatio);
        set(gca,'YLim',log2([min(period),max(period)]), ...
        'YDir','reverse', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',num2str(Fticks'), ...
        'layer','top')
        ylabel('Frequency')
        caxis([0 1])
    end
catch
end

%nwavelet = 
ntime = size(data.Args.coix,2)
nlayer = size(data.Args.fwin,2)
%Affiche un cluster
axes(handles.axes_CLUSTER)
val=get(handles.popup_Cluster,'value')
IDX = data.C{val}.IDX
IDtF =reshape(IDX,nlayer,ntime)
if isfield(data.Args,'layer')
 imagesc(t ,data.Args.layer,IDtF)
 axis xy
 ylabel('Frequency Hz')
 xlabel('Time (s)')
else
        Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period)))); % period ticks
        period=data.Args.period;
        Fticks = 1./(Yticks.*data.Args.tRatio);
%         figure
        imagesc(t ,log2(data.Args.period),IDtF);
        set(gca,'YLim',log2([min(period),max(period)]), ...
    'YDir','reverse', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(Fticks'), ...
    'layer','top')
    ylabel('Frequency')
end
 val=get(handles.popup_Cluster,'value');
 maptmp = jet(max(IDX));
 colormap(maptmp);
 colorbar;
 labelclusterid = [];
for i=1:numel(data.C{val}.sumd)
    labelclusterid =  [labelclusterid;{['IDx',num2str(i)]}];  
end
set(handles.popup_clusterid,'string', labelclusterid)

%Afficher 
%La matrice pour un cluster
clusterid = get(handles.popup_clusterid,'val')
axes(handles.axes_MATRICES)

 id = 1;
NC=numel(data.Args.ZoneList)
nele = NC
matid = zeros(nele,nele);
for ielex=2:nele
    ielex
    ieley = 1;
    while ieley < ielex %| ieley==1
        %                     option{id}.ele1 =  listelectrode{ielex};
        %                     option{id}.ele2 =listelectrode{ieley};
        
        option{id}.matposition = [ielex,ieley];
        matid(ielex,ieley)=id;
        ieley = ieley + 1;
        id = id + 1;        
    end
end
idhalf = find(matid); %index dans la matrice 99x99
idoption = matid(idhalf); %index dans la matrice option
% 
data.C{val}.mean(clusterid,:)
matgr1 = zeros(nele,nele);
matgr1(idhalf)=data.C{val}.mean(clusterid,:);
matgr1 = matgr1 +flipud(rot90(matgr1))
imagesc(matgr1);
caxis([0 1])


% --- Executes on selection change in popup_clusterid.
function popup_clusterid_Callback(hObject, eventdata, handles)
% hObject    handle to popup_clusterid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_clusterid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_clusterid

updateview(handles)
% --- Executes during object creation, after setting all properties.
function popup_clusterid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_clusterid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_SaveMatrix.
function btn_SaveMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to btn_SaveMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiCOHc= getappdata(0,'gui_COHCluster');
data = get(guiCOHc,'UserData');
val=get(handles.popup_Cluster,'value')
clusterid = get(handles.popup_clusterid,'val') 
id = 1;
NC=numel(data.Args.ZoneList)
nele = NC
matid = zeros(nele,nele);
for ielex=2:nele
    ielex
    ieley = 1;
    while ieley < ielex %| ieley==1
        %                     option{id}.ele1 =  listelectrode{ielex};
        %                     option{id}.ele2 =listelectrode{ieley};
        
        option{id}.matposition = [ielex,ieley];
        matid(ielex,ieley)=id;
        ieley = ieley + 1;
        id = id + 1;        
    end
end
idhalf = find(matid); %index dans la matrice 99x99
idoption = matid(idhalf); %index dans la matrice option
% 
data.C{val}.mean(clusterid,:)
meancorr = zeros(nele,nele);
meancorr(idhalf)=data.C{val}.mean(clusterid,:);
meancorr = meancorr +flipud(rot90(meancorr))
   matcorr = meancorr;
   ZoneList = data.Args.ZoneList
   pathout = data.Args.pathout
   filloutput = data.Args.filloutput
    ntime = size(data.Args.coix,2)
    nlayer = size(data.Args.fwin,2)
   IDX = data.C{val}.IDX==clusterid;
   IDtF =reshape(IDX,nlayer,ntime);
   freqmean = round(mean(find(mean(IDtF'))))
   CenterFreq =1/data.Args.period(freqmean);
   label = sprintf('%s%02.0f%s%02.0f%s%f%','cluster',val,'idx',clusterid,'mfreq',CenterFreq )
   save(fullfile(pathout,[filloutput,'COH_HBO',label,'.mat']),'ZoneList','matcorr','meancorr');
