function varargout = GUI_ERLCOH_VIEW_Kmean(varargin)
% GUI_ERLCOH_VIEW_KMEAN M-file for GUI_ERLCOH_VIEW_Kmean.fig
%      GUI_ERLCOH_VIEW_KMEAN, by itself, creates a new GUI_ERLCOH_VIEW_KMEAN or raises the existing
%      singleton*.
%
%      H = GUI_ERLCOH_VIEW_KMEAN returns the handle to a new GUI_ERLCOH_VIEW_KMEAN or the handle to
%      the existing singleton*.
%
%      GUI_ERLCOH_VIEW_KMEAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ERLCOH_VIEW_KMEAN.M with the given input arguments.
%
%      GUI_ERLCOH_VIEW_KMEAN('Property','Value',...) creates a new GUI_ERLCOH_VIEW_KMEAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ERLCOH_VIEW_Kmean_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ERLCOH_VIEW_Kmean_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ERLCOH_VIEW_Kmean

% Last Modified by GUIDE v2.5 08-Oct-2019 14:26:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ERLCOH_VIEW_Kmean_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ERLCOH_VIEW_Kmean_OutputFcn, ...
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


% --- Executes just before GUI_ERLCOH_VIEW_Kmean is made visible.
function GUI_ERLCOH_VIEW_Kmean_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ERLCOH_VIEW_Kmean (see VARARGIN)

% Choose default command line output for GUI_ERLCOH_VIEW_Kmean
handles.output = hObject;
setappdata(0,'GUI_ERLCOH_VIEW',handles.GUI_ERLCOH_VIEW_Kmean);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_ERLCOH_VIEW_Kmean wait for user response (see UIRESUME)
% uiwait(handles.GUI_ERLCOH_VIEW_Kmean);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ERLCOH_VIEW_Kmean_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_BrowseSubjectList.
function btn_BrowseSubjectList_Callback(hObject, eventdata, handles)
% hObject    handle to btn_BrowseSubjectList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,filepath] = uigetfile({'*.xlsx;*.xls; Excel file';...
    ;'*.csv';'*.*'});
if ~file
    return
end

set(handles.edit_subjectList,'string',[filepath,file]);
[~,~,ext] =fileparts([filepath,file]);
if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
    [~,~, info]=xlsread(fullfile(filepath, file));
elseif strcmp(ext,'.csv')
    info = table2cell(readtable(fullfile(filepath, file),'Delimiter',';','ReadVariableNames',0));
end

optionERLCOH ='abs'; %'imag'
%elefile = [pathinfo,'Electrode_99ele.ele']
%check all zone list have the same elefile 
tic
groupid = cell2mat(info(2:end,3));
pathdata = info(2:end,1);
listfile = info(2:end,2); 
readflag = cell2mat(info(2:end,4));
optionERLCOH = info(2:end,5)
id = isnan(readflag);
if sum(id) %set missing to 0
readflag(find(id))=0;
end
 listsubject = [];
 idsubject = 1;
for ifile = 1:numel(pathdata)
 if readflag(ifile)
    load(fullfile(pathdata{ifile},[listfile{ifile},'.mat']),'-mat');
 
    if idsubject==1
        listelectrode = Args.ZoneList; %To ensure all subject have the same electrode order reference first one
    end
    layer = Args.layer;
    nlayer = numel(Args.layer);
    ntime =numel(Args.time);
    nele = numel(Args.ZoneList);
    nlink = (nele*nele-nele)/2;    
    load(fullfile(pathdata{ifile},[listfile{ifile},'.mat']),'-mat');
    a = reshape(MATCORR,nlayer*ntime, nlink);
      if numel(Args.ZoneList)~=numel(listelectrode)
         msgbox(['Check electrode compatibility ']);
         return
      end
          listsubject =  [listsubject;listfile(ifile)];
          caseERLCOH =  optionERLCOH{ifile}
      switch caseERLCOH
          case 'abs'              
                DATA{idsubject}.MATCORR = abs(a(:,:));
                DATA{idsubject}.Args=Args;
             	DATA{idsubject}.C = C;
                DATA{idsubject}.GR = groupid(ifile);
                DATA{idsubject}.name = listfile{ifile};
          case 'fisher'
                  a = abs(a(:,:));
                  min(a(:))
                  max(a(:))
                  id = find(a==1)
                  a(id) = 0.99; %avoid infinity response
                  allG1=abs(1/2.*(log((1+a)./(1-a))));
                  min(allG1(:))
                  max(allG1(:))
                  DATA{idsubject}.MATCORR =allG1;
                  DATA{idsubject}.Args=Args;
                  DATA{idsubject}.C = C;
                 DATA{idsubject}.GR = groupid(ifile);
                 DATA{idsubject}.name = listfile{ifile};
                          
          case 'imag'               
                DATA{idsubject}.MATCORR = imag(a(:,:));
                DATA{idsubject}.Args=Args;
             	DATA{idsubject}.C = C;
                DATA{idsubject}.GR = groupid(ifile);
                DATA{idsubject}.name = listfile{ifile};            

      end
       idsubject =  idsubject+1;
 end
end


toc
set(handles.popupmenu_subject,'string',listsubject);
set(handles.GUI_ERLCOH_VIEW_Kmean,'UserData',DATA);
%updateNetAllView(handles)
guidata(handles.GUI_ERLCOH_VIEW_Kmean, handles);


function edit_cmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cmax as text
%        str2double(get(hObject,'String')) returns contents of edit_cmax as a double

update_ERLCOHView(handles)
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



function edit_cmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cmin as text
%        str2double(get(hObject,'String')) returns contents of edit_cmin as a double

update_ERLCOHView(handles)
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



function edit_treshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_treshold as text
%        str2double(get(hObject,'String')) returns contents of edit_treshold as a double
update_ERLCOHView(handles)

% --- Executes during object creation, after setting all properties.
function edit_treshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_treshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xlszoneorder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xlszoneorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xlszoneorder as text
%        str2double(get(hObject,'String')) returns contents of edit_xlszoneorder as a double


% --- Executes during object creation, after setting all properties.
function edit_xlszoneorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xlszoneorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_subjectList_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subjectList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subjectList as text
%        str2double(get(hObject,'String')) returns contents of edit_subjectList as a double


% --- Executes during object creation, after setting all properties.
function edit_subjectList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subjectList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_subject.
function popupmenu_subject_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_subject contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_subject

update_ERLCOHView(handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_subject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_BrowseZoneOrder.
function btn_BrowseZoneOrder_Callback(hObject, eventdata, handles)
% hObject    handle to btn_BrowseZoneOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,filepath] = uigetfile({'*.xlsx;*.xls; Excel file';...
    ;'*.csv';'*.*'});
if ~file
    return
end

set(handles.edit_xlszoneorder,'string',fullfile(filepath,file))
update_ERLCOHView(handles)

% set(handles.edit_subjectList,'string',[filepath,file]);
% [~,~,ext] =fileparts([filepath,file]);
% if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
%     [~,~, info]=xlsread(fullfile(filepath, file));
% elseif strcmp(ext,'.csv')
%     info = table2cell(readtable(fullfile(filepath, file),'Delimiter',';','ReadVariableNames',0));
% end


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_clusterID.
function popupmenu_clusterID_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_clusterID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_clusterID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_clusterID
update_ERLCOHView(handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_clusterID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_clusterID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_modedisplay.
function popup_modedisplay_Callback(hObject, eventdata, handles)
% hObject    handle to popup_modedisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_modedisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_modedisplay

if get(handles.popup_modedisplay,'value')==1;
    set(handles.edit_definetime,'visible','on');
    set(handles.edit_frequency,'visible','on');
    set(handles.text_definetime,'string','Time');
    set(handles.text_definefreq,'string','Frequency');
    set(handles.popup_clusternb,'visible','off');
    set(handles.popupmenu_clusterID,'visible','off');
else
    set(handles.edit_definetime,'visible','off');
    set(handles.edit_frequency,'visible','off');
     set(handles.text_definetime,'string','Cluster nb');
     set(handles.text_definefreq,'string','Cluster id');
     set(handles.popup_clusternb,'visible','on');
     set(handles.popupmenu_clusterID,'visible','on');
     idsubject = get(handles.popupmenu_subject,'value');
     DATA = get(handles.GUI_ERLCOH_VIEW_Kmean,'UserData');
     if ~isempty( DATA{idsubject}.C)
         ClusterNb = numel(DATA{idsubject}.C)
         ClusterId = 1
         label  =[];
         for i=1:ClusterNb
         label = [label;{['Nb cluster',num2str(i)]}]
         end
        set(handles.popup_clusternb,'string',label )        
        set(handles.popupmenu_clusterID,'string','1')
     end
    
end
update_ERLCOHView(handles);

% --- Executes during object creation, after setting all properties.
function popup_modedisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_modedisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_definetime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_definetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_definetime as text
%        str2double(get(hObject,'String')) returns contents of edit_definetime as a double

update_ERLCOHView(handles)
% --- Executes during object creation, after setting all properties.
function edit_definetime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_definetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frequency as text
%        str2double(get(hObject,'String')) returns contents of edit_frequency as a double
update_ERLCOHView(handles)

% --- Executes during object creation, after setting all properties.
function edit_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_XLSconnectogram_Callback(hObject, eventdata, handles)
% hObject    handle to edit_XLSconnectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_XLSconnectogram as text
%        str2double(get(hObject,'String')) returns contents of edit_XLSconnectogram as a double


% --- Executes during object creation, after setting all properties.
function edit_XLSconnectogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_XLSconnectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_BrowseXLSconnectogram.
function btn_BrowseXLSconnectogram_Callback(hObject, eventdata, handles)
% hObject    handle to btn_BrowseXLSconnectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,filepath] = uigetfile({'*.xlsx;*.xls; Excel file';...
    ;'*.csv';'*.*'});
if ~file
    return
end

set(handles.edit_XLSconnectogram,'string',fullfile(filepath,file))
update_ERLCOHView(handles)


% --------------------------------------------------------------------
function Context_axes_matrices_Callback(hObject, eventdata, handles)
% hObject    handle to Context_axes_matrices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_axes_matrices_newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to context_axes_matrices_newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 GUI_LookMatrice_handles = plot_axes_MATRICE(handles,1);


% --------------------------------------------------------------------
function context_axes_connectogram_newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to context_axes_connectogram_newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 GUI_LookMatrice_handles = plot_axes_connectogramEEG(handles,1);

% --------------------------------------------------------------------
function Context_axes_connectogram_Callback(hObject, eventdata, handles)
% hObject    handle to Context_axes_connectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popup_clusternb.
function popup_clusternb_Callback(hObject, eventdata, handles)
% hObject    handle to popup_clusternb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_clusternb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_clusternb
1
id = get(handles.popup_clusternb,'value')
label = []
for i = 1:id
    label = [label;{['Cluster id', num2str(i)]}]
    end
set(handles.popupmenu_clusterID,'string',label)
set(handles.popupmenu_clusterID,'value',1)
update_ERLCOHView(handles)
% --- Executes during object creation, after setting all properties.
function popup_clusternb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_clusternb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function context_axes_timefreq_newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to context_axes_timefreq_newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 GUI_LookMatrice_handles = plot_axes_timefreq(GUI_LookMatrice_handles,1);

% --------------------------------------------------------------------
function Context_axes_timefreq_Callback(hObject, eventdata, handles)
% hObject    handle to Context_axes_timefreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_KMEAN.
function btn_KMEAN_Callback(hObject, eventdata, handles)
% hObject    handle to btn_KMEAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nbclustermax = str2num(get(handles.edit_cluster,'string'))
%mean all subject to compute Kmean
DATA = get(handles.GUI_ERLCOH_VIEW_Kmean,'UserData');
sumMATCORR = [];

for isubject=1:numel(DATA)
    sumMATCORR = cat(3,sumMATCORR,DATA{isubject}.MATCORR);   
end
nlayer = numel(DATA{isubject}.Args.layer);
ntime =numel(DATA{isubject}.Args.time);
nlink = size(DATA{isubject}.MATCORR,2)
X = nanmean(sumMATCORR,3);
for icluster= 1:nbclustermax
    icluster
    [C{icluster}.IDX,C{icluster}.mean,sumd,D] = kmeans(X,icluster);
    distsum(icluster) = nansum(sumd);
    C{icluster}.sumd = sumd;
end
Args = DATA{isubject}.Args
MATCORR = reshape(X, nlayer,ntime,nlink);
[pathin,filein,ext]=fileparts(get(handles.edit_subjectList,'string'))
[fileCOH,pathCOH] = uiputfile(fullfile(pathin,['kmean',filein,'.mat']));
 save(fullfile(pathCOH,fileCOH) ,'MATCORR','C','Args') ;




function edit_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cluster as text
%        str2double(get(hObject,'String')) returns contents of edit_cluster as a double


% --- Executes during object creation, after setting all properties.
function edit_cluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_displayoption.
function popup_displayoption_Callback(hObject, eventdata, handles)
% hObject    handle to popup_displayoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_displayoption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_displayoption
update_ERLCOHView(handles)

if get(handles.popup_displayoption,'value')==1
   set(handles.text_val,'visible','off')
    set(handles.edit_ptr,'visible','off')
elseif  get(handles.popup_displayoption,'value')==2
    set(handles.text_val,'visible','on')
    set(handles.edit_ptr,'visible','on')
elseif  get(handles.popup_displayoption,'value')==3
       set(handles.text_val,'visible','off')
    set(handles.edit_ptr,'visible','off')
elseif  get(handles.popup_displayoption,'value')==4
    set(handles.text_val,'visible','on')
    set(handles.edit_ptr,'visible','on')
end
   
    
% --- Executes during object creation, after setting all properties.
function popup_displayoption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_displayoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_FLIPSIGN.
function checkbox_FLIPSIGN_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_FLIPSIGN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_FLIPSIGN
update_ERLCOHView(handles)



function edit_ptr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ptr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ptr as text
%        str2double(get(hObject,'String')) returns contents of edit_ptr as a double
update_ERLCOHView(handles)

% --- Executes during object creation, after setting all properties.
function edit_ptr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ptr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_refresh.
function btn_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to btn_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
update_ERLCOHView(handles)
