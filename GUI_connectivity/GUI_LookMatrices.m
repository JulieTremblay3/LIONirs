function varargout = GUI_LookMatrices(varargin)
% GUI_LOOKMATRICES M-file for GUI_LookMatrices.fig
%      GUI_LOOKMATRICES, by itself, creates a new GUI_LOOKMATRICES or raises the existing
%      singleton*.
%
%      H = GUI_LOOKMATRICES returns the handle to a new GUI_LOOKMATRICES or the handle to
%      the existing singleton*.
%
%      GUI_LOOKMATRICES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_LOOKMATRICES.M with the given input arguments.
%
%      GUI_LOOKMATRICES('Property','Value',...) creates a new GUI_LOOKMATRICES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_LookMatrices_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_LookMatrices_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_LookMatrices

% Last Modified by GUIDE v2.5 01-Jul-2020 11:18:42

% Begin initialization code - DO NOT EDITspm
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_LookMatrices_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_LookMatrices_OutputFcn, ...
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


% --- Executes just before GUI_LookMatrices is made visible.
function GUI_LookMatrices_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_LookMatrices (see VARARGIN)

% Choose default command line output for GUI_LookMatrices
handles.output = hObject;
setappdata(0,'GUI_LookMat',handles.GUI_LookMat)


newmfile = mfilename('fullpath');
 [path,name,ext1] = fileparts(newmfile);
           listtemplate = [];

listing = dir([path,filesep,'Template']);
 for ilisting=1:numel(listing)
     if numel(listing(ilisting).name)>3
        if strcmp(listing(ilisting).name(end-3:end),'.mat')
          listtemplate = [listtemplate,{listing(ilisting).name(1:end-4)}];
        end
     end
 end
set(handles.popup_template2dmap,'string', listtemplate)
% Update handles structure
guidata(handles.GUI_LookMat, handles);

% UIWAIT makes GUI_LookMatrices wait for user response (see UIRESUME)
% uiwait(handles.GUI_LookMat);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_LookMatrices_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_subjetxls_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subjetxls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subjetxls as text
%        str2double(get(hObject,'String')) returns contents of edit_subjetxls as a double


% --- Executes during object creation, after setting all properties.
function edit_subjetxls_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subjetxls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Btn_browseList.
function Btn_browseList_Callback(hObject, eventdata, handles)
% hObject    handle to Btn_browseList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,filepath] = uigetfile({'*.xlsx;*.xls; Excel file';...
    ;'*.csv';'*.txt';, '*.*'});
if ~file
    return
end
set(handles.edit_subjetxls,'string',[filepath,file]);
[~,~,ext] =fileparts([filepath,file]);
if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
[~,~, info]=xlsread(fullfile(filepath, file));
elseif strcmp(ext,'.csv')
    info = table2cell(readtable(fullfile(filepath, file),'Delimiter',';','ReadVariableNames',0));
elseif strcmp(ext,'.txt')
   [~,~, info]= readtxtfile_asxlsread([filepath,file]);
end


groupeall = [];
for isubject=2:size(info,1)

    id = isubject-1;
     try
        MAT = load(fullfile(info{isubject,1}, [info{isubject,2}]));
    catch
        disp(['ERROR loading' ,fullfile(info{isubject,1}, info{isubject,2})]);
     end
    if isfield( MAT,'ZoneList')
         DATA{id}.ZoneList = MAT.ZoneList;
    end

%     idbad =find(isnan(MAT.meancorr))
%     MAT.meancorr(idbad)=0;
    DATA{id}.MAT = MAT.meancorr;
    
%     idbad =find(isnan(MAT.meancorrPearsonFisher))
%     MAT.meancorrPearsonFisher(idbad)=0;
%     DATA{id}.MAT = MAT.meancorrPearsonFisher
    DATA{id}.path = info{isubject,1};
    DATA{id}.name = info{isubject,2};
    DATA{id}.MATtrial =  MAT.matcorr;
    DATA{id}.GR = info{isubject,4};
    infocov = [];
    if size(info,2) > 4
        %Define Additional covariable 
        idcov = 1 ;       
        for icolumn = 5:size(info,2)
            try
                 eval( ['DATA{id}.cov',num2str(idcov),' = ', num2str(info{isubject,icolumn}),';']);
            catch
                 eval( ['DATA{id}.cov',num2str(idcov),' = nan;']);
            end
            idcov=idcov+1;
            if isubject==2
            infocov = [infocov ;info(1,icolumn)];
            end
        end
    end
      if isubject==2
        set(handles.popupmenu_covariable,'string',infocov );
      end
      %If format A a1b2... or D01 for NIRx

  try
        load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');
    catch
        disp(['ERROR loading' ,fullfile(info{isubject,1}, info{isubject,3})]);
  end
        name = MAT.ZoneList{1};
       if strcmp(name(1:2),'D0')
          DATA{id}.System = 'NIRx';
      else
           DATA{id}.System = 'ISS Imagent';
   %        DATA{id}.System = 'ISS' ;           
      end  
  
  
    names = fieldnames(zone);
    for iname = 1:numel(names)
        eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname},';']);
    end
    list_subject{id} =DATA{id}.name;
    groupeall = [groupeall; info{isubject,4}];
end
%avg zone groupe 
idnew = size(info,1)-1; %placer les moyennes a la fin des sujets existants. 
for igroupe = 1:max(groupeall)

   idsubject = find(groupeall==igroupe);
   if ~isempty(idsubject)
   idlabelall= DATA{idsubject(1)}.zone.label; %zone premier sujet du groupe                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
   MATAVGALL = zeros(numel(idlabelall),numel(idlabelall),numel(idsubject));
   for isubject=1:numel(idsubject)
        ML=DATA{idsubject(isubject)}.zone.ml;
         
        List=strvcat(DATA{idsubject(isubject)}.ZoneList);
        MATAVG = zeros(numel(idlabelall));
        MAT = DATA{idsubject(isubject)}.MAT;
        idlist = [];
        idlabel=[];
        idzone =[];
        labelzone = DATA{idsubject(isubject)}.zone.label;
      
        for adji = 1:numel(idlabelall)
            for adjj = 1:numel(idlabelall)   
                labelzone = idlabelall{adji};
                x = strmatch({labelzone} ,idlabelall, 'exact');             
                labelzone = idlabelall{adjj};
                y = strmatch({labelzone} ,idlabelall, 'exact');
                if isempty(x)|isempty(y)
                    msgbox('problem zone in subject')
                end
                chzone = DATA{isubject}.zone.plotLst{x};
                idlisti = [];
                for ichzone = 1:numel(chzone);
                    ich = chzone(ichzone);
                    switch  DATA{id}.System
                        case 'ISS Imagent'
                            strDet = SDDet2strboxy_ISS(ML(ich,2));
                            strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                            idch = strmatch([strDet, ' ',strSrs ],List,'exact');                       
                        case 'NIRx'                          
                            strDet = SDDet2strboxy(ML(ich,2));
                            strSrs = SDPairs2strboxy(ML(ich,1));
                            idch = strmatch([strDet, ' ',strSrs ],List,'exact');                                                            
                        otherwise                            
                            strDet = SDDet2strboxy_ISS(ML(ich,2));
                            strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                            idch = strmatch([strDet, ' ',strSrs ],List,'exact');                          
                    end
                      idlisti = [idlisti, idch];
                end
                if numel(x)>1
                    x = x(1);
                    msgbox('Attention the zone are duplicated')
                else
                    DATA{idsubject(isubject)}.zone.chMAT{x} = idlisti;
                end
                
                chzone = DATA{id}.zone.plotLst{y};
                idlistj = [];
                for ichzone = 1:numel(chzone)
                    ich = chzone(ichzone);
                    switch  DATA{id}.System
                        case 'ISS Imagent'
                            strDet = SDDet2strboxy_ISS(ML(ich,2));
                            strSrs = SDPairs2strboxy_ISS(ML(ich,1));      
                            idch = strmatch([strDet, ' ',strSrs ],List,'exact');   
                      case 'NIRx'   
                              strDet = SDDet2strboxy(ML(ich,2));
                            strSrs = SDPairs2strboxy(ML(ich,1));
                            idch = strmatch([strDet, ' ',strSrs ],List,'exact');  
                        otherwise
                           strDet = SDDet2strboxy_ISS(ML(ich,2));
                            strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                            idch = strmatch([strDet, ' ',strSrs ],List,'exact'); 
                    end
                    idlistj = [idlistj, idch];
                end                       
              if isempty(idlisti)|isempty(idlistj)  
                   MATAVG(adji,adjj)=nan;
              else
                  try
                  temp = MAT(idlisti, idlistj);
                   MATAVG(adji,adjj)= nanmean(temp(:)); 
                  catch
                      1
                  end
              end         
        end
 end
    MATAVGALL(:,:,isubject) = MATAVG;
   end
       %Entrer les valeurs moyenne du groupe comme un sujet, les nouveaux canaux sont faux, 
       %il sont la pour l'uniformiter des données afin référer au zone.           
       ZoneList = [];
       plot= [];
       plotLst = [];
       idnew =  idnew + 1;  
       for izoneList = 1:size(MATAVGALL,1)
           MLfake(izoneList,1) = izoneList;%source
           MLfake(izoneList,2) = 1; %detecteur
           MLfake(izoneList,3) = 1;
           MLfake(izoneList,4) = 1;        
           strDet = SDDet2strboxy_ISS(MLfake(izoneList,2));
           strSrs = SDPairs2strboxy_ISS(MLfake(izoneList,1)); 
           ZoneList{izoneList,1}=[strDet,' ', strSrs];
           plot{izoneList} = [izoneList,1];
           plotLst{izoneList} = [izoneList];
         
       end
       
       DATA{idnew}.ZoneList = ZoneList;
       DATA{idnew}.MAT =  nanmean(MATAVGALL,3);
       DATA{idnew}.name = ['AVG groupe ',num2str(igroupe)];
       DATA{idnew}.MATtrial =  MATAVGALL;
       DATA{idnew}.GR = 0;     
       DATA{idnew}.System = 'ISS';
       DATA{idnew}.zone.plot = plot;
       DATA{idnew}.zone.plotLst = plotLst;
       DATA{idnew}.zone.label = idlabelall;
       DATA{idnew}.zone.color = DATA{idsubject(1)}.zone.color; 
       DATA{idnew}.zone.ml = MLfake;
       DATA{idnew}.zone.chMAT = plotLst;
       list_subject{idnew} =DATA{idnew}.name;
    
       clear ZoneList       
       clear MATAVGALL
       clear MLfake
       clear MAT    
      end
end
for  igroupe = 1:max(groupeall)
    	idsubject = find(groupeall==igroupe);
        if ~isempty(idsubject)
              idnew = idnew +1;
        for isubject = 1:numel(idsubject)
             MATall(:,:,isubject) = DATA{idsubject(isubject)}.MAT;
        end    
       DATA{idnew}.ZoneList = DATA{1}.ZoneList;
       DATA{idnew}.MAT =  nanmean(MATall,3);
       DATA{idnew}.name = ['AVG ch groupe ',num2str(igroupe)];
       DATA{idnew}.MATtrial =  MATall;
       DATA{idnew}.GR = 0;       
       DATA{idnew}.System = 'ISS';
       DATA{idnew}.zone.plot =DATA{1}.zone.plot;
       DATA{idnew}.zone.plotLst = DATA{1}.zone.plotLst;
       DATA{idnew}.zone.label =  DATA{1}.zone.label;
       DATA{idnew}.zone.color = DATA{idsubject(1)}.zone.color; 
       DATA{idnew}.zone.ml = DATA{idsubject(1)}.zone.ml;
       DATA{idnew}.zone.pos = DATA{idsubject(1)}.zone.pos;
       DATA{idnew}.zone.chMAT =  DATA{idsubject(1)}.zone.chMAT;
       list_subject{idnew} =DATA{idnew}.name;
        end
        
end%AVG groupe all
  clear MATall 

%Update data
set(handles.GUI_LookMat,'UserData',DATA);
set(handles.popup_listsujet,'string',list_subject);
set(handles.popup_listsujet,'value',1);
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles);
msgbox('Load Data Completed')

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have ahwhite background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_Next.
function btn_Next_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.listbox_zone,'string');
value = get(handles.listbox_zone,'value');
listselected = get(handles.listbox_selectedzone,'string');
listselected = [listselected;{list{value}} ];
 set(handles.listbox_selectedzone,'string',listselected);
 if value<numel(list)
     set(handles.listbox_zone,'value',value+1);
 end
 
alllist = [];
for ik = 1:numel(listselected)
    for jk = 1:numel(listselected)
        alllist = [alllist, {[listselected{ik},'/' ,listselected{jk}]}]
    end    
end
set(handles.listbox_subzone,'string',alllist)
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles);

% --- Executes on selection change in popup_listsujet.
function popup_listsujet_Callback(hObject, eventdata, handles)
% hObject    handle to popup_listsujet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_listsujet contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_listsujet
guidata(handles.GUI_LookMat, handles);
%plot_axesAij(handles);
updateNetAllView(handles)
% --- Executes during object creation, after setting all properties.
function popup_listsujet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_listsujet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in listbox_zone.
function listbox_zone_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_zone contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_zone


% --- Executes during object creation, after setting all properties.
function listbox_zone_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Pre.
function btn_Pre_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listselected = get(handles.listbox_selectedzone,'string');
value = get(handles.listbox_selectedzone,'value');
listselected(value)=[];
if value>1
    set(handles.listbox_selectedzone,'value',value-1);
end
set(handles.listbox_selectedzone,'string',listselected);
alllist = [];
for ik = 1:numel(listselected)
    for jk = 1:numel(listselected)
        alllist = [alllist, {[listselected{ik},'/' ,listselected{jk}]}]
    end    
end
set(handles.listbox_subzone,'string',alllist)
guidata(handles.GUI_LookMat, handles);
%plot_axesAij(handles);
updateNetAllView(handles)

% --- Executes on selection change in listbox_selectedzone.
function listbox_selectedzone_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_selectedzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_selectedzone contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_selectedzone


% --- Executes during object creation, after setting all properties.
function listbox_selectedzone_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_selectedzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function btn_Next_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btn_Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu_view.
function popupmenu_view_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_view contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_view
guidata(handles.GUI_LookMat, handles);
%plot_axesAij(handles)
updateNetAllView(handles)
% --- Executes during object creation, after setting all properties.
function popupmenu_view_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_viewhistogramme.
function popupmenu_viewhistogramme_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_viewhistogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_viewhistogramme contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_viewhistogramme
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_viewhistogramme_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_viewhistogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)

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
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)
%plot_axesAij(handles)

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


% --- Executes on button press in btn_VIEWNet.
function btn_VIEWNet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_VIEWNet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% h=GUI_ViewNetwork
% setappdata(0,'GUI_ViewNetwork',h);
updateNetAllView(handles);
plot_axesviewNET(handles);


% --------------------------------------------------------------------
function menu_Settings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Export_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edit_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_threshold as a double

updateNetAllView(handles)

% --- Executes during object creation, after setting all properties.
function edit_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_Cor_Back.
function radio_Cor_Back_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Cor_Back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Cor_Back


% --- Executes on button press in radio_Cor_Front.
function radio_Cor_Front_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Cor_Front (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Cor_Front


% --- Executes on button press in radio_Sag_Right.
function radio_Sag_Right_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Sag_Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Sag_Right


% --- Executes on button press in radio_Sag_Left.
function radio_Sag_Left_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Sag_Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Sag_Left


% --- Executes on button press in radio_Axi_Down.
function radio_Axi_Down_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Axi_Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Axi_Down


% --- Executes on button press in radio_Axi_Up.
function radio_Axi_Up_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Axi_Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Axi_Up


% --- Executes on selection change in listbox_subzone.
function listbox_subzone_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_subzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_subzone contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_subzone


% --- Executes during object creation, after setting all properties.
function listbox_subzone_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_subzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 

% --- Executes on button press in btn_next_subzone.
function btn_next_subzone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_next_subzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.listbox_subzone,'string');
value = get(handles.listbox_subzone,'value');
listselected = get(handles.listbox_subzone_selected,'string');
listselected = [listselected;{list{value}} ];
set(handles.listbox_subzone_selected,'string',listselected);
 if value<numel(listselected)
     set(handles.listbox_subzone_selected,'value',value+1);
 end
guidata(handles.GUI_LookMat, handles);

% --- Executes on button press in btn_pre_subzone.
function btn_pre_subzone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_pre_subzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listselected = get(handles.listbox_subzone_selected,'string');
value = get(handles.listbox_subzone_selected,'value');
listselected(value)=[];
if value>1
    set(handles.listbox_subzone_selected,'value',value-1);
end
set(handles.listbox_subzone_selected,'string',listselected);

% --- Executes on selection change in listbox_subzone_selected.
function listbox_subzone_selected_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_subzone_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_subzone_selected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_subzone_selected


% --- Executes during object creation, after setting all properties.
function listbox_subzone_selected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_subzone_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_individualplot.
function radio_individualplot_Callback(hObject, eventdata, handles)
% hObject    handle to radio_individualplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_individualplot


% --- Executes on button press in btn_graph_theory_measures.
function btn_graph_theory_measures_Callback(hObject, eventdata, handles)
% hObject    handle to btn_graph_theory_measures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try %Installation BNC message
DATA = get(handles.GUI_LookMat,'UserData');
thr = str2num(get(handles.edit_threshold,'string'));
if 1
idzone = [];
idlist = [];
idlabelall = [];


%POUR IDENTIFICATION DE TOUTE LES ZONES
id =1;
List=   strvcat(DATA{id}.ZoneList);
ML = DATA{id}.zone.ml;
  listok = get(handles.listbox_selectedzone,'string');
    idlist = [];
    idlabel=[];
    idzone =[];
    
    for ilistzone = 1:numel(listok)
    for izone = 1:numel(DATA{id}.zone.plotLst)
        chzone = DATA{id}.zone.plotLst{izone}
        labelzone = DATA{id}.zone.label{izone}
         x = strmatch({labelzone} , {listok{ilistzone}}, 'exact')
         if ~isempty(x)
%             for ichzone = 1:numel(chzone)
%                 ich = chzone(ichzone);
%                 if strcmp(DATA{id}.System,'ISS')
%                     strDet = SDDet2strboxy_ISS(ML(ich,2));
%                     strSrs = SDPairs2strboxy_ISS(ML(ich,1));
%                     idch = strmatch([strDet, ' ',strSrs ],List,'exact');
%                 end
                 idch=DATA{id}.zone.chMAT{izone};
                 idlist = [idlist, idch];
               % if ichzone==1
               %     idzone =[idzone, izone];
               % else
                    idzone =[idzone,izone, zeros(1,numel(idch)-1)];
              %  end
%             end
            idlabel = [idlabel, {[DATA{id}.zone.label{izone}, sprintf('_%03.0f',ilistzone)]}];
         end
    end
    end

    %Global Efficiency
    MatCTL =  DATA{7}.MAT(idlist,idlist);
    MatCHD =  DATA{8}.MAT(idlist,idlist);
    labelCTL = [DATA{7}.name,' CTL']
    labelCHD = [DATA{8}.name,' CHD']
    
   topthr = str2num(get(handles.edit_cmax,'string'))
    tr = 0.01:0.01:topthr;
    id = 1;
    listoption  = get(handles.popupmenu_GraphTheoryOption,'string');
    valoption = get(handles.popupmenu_GraphTheoryOption,'value');
    metric =listoption{valoption};
%     metric = 'Global Efficiency'
%     metric ='Clustering coefficient'
%     metric = 'Characteristic path length'
%     metric = 'Local Efficiency'




    for itr=tr
  

        
        %global efficiency
        switch metric
            case 'Global Efficiency'
                  
                Mat = threshold_absolute( MatCTL, itr);    
              
                   GE_CTL(id)=mean(efficiency_bin(Mat,1));
                   Mat = threshold_absolute(   MatCHD, itr);
                   GE_CHD(id)=mean(efficiency_bin(Mat,1));
%                    Mat = threshold_absolute( MatCTL, itr);           
%                    GE_CTL(id)=efficiency_bin(Mat,0);
%                    Mat = threshold_absolute(   MatCHD, itr);
%                    GE_CHD(id)=efficiency_bin(Mat,0);
                   
                   id = id +1;
            case 'Clustering coefficient'
                   Mat = threshold_absolute( MatCTL, itr);           
                   GE_CTL(id)=mean(clustering_coef_bu(Mat));
                   Mat = threshold_absolute(   MatCHD, itr);
                   GE_CHD(id)=mean(clustering_coef_bu(Mat));
                   id = id +1;
            case 'Characteristic path length'
                   Mat = threshold_absolute( MatCTL, itr);           
                   GE_CTL(id)= charpath(distance_bin(Mat));
                   Mat = threshold_absolute(   MatCHD, itr);
                   GE_CHD(id)  = charpath(distance_bin(Mat));
                   id = id +1;      
            case 'Local Efficiency'
                   Mat = threshold_absolute( MatCTL, itr);           
                   GE_CTL(:,id)=efficiency_bin(Mat,1);
                   Mat = threshold_absolute(   MatCHD, itr);
                   GE_CHD(:,id)=efficiency_bin(Mat,1);
                   id = id +1;
        end 

    end
    figure;hold on
    plot( tr,GE_CTL,'displayname',labelCTL,'linewidth',4,'color','b')
    plot( tr,GE_CHD,'displayname', labelCHD ,'linewidth',4,'color','r')
    xlabel('Threshold','fontsize',20)
    ylabel(metric,'fontsize',20)
    set(gca,'fontsize',20) 
    xlim([0,max(tr)])
    ylim([0,max([GE_CTL,GE_CHD])])
    clear GE_CHD;clear GE_CTL;
    if strcmp(metric,'Local Efficiency')
    figure;hold on
    plot( tr,GE_CTL-GE_CHD,'displayname',['CTL-CHD'],'linewidth',1,'color','b')
     xlabel('Threshold','fontsize',20)
    ylabel(metric,'fontsize',20)
    set(gca,'fontsize',20) 
    xlim([0,max(tr)])
    ylim([min([GE_CTL(:);GE_CHD(:)]),max([GE_CTL(:);GE_CHD(:)])])
    clear GE_CHD;clear GE_CTL;
    end
%     plot( tr,GE_CHD,'displayname', labelCHD ,'linewidth',4,'color','r')
    %

    
end

for isubject = 1:numel(DATA)
        Mat =  DATA{isubject}.MAT;
	    Mat = threshold_absolute(Mat, thr);
        id = find(isnan(Mat));
        if ~isempty(id);Mat(id)=0;end;
	    DATA{isubject}.measure.GlobalEfficiency  = efficiency_bin(Mat,0);            
        DATA{isubject}.measure.ClusteringCoef = mean(clustering_coef_bu(Mat));        
        DATA{isubject}.measure.node.LocalEfficiency = efficiency_bin(Mat,1); 
        DATA{isubject}.measure.node.Degree = degrees_und(Mat);     
        DATA{isubject}.measure.Mean = median(Mat(:));    
end
histogramme_name = get(handles.popupmenu_viewhistogramme,'string')
if ~iscell(histogramme_name)
histogramme_name = [histogramme_name;{'subject pearson'};{'GlobalEfficiency'};{'ClusteringCoef'}]
end
set(handles.popupmenu_viewhistogramme,'string', histogramme_name);

%submodule calculation
for isubject = 1:numel(DATA)

    for izone = 1:numel(DATA{isubject}.zone.plotLst)
        Mat = DATA{isubject}.MAT;
        Mat = threshold_absolute(Mat,thr);
        id = find(isnan(Mat));
        if ~isempty(id);Mat(id)=0;end;
        idzone = DATA{isubject}.zone.chMAT{izone};
        submat = Mat(idzone,idzone);
        
        DATA{isubject}.zone.GlobalEfficiency(izone)  = efficiency_bin(submat,0);
        DATA{isubject}.zone.ClusteringCoef(izone) = mean(clustering_coef_bu(submat));
        DATA{isubject}.zone.Mean(izone)= median(submat(:));
    end
end
 %module_degree_zscore
 
set(handles.GUI_LookMat,'UserData',DATA);
%EXCEL EXPORT 
DATA = get(handles.GUI_LookMat,'UserData');
thr = get(handles.edit_threshold,'string');
filename = ['GraphTheory Measures ', thr,'.xls' ]
[file,path] = uiputfile(filename,'Save file name');
  
for isubject = 1:numel(DATA)
       if strmatch(DATA{isubject }.name,'AVG groupe 1','exact')
        break
    end
     nbsubject = isubject;
end
xlssheet{1,1} = 'Subject'
xlssheet{1,2} = 'Groupe'
xlssheet{1,3} = 'Global Connectivity'
xlssheet{1,4} = 'Clustering coefficient'

icol = 1
for isubject=1:nbsubject
    xlssheet{isubject+1,1} = DATA{isubject}.name;
    xlssheet{isubject+1,2} = DATA{isubject}.GR;
    xlssheet{isubject+1,3} = DATA{isubject}.measure.GlobalEfficiency;
    xlssheet{isubject+1,4} = DATA{isubject}.measure.ClusteringCoef;
    xlssheet{isubject+1,5} = DATA{isubject}.measure.Mean;
end   
    listok = get(handles.listbox_selectedzone,'string')
    icol = 5;
for ilistzone = 1:numel(listok)
    for isubject=1:nbsubject       
        idlabelall = [];
        for izone = 1:numel(DATA{isubject}.zone.plotLst)
           idlabelall = [idlabelall, {DATA{isubject}.zone.label{izone}}]
        end
        x = strmatch(listok{ilistzone},idlabelall ,  'exact')        
        xlssheet{1,icol} = [listok{ilistzone}, ' GlobalEfficiency'];
         
         if ~isempty(x)
             xlssheet{isubject+1,icol} = DATA{isubject}.zone.GlobalEfficiency(x);
         end
    end
    icol = icol + 1;
end

for ilistzone = 1:numel(listok)
    for isubject=1:nbsubject       
        idlabelall = [];
        for izone = 1:numel(DATA{isubject}.zone.plotLst)
           idlabelall = [idlabelall, {DATA{isubject}.zone.label{izone}}]
        end
        x = strmatch(listok{ilistzone},idlabelall ,  'exact')        
        xlssheet{1,icol} = [listok{ilistzone}, ' Clustering coefficient'];
         
         if ~isempty(x)
             xlssheet{isubject+1,icol} = DATA{isubject}.zone.ClusteringCoef(x);
         end
    end
    icol = icol + 1;
end

for ilistzone = 1:numel(listok)
    for isubject=1:nbsubject       
        idlabelall = [];
        for izone = 1:numel(DATA{isubject}.zone.plotLst)
           idlabelall = [idlabelall, {DATA{isubject}.zone.label{izone}}]
        end
        x = strmatch(listok{ilistzone},idlabelall ,  'exact')        
        xlssheet{1,icol} = [listok{ilistzone}, ' Mean'];
         
         if ~isempty(x)
             xlssheet{isubject+1,icol} = DATA{isubject}.zone.Mean(x);
         end
    end
    icol = icol + 1;
end

xlswrite(fullfile(path,file), xlssheet)
  catch
                    msgbox({'Please Install Brain Connectivity Toolbox' 
                    'Complex network measures of brain connectivity: Uses and interpretations.'
                    'Rubinov M, Sporns O (2010) NeuroImage 52:1059-69. '})
   end
%msgbox('DONE')

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_Histogramme_Callback(hObject, eventdata, handles)
% hObject    handle to context_Histogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NewFigure_Callback(hObject, eventdata, handles)
% hObject    handle to NewFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_histogrammeAij(handles,1);


% --- Executes on button press in btn_exportGraphMeasures.
function btn_exportGraphMeasures_Callback(hObject, eventdata, handles)
% hObject    handle to btn_exportGraphMeasures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function context_axes_AJD_Callback(hObject, eventdata, handles)
% hObject    handle to context_axes_AJD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 
% --------------------------------------------------------------------
function Make_newfigure_axesAJD_Callback(hObject, eventdata, handles)
% hObject    handle to Make_newfigure_axesAJD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_axesAij(handles,1);


% --- Executes on button press in radio_drawlink.
function radio_drawlink_Callback(hObject, eventdata, handles)
% hObject    handle to radio_drawlink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_drawlink


% --------------------------------------------------------------------
function Menu_exportSeedFigure_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_exportSeedFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%  

label = get(handles.Menu_exportSeedFigure,'label')
[token, remain] = strtok(label,' ')

DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value');
filename = get(handles.popup_listsujet, 'string')
filename{id}
MAT = DATA{id}.MAT;
zonelist = DATA{id}.ZoneList;
A = [MAT(:,str2num(token))]';
[file,path]= uiputfile([filename{id},remain,'.mat'])
save([path,file],'A','zonelist','-mat')



% --------------------------------------------------------------------
function context_mat_Callback(hObject, eventdata, handles)
% hObject    handle to context_mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% posFig = get(hObject, 'Position')
% posFig = get(gca, 'Position')
% h=get(gca,'children')
% get(h(5))

d1=get(gco,'displayname')
set(handles.Menu_exportSeedFigure,'label',d1)


% --- Executes on button press in btn_STat.
function btn_STat_Callback(hObject, eventdata, handles)
% hObject    handle to btn_STat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_StatMatrices


% --------------------------------------------------------------------
function menu_NANmat_Callback(hObject, eventdata, handles)
% hObject    handle to menu_NANmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

label = get(handles.Menu_exportSeedFigure,'label')
[token, remain] = strtok(label,' ')

DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value');
filename = get(handles.popup_listsujet, 'string')
filename{id}
DATA{id}.MAT(:,str2num(token))=nan;
DATA{id}.MAT(str2num(token),:)=nan;
set(handles.GUI_LookMat,'UserData',DATA);
updateNetAllView(handles);




% --- Executes on button press in btn_savematrices.
function btn_savematrices_Callback(hObject, eventdata, handles)
% hObject    handle to btn_savematrices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
id = get(handles.popup_listsujet, 'value');
filename = get(handles.popup_listsujet, 'string')
filename{id}

DATA = get(handles.GUI_LookMat,'UserData');

matcorr=DATA{id}.MAT;
meancorr=DATA{id}.MAT;
ZoneList=DATA{id}.ZoneList;
save(fullfile(DATA{id}.path,['nan_',DATA{id}.name,'.mat']),'ZoneList','matcorr','meancorr');


% --- Executes on selection change in popupmenu_GraphTheoryOption.
function popupmenu_GraphTheoryOption_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_GraphTheoryOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_GraphTheoryOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_GraphTheoryOption


% --- Executes during object creation, after setting all properties.
function popupmenu_GraphTheoryOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_GraphTheoryOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Connectogram.
function btn_Connectogram_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Connectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%plot the link of the current matrice
%option connectogram or 2d map representation. 

if get(handles.popupmenu_linkoption,'value')==1 %connectogramme 
    % Check zone order for the connectogram. 
    fileorderconnectogramme  = get(handles.edit_linkSettingmapConnectogram,'string')
    fid = fopen( fileorderconnectogramme,'r')
    id = 1
    try
    get(handles.listbox_selectedzone,'string')
    while ~feof(fid)
        tline = fgetl(fid)
        listselected{id} = tline
        id = id+1;
    end
        fclose(fid)
    catch
        disp('Enter a list order for connectogram')
    end
    
  
    
elseif get(handles.popupmenu_linkoption,'value')==3 %2d map 
    
end

DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value')
MAT = DATA{id}.MAT;



List=   strvcat(DATA{id}.ZoneList);
ML =  DATA{id}.zone.ml;
idzone = [];
idlist = [];
idlabelall = [];
%POUR IDENTIFICATION DE TOUTE LES ZONES
for izone = 1:numel( DATA{id}.zone.plotLst)
    chzone =  DATA{id}.zone.plotLst{izone};
    for ichzone = 1:numel(chzone)
        ich = chzone(ichzone);
        strDet = SDDet2strboxy_ISS(ML(ich,2));
        strSrs = SDPairs2strboxy_ISS(ML(ich,1));
        idch = strmatch([strDet, ' ',strSrs ],List,'exact');
     
        
        idlist = [idlist, idch];
        if ichzone==1
            idzone =[idzone, izone];
        else
            idzone =[idzone, 0];
        end
    end
    idlabelall = [idlabelall, {[ DATA{id}.zone.label{izone}]}];
end
  for adji = 1:numel(idlabelall)
            for adjj = 1:numel(idlabelall)   
                labelzone = idlabelall{adji};
                x = strmatch({labelzone} ,idlabelall, 'exact');             
                labelzone = idlabelall{adjj};
                y = strmatch({labelzone} ,idlabelall, 'exact');
                if isempty(x)|isempty(y)
                    msgbox('problem zone in subject')
                end
                chzone = zone.plotLst{x};
                idlisti = [];
                for ichzone = 1:numel(chzone);
                    ich = chzone(ichzone);
                    strDet = SDDet2strboxy_ISS(ML(ich,2));
                    strSrs = SDPairs2strboxy_ISS(ML(ich,1));      
                    idch = strmatch([strDet, ' ',strSrs ],List,'exact');   
                    idlisti = [idlisti, idch];
                end       
                if numel(x)>1
                    x = x(1);
                    msgbox('Attention the zone are duplicated')
                else
                   DATA{id}.zone.chMAT{x} = idlisti;
                end
            end
  end





%%ICI
%     listok = {'Fp1-F7' 'F7-T3' 'T3-T5' 'Fp1-F3' 'Fp2-F8' 'T4-T6' 'Fp2-F4' 'F4-C4' 'C4-P4' 'P4-O2'}
%     listok = {'T3-T5' 'F3-C3' 'C3-P3' 'P3-O1' 'P4-O2' 'C4-P4' 'Fp2-F4' 'T4-T6''F8-T4' 'Fp2-F8' 'Fp1-F7' 'F7-T3' }
%     %zone.label; %ENTER HERE DIFFERENT ORDER{'FT_R', 'TR','TL', 'FT_L','CL','FL', 'FTL','CR'}
%     listok  = {'Fp1-F7' 'Fp1-F3'  'F7-T3' 'T3-T5'  'F3-C3' 'C3-P3' 'P3-O1' 'P4-O2','C4-P4','F4-C4'   'T4-T6' 'F8-T4' 'Fp2-F4'  'Fp2-F8'}
%     listselected = get(handles.listbox_selectedzone,'string');
    try
    listok = listselected;
    catch
        disp('Please enter a list for the connectogram')
    end
    idlist = [];
    idlabel=[];
    idzone =[];
    for ilistzone = 1:numel(listok)
    for izone = 1:numel(DATA{id}.zone.plotLst)
        chzone = DATA{id}.zone.plotLst{izone}
        labelzone = DATA{id}.zone.label{izone}
         x = strmatch({labelzone} , {listok{ilistzone}}, 'exact')
         if ~isempty(x)
%             for ichzone = 1:numel(chzone)
%                 ich = chzone(ichzone);
%                 if strcmp(DATA{id}.System,'ISS')
%                     strDet = SDDet2strboxy_ISS(ML(ich,2));
%                     strSrs = SDPairs2strboxy_ISS(ML(ich,1));
%                     idch = strmatch([strDet, ' ',strSrs ],List,'exact');
%                 end
                 idch=DATA{id}.zone.chMAT{izone};
                 idlist = [idlist, idch];
               % if ichzone==1
               %     idzone =[idzone, izone];
               % else
                    idzone =[idzone,izone, zeros(1,numel(idch)-1)];
              %  end
%             end
            idlabel = [idlabel, {[DATA{id}.zone.label{izone}, sprintf('_%03.0f',ilistzone)]}];
         end
    end
    end
    idline = [find(idzone)-0.5,numel(idzone)+0.5];
%FIN


colorMap = jet(100);
cmin = str2num(get(handles.edit_cmin,'string'))
cmax = str2num(get(handles.edit_cmax,'string'))
cstep = (cmax-cmin)/100;
cf=cmin:cstep:cmax-cstep;
for i=1:size(MAT,1)
    for j = 1:size(MAT,2)
    colorMatrix(i,j) = sum(cf<MAT(i,j));
    end
end



colorlistline = zeros(size(MAT,1),3)

tmp = [find(idzone), numel(idzone)]
tmphald = floor(find(idzone) + (tmp(2:end) - tmp(1:end-1))/2)
idmiddle = zeros(1,numel(idzone))
idmiddle(tmphald ) = idzone(find(idzone));
% Create custom node labels
myLabel = cell(length(x));
%line(numel(listok))
idcolor = 1;
if get(handles.popup_connectogramlabel,'value')==1
for i = 1:length(idzone)
    if idmiddle(i) %idzone(i)
        myLabel{i,1} =  DATA{id}.zone.label{idmiddle(i)};
        idcolor = idcolor + 1
    else
        myLabel{i,1} =  ' '; 
    end
    colorhomemade(i,:) =  colorlistline(idcolor,:)
end

for i = 1:length(idzone)    
    colorhomemade(i,:) =  [0,0,0];
   % idcolor = idcolor + 1
end

elseif get(handles.popup_connectogramlabel,'value')==2
for i=1:numel(idzone)
     myLabel{i,1} = [sprintf('%03.0f', idlist(i)), List(idlist(i),:)]
    
end
end
%
figure;hold on
%%
% Create custom colormap

axes(handles.axes_viewlink)
cla
hold on
x = MAT(idlist,idlist)
id0 = find(isnan(x))
x(id0)=0;
id0 = find(isinf(x))
x(id0)=0;
thresh = str2num(get(handles.edit_threshold,'string'));
x(x >  thresh) = 1;
x(x <= thresh) = 0;
myconnectogram(x,myLabel,colorMatrix,colorMap,idlist)
% myColorMap = jet(length(x));
% myColorMap(:,:)=0

%circularGraph(x,'Label',myLabel);


% --- Executes on button press in btn_2DMAP.
function btn_2DMAP_Callback(hObject, eventdata, handles)
% hObject    handle to btn_2DMAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newmfile = mfilename('fullpath')
[path,name,ext1] = fileparts(newmfile)
pathin = [path,'Template']
 
load(fullfile(pathin,'DoubleBanana_M00.zone'),'-mat')
Loc= load(fullfile(pathin,'1020doubleBananapositionmap.mat'));


listelectrode={'Fp1F7',
    'F7T3',
    'T3T5',    
    'Fp1F3',
    'F3C3',
    'C3P3',
    'Fp2F8',
    'F8T4',
    'T4T6',
    'Fp2F4',
    'F4C4',
    'C4P4'};
id = 1;
matid = zeros(numel(zone.label),numel(zone.label));
    for ielex=1:numel(zone.label)
        ielex
        ieley = 1;
        while ieley <= ielex | ieley==1
            option{id}.ele1 =  listelectrode{ielex};
            option{id}.ele2 =listelectrode{ieley};
            option{id}.matposition = [ielex,ieley];      
            matid(ielex,ieley)=id;
            ieley = ieley + 1;
            id = id + 1;
        end
    end
idhalf = find(matid); %index dans la matrice 99x99
idoption = matid(idhalf); %index dans la matrice option 
figure
map = imread(fullfile(pathin,'DoubleBanana10-20_vr2.png'));
image(map(:,:,1:3));hold on;
set(gca,...
    'XDir','normal',...
    'YDir','reverse');
colormap = jet(100)
valeur = 0:0.01:1
valeur = valeur./3
for io=1:numel(option)
    val = MAT.meancorr(option{io}.matposition(1),option{io}.matposition(2))
    id =sum(val>valeur);
    if id>100
        id=100;
    end
    if val>0
        p1= eval(['Loc.',option{io}.ele1,'.Position']);
        p2= eval(['Loc.',option{io}.ele2,'.Position']);
        plot([p1(1), p2(1)],[p1(2),p2(2)],'color','r','linewidth',2,'displayname', [option{io}.ele1,option{io}.ele2,' ',' idnb ',num2str(io),]);
       % plot([p1(1), p2(1)],[p1(2),p2(2)],'color',colormap(id,:),'marker','x','markersize',12,'linewidth',2,'displayname', [option{io}.ele1,option{io}.ele2,' ',' COH', num2str(  val)]);

    end
    
end


% --- Executes on selection change in popup_template2dmap.
function popup_template2dmap_Callback(hObject, eventdata, handles)
% hObject    handle to popup_template2dmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_template2dmap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_template2dmap


% --- Executes during object creation, after setting all properties.
function popup_template2dmap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_template2dmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_connectogramlabel.
function popup_connectogramlabel_Callback(hObject, eventdata, handles)
% hObject    handle to popup_connectogramlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_connectogramlabel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_connectogramlabel
updateNetAllView(handles)

% --- Executes during object creation, after setting all properties.
function popup_connectogramlabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_connectogramlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_ConnectogramColor.
function popup_ConnectogramColor_Callback(hObject, eventdata, handles)
% hObject    handle to popup_ConnectogramColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_ConnectogramColor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_ConnectogramColor

guidata(handles.GUI_LookMat, handles);



updateNetAllView(handles)
% --- Executes during object creation, after setting all properties.
function popup_ConnectogramColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_ConnectogramColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_loadZoneselected.
function btn_loadZoneselected_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadZoneselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile('.txt')
fid = fopen([path,file],'r')
id = 1
get(handles.listbox_selectedzone,'string')
while ~feof(fid)
    tline = fgetl(fid)
    listok{id} = tline
    id = id+1;
end
fclose(fid)
set(handles.listbox_selectedzone,'string',listok)
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_EntersettingLink.
function btn_EntersettingLink_Callback(hObject, eventdata, handles)
% hObject    handle to btn_EntersettingLink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if get(handles.popupmenu_linkoption,'value')==1
  [file,path] =  uigetfile({'*.xlsx';'*.xls'; '*.txt';'*.*'})
   set(handles.edit_linkSettingmapConnectogram,'string',[path,file]) %connectogramme
elseif get(handles.popupmenu_linkoption,'value')==2
    [file,path] =  uigetfile({ '*.xlsx';'*.xls';'*.txt';'*.*'})
   set(handles.edit_linkSettingmapConnectogram,'string',[path,file])
elseif get(handles.popupmenu_linkoption,'value')==3
   [file,path] =  uigetfile('.mat')
   set(handles.edit_linkSettingmap2dmap,'string',[path,file])
end
    
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)


function edit_linkSettingmapConnectogram_Callback(hObject, eventdata, handles)
% hObject    handle to edit_linkSettingmapConnectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_linkSettingmapConnectogram as text
%        str2double(get(hObject,'String')) returns contents of edit_linkSettingmapConnectogram as a double


% --- Executes during object creation, after setting all properties.
function edit_linkSettingmapConnectogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_linkSettingmapConnectogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_maskoption.
function popupmenu_maskoption_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_maskoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_maskoption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_maskoption
if get(handles.popupmenu_maskoption,'value')==1
    set(handles.edit_filenameMaskMatrice,'visible','off')
    set(handles.btn_BrowseMask,'visible','off')
elseif get(handles.popupmenu_maskoption,'value')==2
    set(handles.edit_filenameMaskMatrice,'visible','on')
    set(handles.btn_BrowseMask,'visible','on')
end
  guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_maskoption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_maskoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_linkoption.
function popupmenu_linkoption_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_linkoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_linkoption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_linkoption

if get(handles.popupmenu_linkoption,'value')==1 %Connectogram
     set(handles.edit_linkSettingmapConnectogram,'visible','on')
    set(handles.edit_linkSettingmap2dmap,'visible','off')
elseif get(handles.popupmenu_linkoption,'value')==2 %notting
    set(handles.edit_linkSettingmapConnectogram,'visible','off')
    set(handles.edit_linkSettingmap2dmap,'visible','off')
elseif get(handles.popupmenu_linkoption,'value')==3
    set(handles.edit_linkSettingmapConnectogram,'visible','off')
    set(handles.edit_linkSettingmap2dmap,'visible','on')   
end %2d map
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_linkoption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_linkoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_linkSettingmap2dmap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_linkSettingmap2dmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_linkSettingmap2dmap as text
%        str2double(get(hObject,'String')) returns contents of edit_linkSettingmap2dmap as a double


% --- Executes during object creation, after setting all properties.
function edit_linkSettingmap2dmap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_linkSettingmap2dmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filenameMaskMatrice_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filenameMaskMatrice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filenameMaskMatrice as text
%        str2double(get(hObject,'String')) returns contents of edit_filenameMaskMatrice as a double


% --- Executes during object creation, after setting all properties.
function edit_filenameMaskMatrice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filenameMaskMatrice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_BrowseMask.
function btn_BrowseMask_Callback(hObject, eventdata, handles)
% hObject    handle to btn_BrowseMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)

% --------------------------------------------------------------------
function menu_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to menu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function menu_colormap_Jet_Callback(hObject, eventdata, handles)
% hObject    handle to menu_colormap_Jet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_axesviewlink_Callback(hObject, eventdata, handles)
% hObject    handle to context_axesviewlink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Make_newfigure_axesviewlink_Callback(hObject, eventdata, handles)
% hObject    handle to Make_newfigure_axesviewlink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = plot_axes_linkConnectogramm(handles,1)

% --------------------------------------------------------------------
function context_link_Callback(hObject, eventdata, handles)
% hObject    handle to context_link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
namech = get(gco,'displayname');
set(handles.context_link_name,'label',namech);

% --------------------------------------------------------------------
function context_link_name_Callback(hObject, eventdata, handles)
% hObject    handle to context_link_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_link_histogramme_Callback(hObject, eventdata, handles)
% hObject    handle to context_link_histogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value');
linkselected = get(handles.context_link_name,'label')
[tok,rem] = strtok(linkselected,'=')
[linkij,rem] = strtok(rem,'=')
[linkname,rem] = strtok(rem,'=')
valG1 = [];
valG2 = [];
groupeall = [];
for i=1:numel(DATA)
    groupeall=[groupeall,DATA{1,i}.GR];
end

figure;hold on
for igr = 1:max(groupeall)
    idgroupe =  find(groupeall==igr);
    valG1 = []
for id = 1:numel(idgroupe)
      idsubjet = idgroupe(id);
    if DATA{idsubjet}.GR == igr
        eval(['new=DATA{',num2str(idsubjet),'}.MAT',linkij,';'])
        if get(handles.radio_fisher,'value')
            new =1/2*(log((1+new )./(1-new )));
        end
        valG1 = [valG1,new];
        eval(['namesubject = DATA{',num2str(idsubjet),'}.name;']);
        plot(igr,new,'x','displayname',[namesubject, '=' num2str(new) ]);
    end

       
%     if DATA{id}.GR == 2
%         eval(['new=DATA{',num2str(id),'}.MAT',linkij])
%         valG2 = [valG2,new];
%           eval(['namesubject = DATA{',num2str(id),'}.name;'])
%         plot(2,new,'x','displayname',[namesubject, '=' num2str(new) ])
%     end
end
 bar(igr,nanmean(valG1),'facealpha',0.5,'displayname',['mean group ',num2str(igr) ]);
end

xlabel('Groups')
   if get(handles.radio_fisher,'value')==0
        ylabel('Connectivity score')
   elseif get(handles.radio_fisher,'value')==1
       ylabel('Connectivity score (Fisher)')
   end
%bar(2,nanmean(valG2),'facealpha',0.5)
[filepath,name,ext] =fileparts(get(handles.edit_subjetxls,'string'))

title([linkname,' ' , name])



% --------------------------------------------------------------------
function context_link_Covariable_Callback(hObject, eventdata, handles)
% hObject    handle to context_link_Covariable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value');
linkselected = get(handles.context_link_name,'label');

covnumber = get(handles.popupmenu_covariable,'value');
covname = get(handles.popupmenu_covariable,'string');
[tok,rem] = strtok(linkselected,'=');
[linkij,rem] = strtok(rem,'=');
[linkname,rem] = strtok(rem,'=');
valG1 = [];
valG2 = [];
groupeall = [];
for i=1:numel(DATA)
    groupeall=[groupeall,DATA{1,i}.GR];
end
colorlink = lines(max(groupeall));

figure;hold on
for id = 1:numel(DATA)
    if DATA{id}.GR > 0
    eval(['COH=DATA{',num2str(id),'}.MAT',linkij,';'])    
    eval(['namesubject = DATA{',num2str(id),'}.name;'])
    eval(['COV=DATA{',num2str(id),'}.cov',num2str(covnumber),';'])   
    hplot= plot(COV,COH,'x','displayname',['G',num2str(DATA{id}.GR),namesubject]);
    set(hplot,'color',colorlink(DATA{id}.GR ,:));
  
    end
end


xlabel(covname{covnumber});
ylabel('COH');
title(linkname);



% --- Executes on selection change in popupmenu_covariable.
function popupmenu_covariable_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_covariable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_covariable contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_covariable


% --- Executes during object creation, after setting all properties.
function popupmenu_covariable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_covariable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function context_mat_copylabel_Callback(hObject, eventdata, handles)
% hObject    handle to context_mat_copylabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1
label = get(handles.Menu_exportSeedFigure,'label')
clipboard('copy', label);



function editchnb_Callback(hObject, eventdata, handles)
% hObject    handle to editchnb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editchnb as text
%        str2double(get(hObject,'String')) returns contents of editchnb as a double


% --- Executes during object creation, after setting all properties.
function editchnb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editchnb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_fisher.
function radio_fisher_Callback(hObject, eventdata, handles)
% hObject    handle to radio_fisher (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_fisher
updateNetAllView(handles)

% --- Executes on button press in radio_negativemap.
function radio_negativemap_Callback(hObject, eventdata, handles)
% hObject    handle to radio_negativemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_negativemap
updateNetAllView(handles)


% --- Executes on button press in btn_savelisttxt.
function btn_savelisttxt_Callback(hObject, eventdata, handles)
% hObject    handle to btn_savelisttxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uiputfile('.txt')
fid = fopen([path,file],'w')
tmp = get(handles.listbox_selectedzone,'string')
for i=1:size(tmp,1)
    fprintf(fid,'%s\n',tmp{i})
end
fclose(fid)



% --- Executes on button press in btn_autoscale.
function btn_autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to btn_autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
%data = get(handles.GUI_LookMat,'UserData')
guidata(handles.GUI_LookMat, handles);
handles.GUI_LookMat

DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value')
MAT = DATA{id}.MAT

cmax = max(MAT)
cmin = min(MAT)
cmin = sprintf('%0.2g',cmin);
cmax = sprintf('%0.2g',cmax);
set(handles.edit_cmin,'string',[cmin]);
set(handles.edit_cmax,'string',[cmax]);

guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles)
 

msgbox('Not code yet !')


