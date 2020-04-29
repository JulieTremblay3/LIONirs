function varargout = gui_icaIO(varargin)
% GUI_ICAIO M-file for gui_icaIO.fig
%      GUI_ICAIO, by itself, creates a new GUI_ICAIO or raises the existing
%      singleton*.
%
%      H = GUI_ICAIO returns the handle to a new GUI_ICAIO or the handle to
%      the existing singleton*.
%
%      GUI_ICAIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ICAIO.M with the given input arguments.
%
%      GUI_ICAIO('Property','Value',...) creates a new GUI_ICAIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_icaIO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_icaIO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_icaIO

% Last Modified by GUIDE v2.5 16-Oct-2017 14:59:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_icaIO_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_icaIO_OutputFcn, ...
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


% --- Executes just before gui_icaIO is made visible.
function gui_icaIO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_icaIO (see VARARGIN)

% Choose default command line output for gui_icaIO
%CALCUL INTERACTIF PARAFAC SELECTIONNNER POUR LA FENETRE DE TEMPS DANS LE
%GUI
handles.output = hObject;
tstart = varargin{1};
tstop = varargin{2};
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
d = PMI{currentsub}.data(cf).HRF.AvgC;
indt = tstart(end):tstop(end);
intensnorm = d(indt,:);
  
%Detrent DATA segment for centrering 
X = 1:1:size(intensnorm,1);
Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
Mb2 =  intensnorm(1,:)'; %offset    
A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
spar = intensnorm - A;

listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
MeasListActplotLst = zeros(size(A,2),1);
MeasListActplotLst(PMI{1}.plotLst,:)=1;
if min(PMI{1}.plotLst-size(A,2)/2) >0
    MeasListActplotLst(PMI{1}.plotLst-size(A,2)/2,:)=1;
end 
MeasListActplotLst((size(A,2)/2+1):end,1)=0; %forcer la fin a zeros
nc830 = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
listgood = PMI{currentsub}.plotLst;
%ensure both wavelength are selected
nc830 = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
all830 = zeros(numel(PMI{currentsub}.data(cf).MeasListAct),1);
all830(listgood)=1;
tmp = reshape(all830,nc830 ,2);
t = sum(tmp');
listgood = find([t,t]);
PMI{currentsub}.plotLst = listgood;

% PMI{currentsub}.plotLst

t = PMI{currentsub}.data(cf).HRF.tHRF;
%spar = cat(3,spar(:,1:end/2),spar(:,end/2+1:end));
spartmp = spar(:,listgood);
PMI{currentsub}.tmpICA.spar = [];
PMI{currentsub}.tmpICA.listgood =[];
PMI{currentsub}.tmpICA.indt =[];
PMI{1}.tmpICA.selected =1;     
PMI{currentsub}.tmpICA.spar = spartmp;  %cat(2, spartmp(:,:,1),spartmp(:,:,2));                    %Data initial
PMI{currentsub}.tmpICA.listgood = listgood ; 
PMI{currentsub}.tmpICA.indt = [tstart(end),tstop(end)];%Time indice
set(guiHOMER,'UserData',PMI);
btn_ICA_Callback(hObject, eventdata, handles);
set(handles.edit_nbComponent,'string',num2str(numel(listgood)*2));
nbcomp = str2num(get(handles.edit_nbComponent,'string'));
label = [];
for i=1:nbcomp
    label = [label,{['F',num2str(i)]}];
end
set(handles.popupmenu_ICA,'string',label);
set(handles.popupmenu_ICA,'value',1);
resetview(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_ParafacIO wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% UIWAIT makes gui_icaIO wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_icaIO_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_ICA.
function popupmenu_ICA_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_ICA contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ICA
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
PMI{1}.tmpICA.selected =get(handles.popupmenu_ICA,'value')
set(guiHOMER,'UserData',PMI);
resetview(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_ICA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_ICA.
function btn_ICA_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
DataCanauxSupp = PMI{currentsub}.tmpICA.spar;
Nc = str2num(get(handles.edit_nbComponent,'string'));
% figure 
% plot(spar(:,:,1))
% find(isnan(DataCanauxSupp))
%  [U,s,V] = csvd(DataCanauxSupp');
%  %donot include svd with s=0
% %  figure; plot(DataCanauxSupp)
% %  
% %  idremove = find(s<0.000001)
% %  U(idremove,:) = [];
% %  U(:,idremove) = [];
% sa=s*ones(1,size(V,1));
% % V(:,idremove)=[];
% % 
% for SVDnumber = 1:size(DataCanauxSupp,2)-1
%     try
%         [Ua,rho(SVDnumber),eta(SVDnumber)] = tsvdregtool(U,sa,V,DataCanauxSupp',SVDnumber);
%     catch
%         rho(SVDnumber)=nan;
%         eta(SVDnumber)=nan; 
%     end
% end
%      [reg_corner,rho,eta,reg_param] = l_curve(U,sa,Ua,'tsvd');
%     figure
%     plot(rho)
%     plot(eta)
%   %  [Ua,rho(SVDnumber),eta(SVDnumber)] =
%   %  tsvdregtool(U,s,V,dUa,reg_corner);
% Nc=size(DataCanauxSupp',1)
try
[icasig,A, W] = fastica(DataCanauxSupp',10,'numOfIC' , Nc ,'g', 'gauss','approach', 'symm','verbose','on','only', 'ICA');
catch
    close(gcf)
    msgbox(['Please install FASTICA toolbox for matlab https://research.ics.aalto.fi/ica/fastica/ ',...
    'to perform ICA decomposition Hyvärinen et al. 2000 DOI: 10.1016/S0893-6080(00)00026-5'])       
    return
end
%         figure;plot(icasig')
%         figure;plot(mean(DataCanauxSupp))
PMI{currentsub}.tmpICA.Factors{1}=icasig;
PMI{currentsub}.tmpICA.Factors{2}=W;
PMI{currentsub}.tmpICA.Factors{3}=A;

set(guiHOMER,'UserData',PMI);
guidata(hObject, handles);
resetview(handles)
%set(handles.text_ConcordiaVal,'string', num2str(corcondia))
nbcomp = str2num(get(handles.edit_nbComponent,'string'));
label = [];
for i=1:nbcomp
    label = [label,{['F',num2str(i)]}];
end
set(handles.popupmenu_ICA,'string',label);
set(handles.popupmenu_ICA,'value',1);

function edit_nbComponent_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nbComponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nbComponent as text
%        str2double(get(hObject,'String')) returns contents of edit_nbComponent as a double

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
PMI{1}.tmpICA.selected =get(handles.popupmenu_ICA,'value');
set(guiHOMER,'UserData',PMI);
btn_ICA_Callback(hObject, eventdata, handles);
resetview(handles);

% --- Executes during object creation, after setting all properties.
function edit_nbComponent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nbComponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function resetview(handles)
%All normal 
axes(handles.axes_timecourse);
cla
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');


    FacTime = PMI{currentsub}.tmpICA.Factors{1}';
    FacSpatial = PMI{currentsub}.tmpICA.Factors{2};
    Nc = size(FacTime,2);
    colorlist = jet(size(FacTime,2));
    cf = PMI{currentsub}.currentFile;

    axes(handles.axes_timecourse);
    cla
    hold on
    t = PMI{currentsub}.data(cf).HRF.tHRF(PMI{currentsub}.tmpICA.indt(1):PMI{currentsub}.tmpICA.indt(2));

    for i=1:Nc
         if PMI{1}.tmpICA.selected==i
            plot(t, FacTime(:,i),'color',colorlist(i,:),'linewidth',4 )
         else
            plot(t, FacTime(:,i),'color',colorlist(i,:) )
         end
    end

 
    axes(handles.axes_Topo);
    cla
    hold on
    for i=1:Nc
        if PMI{1}.tmpICA.selected==i
            plot(FacSpatial(i,:),'color',colorlist(i,:),'linewidth',4 )
        else
            plot(FacSpatial(i,:),'color',colorlist(i,:) )
        end
    end

    
