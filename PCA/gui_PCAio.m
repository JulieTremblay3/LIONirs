function varargout = gui_PCAio(varargin)
% GUI_PCAIO M-file for gui_PCAio.fig
%      GUI_PCAIO, by itself, creates a new GUI_PCAIO or raises the existing
%      singleton*.
%
%      H = GUI_PCAIO returns the handle to a new GUI_PCAIO or the handle to
%      the existing singleton*.
%
%      GUI_PCAIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PCAIO.M with the given input arguments.
%
%      GUI_PCAIO('Property','Value',...) creates a new GUI_PCAIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_PCAio_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_PCAio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_PCAio

% Last Modified by GUIDE v2.5 28-Sep-2018 15:46:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_PCAio_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_PCAio_OutputFcn, ...
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


% --- Executes just before gui_PCAio is made visible.
function gui_PCAio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_PCAio (see VARARGIN)

% Choose default command line output for gui_PCAio
handles.output = hObject;


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

indt = tstart:tstop;
intensnorm = d(indt,:);
  
%Detrent DATA segment for centrering 
X = 1:1:size(intensnorm,1);
Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
Mb2 =  intensnorm(1,:)'; %offset    
A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
spar = intensnorm - A;
idnan = find(isnan(sum(spar,1)))
listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
MeasListActplotLst = zeros(size(A,2),1);
MeasListActplotLst(PMI{1}.plotLst,:)=1;
if min(PMI{1}.plotLst-size(A,2)/2) >0
    MeasListActplotLst(PMI{1}.plotLst-size(A,2)/2,:)=1;
end
MeasListActplotLst((size(A,2)/2+1):end,1)=0; %forcer la fin a zeros
 
listgood = PMI{currentsub}.plotLst ; % find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));
nc830 = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
all830 = zeros(numel(PMI{currentsub}.data(cf).MeasListAct),1);
all830(listgood)=1;
tmp = reshape(all830,nc830 ,2);
t = sum(tmp');
listgood = find([t,t]);
PMI{currentsub}.plotLst = listgood;

d1 = spar(:,listgood);

t = PMI{currentsub}.data(cf).HRF.tHRF;

if 0
    AdvOptions.FiltOptions.FilterType=1;  %Default Butterworth
    AdvOptions.FiltOptions.FilterOrder=3;  %Default 3rd order
    lpf = 0.2;
    fs = 1/(t(2)-t(1)); 
    [fb,fa] = MakeFilter(AdvOptions.FiltOptions.FilterType,AdvOptions.FiltOptions.FilterOrder,fs,lpf,'low');
    d1 = filtfilt(fb,fa,d1);    
end

c = d1'*d1;
%figure;imagesc(c)
[v,s,foo]=svd(c);
svs = diag(s);
u = d1*v*inv(s);

%CHECK COV
% figure
% imagesc(d1*d1')
% title('iner product time' )
% figure
% imagesc(d1'*d1)
% title('iner product channel' )
%end


PMI{currentsub}.tmpPCA.selected =1;
PMI{currentsub}.tmpPCA.d1 = d1;                    %Data initial
PMI{currentsub}.tmpPCA.checksumd = sum(intensnorm(:)); 
PMI{currentsub}.tmpPCA.listgood = listgood ; 
PMI{currentsub}.tmpPCA.indt = [tstart(end),tstop(end)];%Time indice
PMI{currentsub}.tmpPCA.v = v; %eigen vector
PMI{currentsub}.tmpPCA.svs = svs;  %singular value
PMI{currentsub}.tmpPCA.s = s;
PMI{currentsub}.tmpPCA.u = u; %PCA time x space

lstSV =1; %composante 1
temp = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';

% figure
% plot(temp)
% plot(u(:,1))
% plot(u(:,2))
% 
% %'time')
% figure;hold on
% plot(u(1,:))
% plot(u(2,:))
% plot(u(3,:))
% 
% %'spatial')
% figure; hold on;
% plot(u(:,1))
% plot(u(:,2))
% plot(u(:,3))
% title('spatial')
% 
% %'% variance explain')
% svs(1)/sum(svs);
% svs(2)/sum(svs);
% svs(3)/sum(svs);

% Norm_Cov = d1- u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)'+1;
% handles.dataPCA(:,plotLst) = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)'-1;%Norm_Cov;
id = find(svs/sum(svs)>0.01);
label = [];
for i=1:numel(id)
    label = [label,{['F',num2str(i), 'var=', num2str(svs(i)/sum(svs))]}];
end

set(guiHOMER,'UserData',PMI);


nbcomp = length(PMI{currentsub}.plotLst); % Get the number of channels used.
if nbcomp > 10; nbcomp = 10; end % If there is too many channel, just display 10 components.
label = [];
for i=1:nbcomp
    label = [label,{['F',num2str(i)]}]; % Create the channels in the dropdown listbox.
end

set(handles.popup_Component,'string',label)
resetview(handles)

guidata(hObject, handles);




% Update handles structure

% UIWAIT makes gui_PCAio wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_PCAio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popup_Component.
function popup_Component_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Component contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Component
resetview(handles)

% --- Executes during object creation, after setting all properties.
function popup_Component_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_viewtopo.
function btn_viewtopo_Callback(hObject, eventdata, handles)
% hObject    handle to btn_viewtopo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
option.conc = 1;
option.SPMgui = 1;
IO_HelmetMTG_Display(handles,option)

function resetview(handles)

newfigure  = 0;
axes(handles.axes_PCA);hold on
cla;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
u =PMI{currentsub}.tmpPCA.u ;
s =PMI{currentsub}.tmpPCA.s ;
v =PMI{currentsub}.tmpPCA.v ;
lstSV = get(handles.popup_Component,'value');
if size(u,2) >= lstSV % Check if the component index is valid (inside array bounds).
    temp = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
    id = 1:size(temp,1);
    t=PMI{currentsub}.data(1).HRF.tHRF(id);
    plot(t,temp);
else
    errordlg('Not enough selected channels.','Operation error') % Inform the user of the error.
end


PMI{currentsub}.tmpPCA.selected =get(handles.popup_Component,'value');
set(guiHOMER,'UserData',PMI);