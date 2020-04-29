function varargout = GUI_GLMidentification(varargin)
% GUI_GLMIDENTIFICATION M-file for GUI_GLMidentification.fig
%      GUI_GLMIDENTIFICATION, by itself, creates a new GUI_GLMIDENTIFICATION or raises the existing
%      singleton*.
%
%      H = GUI_GLMIDENTIFICATION returns the handle to a new GUI_GLMIDENTIFICATION or the handle to
%      the existing singleton*.
%
%      GUI_GLMIDENTIFICATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GLMIDENTIFICATION.M with the given input arguments.
%
%      GUI_GLMIDENTIFICATION('Property','Value',...) creates a new GUI_GLMIDENTIFICATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_GLMidentification_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_GLMidentification_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_GLMidentification

% Last Modified by GUIDE v2.5 10-Jan-2019 15:38:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_GLMidentification_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_GLMidentification_OutputFcn, ...
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


% --- Executes just before GUI_GLMidentification is made visible.
function GUI_GLMidentification_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_GLMidentification (see VARARGIN)

% Choose default command line output for GUI_GLMidentification
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
 set(handles.text_tstart,'string', ['Time :' num2str(PMI{currentsub}.data(cf).HRF.tHRF(tstart)), ' to ',num2str(PMI{currentsub}.data(cf).HRF.tHRF(tstop))] ) 
%Detrent DATA segment for centrering 
X = 1:1:size(intensnorm,1);
sous = intensnorm(end,:) - intensnorm(1,:)
Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
Mb2 =  intensnorm(1,:)'; %offset    
A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
spar = intensnorm - A; %with detrending 
spar = intensnorm;     %without detrending
listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
MeasListActplotLst = zeros(size(A,2),1);
MeasListActplotLst(PMI{1}.plotLst,:)=1;
if min(PMI{1}.plotLst-size(A,2)/2) >0
    MeasListActplotLst(PMI{1}.plotLst-size(A,2)/2,:)=1;
end
MeasListActplotLst((size(A,2)/2+1):end,1)=0; %forcer la fin a zeros
listgood = find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));

isfield( PMI{currentsub}.data(cf),'AUX') 
PMI{currentsub}.data(cf).AUX

t = PMI{currentsub}.data(cf).HRF.tHRF;
% spar = cat(4,spar(:,:,1:end/2),spar(:,:,end/2+1:end));
% spartmp = spar(:,:,listgood,:);
PMI{1}.tmpGLM.selected =1     %used to visualised
PMI{1}.tmpGLM.spar = spar;                    %Data initial
PMI{1}.tmpGLM.listgood = listgood ; 
PMI{1}.tmpGLM.indt = [tstart(end),tstop(end)];%Time indice
load(PMI{1}.NIRSmat{1});
iAUX = 1 %each aux file could have a different sample rate
fsNIRS=NIRS.Cf.dev.fs;
idfile=PMI{1}.idfile;
iRegressor =  2;
if isfield(PMI{currentsub}.tmpGLM,'AUX')
    rmfield(PMI{currentsub}.tmpGLM,'AUX')
end

for iAUX = 1:numel(NIRS.Dt.AUX)
    nameAUX = NIRS.Dt.AUX(iAUX).pp(end).p{idfile};
    if isfield(NIRS.Dt.AUX(iAUX).pp(end),'sync_timesec')
        tstartf = NIRS.Dt.AUX(iAUX).pp(end).sync_timesec{idfile};
    else 
        tstartf = 1;
    end
    tstopf = tstartf+PMI{currentsub}.data(cf).HRF.tHRF(end);
    [pathtmp,filetmp,exttmp]=fileparts(nameAUX);
    [data,infoBV,label,ind_dur_ch] = fopen_EEG(nameAUX, tstartf, tstopf);
      for ich=1:numel(infoBV.name_ele)
            PMI{currentsub}.tmpGLM.AUX.label{iRegressor} =[filetmp,' ',infoBV.name_ele{ich}];           
            fsAUX =1/(infoBV.SamplingInterval/1000000); %Frequence echantillonage Hz
            tmp = data(:,ich);
            [p,q] = rat(fsNIRS/fsAUX,0.0001);
            tmpr=resample( tmp , p, q);
            % we resample aux to the data initial size
            if numel(tmpr)<numel(t)
                nplus = numel(t)-numel(tmpr)
                tmpr = [tmpr;tmpr(end-nplus:end) ]
            elseif numel(tmpr)>numel(t)
                tmpr = tmpr(1:numel(t));
            end
            %we cut to fit the spar selection                      
            PMI{currentsub}.tmpGLM.AUX.fs{iRegressor} = fsNIRS;
            PMI{currentsub}.tmpGLM.AUX.data{iRegressor} = tmpr(PMI{1}.tmpGLM.indt(1):PMI{1}.tmpGLM.indt(end));
            clear tmpr
            PMI{currentsub}.tmpGLM.AUX.view{iRegressor} = 1;
            iRegressor = iRegressor +1;
        end
end

PMI{1}.tmpGLM.AUX.label{1} = 'Constant'
PMI{1}.tmpGLM.AUX.fs{1} = fsNIRS;
PMI{1}.tmpGLM.AUX.data{1} = ones(numel(PMI{1}.tmpGLM.indt(1):PMI{1}.tmpGLM.indt(end) ),1);
PMI{1}.tmpGLM.AUX.view{1} = 1;

        set(handles.listbox_auxilary,'string',PMI{currentsub}.tmpGLM.AUX.label');
        set(handles.listbox_auxilary,'enable','on');
%load and resample the auxiliary for use in the GLM 

% % handles.NIRS.Dt.AUX(iAUX).pp(end).p{idfile}
%     if isfield(handles.NIRS.Dt,'AUX') %LOAD AUX ICI AUXILIARY AND ADD IN PMI structure
%         iAUX = get(handles.popupauxnb,'value')
%         PMI{currentsub}.data(cf).AUX = [];
%         nameAUX = handles.NIRS.Dt.AUX(iAUX).pp(end).p{idfile};
%         [data,infoBV,label,ind_dur_ch] = fopen_EEG(nameAUX);
%         for ich=1:numel(infoBV.name_ele)
%             PMI{currentsub}.data(cf).AUX.label{ich} =[infoBV.name_ele{ich},'  on'];
%             PMI{currentsub}.data(cf).AUX.data{ich} = data(:,ich);
%             PMI{currentsub}.data(cf).AUX.fs{ich} =1/(infoBV.SamplingInterval/1000000); %Frequence echantillonage Hz
%             PMI{currentsub}.data(cf).AUX.view{ich} = 1;
%         end
%         set(handles.listbox_AUXCH,'string',PMI{currentsub}.data(cf).AUX.label')
%         set(handles.listbox_AUXCH,'enable','on')
%     end
%     
% 
% load

PMI{currentsub}.data(cf).AUX
set(guiHOMER,'UserData',PMI);
% btn_MultipleLinearRegression_Callback(hObject, eventdata, handles)
% nbcomp = str2num(get(handles.edit_nbComponent,'string'));
% label = [];
% for i=1:nbcomp
%     label = [label,{['F',num2str(i)]}];
% end
% set(handles.popupmenu_Factor,'string',label)



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_GLMidentification wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_GLMidentification_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_listAuxiliary_Callback(hObject, eventdata, handles)
% hObject    handle to edit_listAuxiliary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_listAuxiliary as text
%        str2double(get(hObject,'String')) returns contents of edit_listAuxiliary as a double


% --- Executes during object creation, after setting all properties.
function edit_listAuxiliary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_listAuxiliary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_MultipleLinearRegression.
function btn_MultipleLinearRegression_Callback(hObject, eventdata, handles)
% hObject    handle to btn_MultipleLinearRegression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1
handlesoutput = hObject;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
y=PMI{1}.tmpGLM.spar; 
Regressorlist = get(handles.listbox_AUXselected,'string')
if ~iscell(Regressorlist)
    Regressorlist = {Regressorlist}
end
%find indice of regressor in the list
idreg = [];
for iReg = 1:numel(Regressorlist)
    for ilist = 1:numel(PMI{1}.tmpGLM.AUX.label)
        if strcmp(PMI{1}.tmpGLM.AUX.label{ilist},Regressorlist{iReg})
        idreg = [idreg,ilist];
        end
    end
    
end
PMI{1}.tmpGLM.idreg = idreg; %regression auxiliary indice used 
idreg
X = [];%ones(size(PMI{1}.tmpGLM.spar,1),1);
for ireg = 1:numel(idreg)
    X = [X,    PMI{1}.tmpGLM.AUX.data{idreg(ireg)}];
end
%figure;plot(X)
%figure;plot(PMI{1}.tmpGLM.AUX.data{idreg(ireg)})
beta = zeros(size(X,2),numel(PMI{1}.tmpGLM.listgood))
bstd = zeros(size(X,2),numel(PMI{1}.tmpGLM.listgood))
for ich = 1:numel(PMI{1}.tmpGLM.listgood)
    idch = PMI{1}.tmpGLM.listgood(ich);
    y = PMI{1}.tmpGLM.spar(:,idch);
   %figure;plot(PMI{1}.tmpGLM.spar)
    if sum(isnan(y)) 
        b = zeros(size(X,2),1)
        beta(:,ich) = 0;
        bstd(:,ich) = 0;

    else
        [b,bint,r,rint,stats]=  regress(y,X);
        beta(:,ich) = b;
        bstd(:,ich) = stats(4);
    end
end
PMI{1}.tmpGLM.beta = beta;
PMI{1}.tmpGLM.std = bstd;
% PMI{1}.tmpGLM.selected = 1;
% % PMI{1}.tmpGLM.Xm = 
% % PMI{1}.tmpGLM.data = 
% set(handles.listbox_AUXselected,'value',1) %PMI{1}.tmpGLM.selected );
set(guiHOMER,'UserData',PMI);

% --- Executes on selection change in popup_factor.
function popup_factor_Callback(hObject, eventdata, handles)
% hObject    handle to popup_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_factor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_factor


% --- Executes during object creation, after setting all properties.
function popup_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_auxilary.
function listbox_auxilary_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_auxilary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_auxilary contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_auxilary


% --- Executes during object creation, after setting all properties.
function listbox_auxilary_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_auxilary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_get_regressor.
function btn_get_regressor_Callback(hObject, eventdata, handles)
% hObject    handle to btn_get_regressor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listchoose =get(handles.listbox_AUXselected,'string');
if isempty(listchoose)
    listchoose = [];
elseif  ~iscell(listchoose)
    listchoose = {listchoose};
end
listaux = get(handles.listbox_auxilary,'string');
valaux = get(handles.listbox_auxilary,'value'); 
listchoose = [listchoose;listaux{valaux}]
set(handles.listbox_AUXselected,'string',listchoose)

% --- Executes on button press in btn_removeregressor.
function btn_removeregressor_Callback(hObject, eventdata, handles)
% hObject    handle to btn_removeregressor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listchoose = get(handles.listbox_AUXselected,'string');
valchoose = get(handles.listbox_AUXselected,'value'); 
listchoose(valchoose) = [];
set(handles.listbox_AUXselected,'string',listchoose);
set(handles.listbox_AUXselected,'value',1); 

% --- Executes on selection change in listbox_AUXselected.
function listbox_AUXselected_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_AUXselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_AUXselected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_AUXselected
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
PMI{1}.tmpGLM.selected = get(handles.listbox_AUXselected,'value');
set(guiHOMER,'UserData',PMI);


% --- Executes during object creation, after setting all properties.
function listbox_AUXselected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_AUXselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
