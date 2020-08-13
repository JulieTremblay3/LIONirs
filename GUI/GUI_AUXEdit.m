function varargout = GUI_AUXEdit(varargin)
% GUI_AUXEDIT M-file for GUI_AUXEdit.fig
%      GUI_AUXEDIT, by itself, creates a new GUI_AUXEDIT or raises the existing
%      singleton*.
%
%      H = GUI_AUXEDIT returns the handle to a new GUI_AUXEDIT or the handle to
%      the existing singleton*.
%
%      GUI_AUXEDIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_AUXEDIT.M with the given input arguments.
%
%      GUI_AUXEDIT('Property','Value',...) creates a new GUI_AUXEDIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_AUXEdit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_AUXEdit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_AUXEdit

% Last Modified by GUIDE v2.5 09-Sep-2019 15:03:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_AUXEdit_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_AUXEdit_OutputFcn, ...
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


% --- Executes just before GUI_AUXEdit is made visible.
function GUI_AUXEdit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_AUXEdit (see VARARGIN)

% Choose default command line output for GUI_AUXEdit
handles.output = hObject;
popupmenuoption_Callback(hObject, eventdata, handles)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_AUXEdit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_AUXEdit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_rawdata_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rawdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rawdata as text
%        str2double(get(hObject,'String')) returns contents of edit_rawdata as a double





% --- Executes on button press in btn_Browserawdata.
function btn_Browserawdata_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Browserawdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listoption = get(handles.popupmenuoption,'string');
idval = get(handles.popupmenuoption,'value');
if strcmp(listoption{idval},'Create HRF using onset duration xls file')
 [file,pathout]= uigetfile({'*.nir';'*.dat';'*.*';},'Selector file to create analogous trig in aux file .nir or .dat EEG BrainVision generic data export');
 set(handles.edit_rawdata,'string', [pathout,file])
elseif strcmp(listoption{idval},'Global regressor NIRS.mat average exclude rejected channel')
 [file,pathout]= uigetfile({'NIRS.mat'},'Find global regressor form NIRS.mat data, use before normalisation and segmentation');
 set(handles.edit_rawdata,'string', [pathout,file])
elseif strcmp(listoption{idval},'Short channel regressor NIRS.mat (channel is measurement in the data)')
 [file,pathout]= uigetfile({'NIRS.mat'},'Find global regressor form NIRS.mat data, use before normalisation and segmentation');
 set(handles.edit_rawdata,'string', [pathout,file])
elseif strcmp(listoption{idval},'Audio .wav (sum mel coefficients)')
 [file,pathout]= uigetfile({'*.nir';'*.dat';'*.*';},'Selector file to create analogous trig in aux file .nir or .dat EEG BrainVision generic data export');
  set(handles.edit_rawdata,'string', [pathout,file]);
elseif strcmp(listoption{idval},'Concatenate AUX/EEG')
 [file,pathout]= uigetfile({'*.dat'},'Selector multiple EEG or AUX create by BrainVision generic data export','MultiSelect','on');
 set(handles.edit_rawdata,'string', [pathout]);
 set(handles.listbox_AUXCONCATENATE,'string',file);
elseif strcmp(listoption{idval},'VIDEO CONVERT and MERGE file CODEC 32 to 64 bit')
 [file,pathout]= uigetfile({'*.avi'},'Select video to convert','MultiSelect','on');
  set(handles.edit_rawdata,'string', [pathout]);
  set(handles.listbox_AUXCONCATENATE,'string',file);
  
end


function edit_eventfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eventfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_eventfile as text
%        str2double(get(hObject,'String')) returns contents of edit_eventfile as a double


% --- Executes during object creation, after setting all properties.
function edit_eventfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eventfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Browseeventfile.
function btn_Browseeventfile_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Browseeventfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listoption = get(handles.popupmenuoption,'string');
idval = get(handles.popupmenuoption,'value');
if strcmp(listoption{idval},'Create HRF using onset duration xls file')
    [file,path]=uigetfile('.xlsx');
elseif strcmp(listoption{idval},'Create HRF convolution .dat AUX signal ')%
    [file,path]=uigetfile('.dat'); 
elseif strcmp(listoption{idval},'Global regressor NIRS.mat average exclude rejected channel')
    path = uigetdir();
    file = [];
elseif strcmp(listoption{idval},'Short channel regressor NIRS.mat (channel is measurement in the data)')
    path = uigetdir();
    file = [];
elseif strcmp(listoption{idval},'Audio .wav (sum mel coefficients)')
    [file,path] = uigetfile('.wav');
end
set(handles.edit_eventfile,'string',[path,file])
% --- Executes on button press in btn_createHRF.
function btn_createHRF_Callback(hObject, eventdata, handles)
% hObject    handle to btn_createHRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileevent = get(handles.edit_eventfile,'string');
filedata = get(handles.edit_rawdata,'string');
label = get(handles.edit_name,'string');
pathout = get(handles.edit_pathout,'string');
%Read onset
listoption = get(handles.popupmenuoption,'string');
idval = get(handles.popupmenuoption,'value');

if strcmp(listoption{idval},'Create HRF using onset duration xls file')
    [~,~,ext] =fileparts(fileevent);
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        [num, txt, raw] = xlsread(fileevent);
    elseif strcmp(ext,'.txt')   
        [num, txt, raw] = readtxtfile_asxlsread(fileevent);
    end
    [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]= fopen_EEG(filedata);
    idtrigger = find(strcmp('trigger',EEG.marker(:,1))& EEG.ind_dur_ch(:,1)>0);
    EEG.marker = EEG.marker(idtrigger,:);
    EEG.ind_dur_ch = EEG.ind_dur_ch(idtrigger,:);
    %HRF convolution
    onset =num(:,1);
    dur = num(:,2);
    weight = num(:,3);
    eventV =zeros(size(EEG.data,1),1);
    fs = 1/(EEG.infoBV.SamplingInterval/1e6)
    t = 1/fs:1/fs:1/fs*size(EEG.data,1);
    idevent = [];
    for ievent = 1:numel(onset)
        idon = find(t>onset(ievent) & t<(onset(ievent)+dur(ievent)));      
        eventV(idon) = weight(ievent,1);
    end
    %eventV( idevent) = 1;
   
    frametimes = 1/fs:1/fs:size(EEG.data,1)*1/fs;
    param1 = str2num(get(handles.edit_parameter_1,'string'));
    param2 = str2num(get(handles.edit_parameter_2,'string'));
    param3 = str2num(get(handles.edit_parameter_3,'string'));
    param4 = str2num(get(handles.edit_parameter_4,'string'));
    param5 = str2num(get(handles.edit_parameter_5,'string'));
    H = fmridesign(frametimes,0,[1 0],[],[param1 param2 param3 param4 param5 0]); %Default
    Sc1 = conv2(H.X(:,1,1), eventV);
    Scn = Sc1(1:end-numel(H.X(:,1,1))+1);

     Scn = Scn/max(Scn);
     Regressor.xdata = t'
     Regressor.ydata = Scn
     Onset.xdata = [onset,onset];
     Onset.ydata = [min(Scn),max(Scn)];
    save(fullfile(pathout,'Regressor.mat'),'Regressor', 'Onset')
    h= figure;plot(t,Scn,'displayname','HRF');
     hold on
     plot([onset,onset],[min(Scn),max(Scn) ],'displayname','onset' );
    filename = 'AUX';
    saveas(h,fullfile(pathout,'Regressor.fig'))
    
    
    
    HRFcurve.data = Scn;
    HRFcurve.ind_dur_ch=EEG.ind_dur_ch;
    HRFcurve.marker=EEG.marker;
    HRFcurve.infoBV = EEG.infoBV;
    HRFcurve.infoBV.DataType = 'TIMEDOMAIN';
    HRFcurve.infoBV.NumberOfChannels = '1';
    HRFcurve.infoBV.name_ele = {deblank(label)};
    HRFcurve.infoBV.coor_r = 1;
    HRFcurve.infoBV.coor_theta  = 90;
    HRFcurve.infoBV.coor_phi = 45;
    HRFcurve.infoBV.DataPoints = numel(Scn);
    [pathstr, name, ext] = fileparts(fileevent);
    fileoutHRFcurve = fullfile(pathout, [label,'.dat']);
    fwrite_EEG(fileoutHRFcurve,HRFcurve,1,HRFcurve.infoBV.DataPoints );
elseif strcmp(listoption{idval},'Create HRF convolution .dat AUX signal ')
    [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]= fopen_EEG(fileevent);
    idtrigger = find(strcmp('trigger',EEG.marker(:,1))& EEG.ind_dur_ch(:,1)>0);
    EEG.marker = EEG.marker(idtrigger,:);
    EEG.ind_dur_ch = EEG.ind_dur_ch(idtrigger,:);
     
    fs = 1/(EEG.infoBV.SamplingInterval/1e6);
    t = 1/fs:1/fs:1/fs*size(EEG.data,1);
    frametimes = 1/fs:1/fs:30;
    
    
    param1 = str2num(get(handles.edit_parameter_1,'string'));
    param2 = str2num(get(handles.edit_parameter_2,'string'));
    param3 = str2num(get(handles.edit_parameter_3,'string'));
    param4 = str2num(get(handles.edit_parameter_4,'string'));
    param5 = str2num(get(handles.edit_parameter_5,'string'));
    param6 = str2num(get(handles.edit_parameter_6,'string'));
     
    figure;
    H = fmridesign(frametimes,0,[1 0],[],[param1 param2 param3 param4 param5 0]); %Default
    if param6 %binarise event
        tr = max(EEG.data)*0.25;
        eventV = EEG.data>tr;
        subplot(2,1,1)
        plot(t,eventV);title('EVENT binarize');
    else
        eventV = EEG.data;      
        subplot(2,1,1)
        plot(t,eventV);title('EVENT analog');
    end
    Sc1 = conv2(H.X(:,1,1), eventV);
    Scn = Sc1(1:end-numel(H.X(:,1,1))+1);
    subplot(2,1,2);plot(t, Scn);title('HRF');
    filename = 'AUX';
    HRFcurve.data = Scn;
    HRFcurve.ind_dur_ch=EEG.ind_dur_ch;
    HRFcurve.marker=EEG.marker;
    HRFcurve.infoBV = EEG.infoBV;
    HRFcurve.infoBV.DataType = 'TIMEDOMAIN';
    HRFcurve.infoBV.NumberOfChannels = '1';
    HRFcurve.infoBV.name_ele = {deblank(label)};
    HRFcurve.infoBV.coor_r = 1;
    HRFcurve.infoBV.coor_theta  = 90;
    HRFcurve.infoBV.coor_phi = 45;
    HRFcurve.infoBV.DataPoints = numel(Scn);
    [pathstr, name, ext] = fileparts(fileevent);
    fileoutHRFcurve = fullfile(pathout, [label,'.dat']);
    fwrite_EEG(fileoutHRFcurve,HRFcurve,1,HRFcurve.infoBV.DataPoints );
elseif strcmp(listoption{idval},'Global regressor NIRS.mat average exclude rejected channel')
    load(filedata); 
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    HbOch = 1:NC/2;
    HbRch = (NC/2+1):NC;
    fs = NIRS.Cf.dev.fs;
    regressch = 'Global';
    for f=1:size(rDtp,1)
        idok = find(NIRS.Cf.H.C.ok(HbOch,f));
        listHBO = HbOch(idok);
        idok = find(NIRS.Cf.H.C.ok(HbRch,f));
        listHBR = HbRch(idok);
        [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]= fopen_EEG(rDtp{f});
        idtrigger = find(strcmp('trigger',EEG.marker(:,1)));
        EEG.marker = EEG.marker(idtrigger,:);
        EEG.ind_dur_ch = EEG.ind_dur_ch(idtrigger,:);
        
        x1HbO=nanmean(EEG.data(:,listHBO),2);
        x1HbR=nanmean(EEG.data(:,listHBR),2);
        AUX.data = [x1HbO,x1HbR];
        AUX.ind_dur_ch=EEG.ind_dur_ch;
        AUX.marker=EEG.marker;
        AUX.infoBV = EEG.infoBV;
        AUX.infoBV.DataType = 'TIMEDOMAIN';
        AUX.infoBV.NumberOfChannels = '2';
        AUX.infoBV.name_ele = {'GlobalHBO', 'GlobalHbR'};
        AUX.infoBV.coor_r = [1 1];
        AUX.infoBV.coor_theta  = [90 90];
        AUX.infoBV.coor_phi = [45 45];
        AUX.infoBV.DataPoints = size(AUX.data,1);  
        fileoutHRFcurve = fullfile(pathout, [label,'.dat']);
        fwrite_EEG(fileoutHRFcurve,AUX,1,AUX.infoBV.DataPoints );
     end 
elseif strcmp(listoption{idval},'Short channel regressor NIRS.mat (channel is measurement in the data)')
    if ~isdir(fileevent)
        mkdir(fileevent)
    end
    [Detstr, Srsstr]=strtok(label);              
    iddet = StrBoxy2SDDet_ISS( deblank(strtrim(Detstr)));
    idsrs = StrBoxy2SDPairs( deblank(strtrim(Srsstr)));
    load(filedata); 
    idReg = find(NIRS.Cf.H.C.id(2,:)== idsrs & NIRS.Cf.H.C.id(3,:)==iddet);
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    HbOch = 1:NC/2;
    HbRch = (NC/2+1):NC;
    fs = NIRS.Cf.dev.fs;
    regressch = 'ShortDistance';
  
    for f=1:size(rDtp,1)
        idok = find(NIRS.Cf.H.C.ok(HbOch,f));
        listHBO = HbOch(idok);
        idok = find(NIRS.Cf.H.C.ok(HbRch,f));
        listHBR = HbRch(idok);
        [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]= fopen_EEG(rDtp{f});
        idtrigger = find(strcmp('trigger',EEG.marker(:,1)));
        EEG.marker = EEG.marker(idtrigger,:);
        EEG.ind_dur_ch = EEG.ind_dur_ch(idtrigger,:);

        idHbO = idReg(1);
        idHbR = idReg(2);
        x1HbO=EEG.data(:,idHbO);
        x1HbR=EEG.data(:,idHbR);
        AUX.data = [x1HbO,x1HbR];%VECTORIZED=ch1,pt1, ch1,pt2
        AUX.ind_dur_ch=EEG.ind_dur_ch;
        AUX.marker=EEG.marker;
        AUX.infoBV = EEG.infoBV;
        AUX.infoBV.DataType = 'TIMEDOMAIN';
        AUX.infoBV.NumberOfChannels = '2';
        AUX.infoBV.name_ele = {[deblank(label),'HBO'],[deblank(label), 'HbR']};
        AUX.infoBV.coor_r = [1 1];
        AUX.infoBV.coor_theta  = [90 90];
        AUX.infoBV.coor_phi = [45 45];
        AUX.infoBV.DataPoints = size(AUX.data,1);  
        fileoutHRFcurve = fullfile(pathout, [label,'.dat']);
        fwrite_EEG(fileoutHRFcurve,AUX,1,AUX.infoBV.DataPoints );
    end    
elseif strcmp(listoption{idval},'Audio .wav (sum mel coefficients)')
    [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]= fopen_EEG(filedata);
    idtrigger = find((strcmp('trigger',EEG.marker(:,1))& EEG.ind_dur_ch(:,1)>0)|(strcmp('Stimulus',EEG.marker(:,1))& EEG.ind_dur_ch(:,1)>0));
    %idtrigger = find(strcmp('Stimulus',EEG.marker(:,1))& EEG.ind_dur_ch(:,1)>0);
    EEG.marker = EEG.marker(idtrigger,:);
    EEG.ind_dur_ch = EEG.ind_dur_ch(idtrigger,:);
    [speech,Fsaudio] = audioread(fileevent);
    [filepath,name,ext] = fileparts(fileevent);
    fileoutMEL = fullfile(pathout,['sumMEL',name,'.dat']);
   
    Tw=str2num(get(handles.edit_parameter_1,'string'));
    Ts=str2num(get(handles.edit_parameter_2,'string'));
    alpha=str2num(get(handles.edit_parameter_3,'string'));
    R=str2num(get(handles.edit_parameter_4,'string'));
    M=str2num(get(handles.edit_parameter_5,'string'));
    C=str2num(get(handles.edit_parameter_6,'string'));
    L=str2num(get(handles.edit_parameter_7,'string'));
    N=str2num(get(handles.edit_parameter_8,'string'));
    hamming = @(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));
   
    [ MFCCs, FBEs, frames ] = mfcc( speech(:,1), Fsaudio, Tw, Ts, alpha, hamming, R, M, C, L );
    sumMel = sum(log10(FBEs));
     
     %exclude inf or nan number
    idbad = find(isinf(sumMel)|isnan(sumMel));
    sumMel(idbad) = 0;  
    msumMel = mean(sumMel);
    sumMel(idbad) = msumMel;
    sumMel = sumMel - msumMel;
    time = 1/Fsaudio:1/Fsaudio:((1/Fsaudio)*numel(speech));
    Fsmel = numel(sumMel)/time(end);
    offsetval_sec = str2num(get(handles.edit_name,'string')); %CHANGE HERE TO ADJUST EEG TRIG AND SOUND DISPLAY
	if offsetval_sec < 0
       % the audio segment is shorter, pad the beggining with zeros to
       % avoid missing trig
       % AddOffsetval in the biggining 
       sampleadd = round(Fsmel*abs(offsetval_sec));
       sumMel
       sumMel = [zeros(1,sampleadd),sumMel]; 
       offsetval_sec = 1/Fsmel;
    end
    tmel = 1/Fsmel:1/Fsmel:(1/Fsmel*numel(sumMel));
    
    EEG.ind_dur_ch
    Fsnirs = 1/(EEG.infoBV.SamplingInterval/1000000);
    tNIRS = size(EEG.data,1)*1/Fsnirs;
    if tmel(end)<tNIRS %tmel too short add sample to the end
       timetoadd = tNIRS-tmel(end);
       sampleadd = round(Fsmel*timetoadd);
       sumMel = [sumMel,zeros(1,sampleadd),] ;
       tmel = 1/Fsmel:1/Fsmel:(1/Fsmel*numel(sumMel));
    end
         

      
     
     %AFFICHAGE TRIG 2 TO SEE OFFSET 
    offsetval = sum(tmel<=offsetval_sec);
    if offsetval ==0
        offsetval = 1;
    end
    figure;hold on
    plot(tmel,sumMel)
    plot([tmel(offsetval),tmel(offsetval)],[-max(abs(sumMel)),max(abs(sumMel))])
    %AJUST TIME FOR ALL TRIG
    Fsdata = 1/(EEG.infoBV.SamplingInterval/1000000);
    trig = EEG.ind_dur_ch(:,1)*1/Fsdata; %in seconds
    trigMel = zeros(numel(sumMel),1);
    outlist = [];     
    
    
    
    
    for i=1:numel(trig)
        id = sum(tmel<trig(i))+offsetval;
        if id>0
        trigMel(id)=1;
        else
            outlist = [ outlist, i];
        end
    end
    %figure;plot(trigMel)
    %AFFICHAGE AVEC OFFSET 
    SUMMEL= EEG;
     if ~isempty(outlist)
        SUMMEL.marker(outlist,:) = [];
        SUMMEL.ind_dur_ch(outlist,:) = [];
     end
    SUMMEL.ind_dur_ch(:,1)=find(trigMel);
     %EEG POUR LE MEL SUM 
    SUMMEL.data = sumMel';
    SUMMEL.infoBV = EEG.infoBV;
    SUMMEL.infoBV.DataType = 'TIMEDOMAIN';
    SUMMEL.infoBV.NumberOfChannels = '1';
    SUMMEL.infoBV.SamplingInterval =(1/Fsmel)*1000000;
    SUMMEL.infoBV.name_ele = {'MelSum'};
    SUMMEL.infoBV.coor_r = 1;
    SUMMEL.infoBV.coor_theta  = 90;
    SUMMEL.infoBV.coor_phi = 45;
    SUMMEL.infoBV.DataPoints = numel(sumMel);
    %%%%%%%%OUTPUT AUX MEL SUM 
    fwrite_EEG(fileoutMEL,SUMMEL,1,SUMMEL.infoBV.DataPoints );
    %WRITE AUX avg
elseif strcmp(listoption{idval},'Concatenate AUX/EEG')
    pathin = get(handles.edit_rawdata,'string');
    offsetsizebloc = 0;
    AUX.data = [];
    AUX.marker = [];
    AUX.ind_dur_ch = [];
    rDtp = get(handles.listbox_AUXCONCATENATE,'string');
    for f=1:size(rDtp,1)
         fileAUX = fullfile(pathin, rDtp{f});
        [data,infoBV,marker,ind_dur_ch]= fopen_EEG(fileAUX); 
         AUX.data=[ AUX.data;data];
         ind_dur_ch(:,1)=ind_dur_ch(:,1)+ offsetsizebloc;      
         offsetsizebloc =  offsetsizebloc + size(data,1);
         AUX.marker =  [ AUX.marker;marker];
         AUX.ind_dur_ch = [AUX.ind_dur_ch ;ind_dur_ch];               
    end
    outfileAUX =fullfile(pathout,['all', rDtp{f}]);
    infoBV.DataPoints =  size(AUX.data,1);
    AUX.infoBV = infoBV;
    fwrite_EEG(outfileAUX, AUX,1,AUX.infoBV.DataPoints );
    set(handles.text_parameterdescription,'visible','on')
    set(handles.text_parameterdescription,'string',['We finish to create ',outfileAUX])
elseif strcmp(listoption{idval},'VIDEO CONVERT and MERGE file CODEC 32 to 64 bit')
    pathin = get(handles.edit_rawdata,'string');
    rDtp = get(handles.listbox_AUXCONCATENATE,'string');
   
%     if iscell(rDtp)
%         outfilevideo=fullfile(pathin,['all', rDtp{1}]);
%     else
%          outfilevideo=fullfile(pathin,['all', rDtp]);
%     end
    set(handles.text_parameterdescription,'visible','on')
    set(handles.text_parameterdescription,'string','Please wait video conversion take a while')
    
   if ~iscell(rDtp)
       rDtp = {rDtp};
   end
   tic
     filename = rDtp{1};
     outfilevideo=fullfile(pathin,['all',filename(1:end-4),'v64',filename(end-3:end)]);
     tmp = VideoWriter(outfilevideo);  
    for f=1:numel(rDtp)  
        v = VideoReader(fullfile(pathin,rDtp{f}));
        if f==1  
            tmp.FrameRate = v.FrameRate;       
            open(tmp)
        end
        id = 1;
        while hasFrame(v)          
            frame = readFrame(v);
            % [frame,Cb,Cr,audio] = videofreader(v);
            writeVideo(tmp,frame);   
            id = id+1;
            if ~mod(id,300)
                id
            end
        end 
      
    end
     toc
    yall =[];   
        
    for f=1:numel(rDtp) 
         
            filename = rDtp{1};
            outfileaudio=fullfile(pathout,['all',filename(1:end-4),'v64','.wav']);
            filename = fullfile(pathin,rDtp{f});
            [y,Fs] = audioread(filename);
            yall = [yall;y];
            audiowrite(outfileaudio,yall,Fs);
    end
    
    close(tmp)
    set(handles.text_parameterdescription,'visible','on')
    set(handles.text_parameterdescription,'string',['We finish to convert file ',outfilevideo,' is created'])

end
%msgbox([fileoutHRFcurve, ' created'])
% HRF.marker = 
% HRF.ind_dur_ch = 
% fwrite_EEG(filename,HRFcurve

function edit_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_name as text
%        str2double(get(hObject,'String')) returns contents of edit_name as a double


% --- Executes during object creation, after setting all properties.
function edit_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuoption.
function popupmenuoption_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuoption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuoption
hideall_handles(handles)
listoption = get(handles.popupmenuoption,'string');
idval = get(handles.popupmenuoption,'value');
if strcmp(listoption{idval},'Create HRF using onset duration xls file')%HRF from event xls
    set(handles.text_rawdata,'visible','on')
    set(handles.text_rawdata,'string','1. Enter whole bloc nir data from the same recording to know trigger info for segmenting') 
    set(handles.edit_rawdata,'visible','on')
    set(handles.btn_Browserawdata,'visible','on')  
    
    set(handles.text_eventfile,'visible','on')
    set(handles.text_eventfile,'string','2. Enter onset xls file, onset (s) and duration (s), weigth Output will be writen here')
    set(handles.edit_eventfile,'visible','on')
    set(handles.btn_Browseeventfile,'visible','on')
 
    set(handles.text_name,'visible','on')    
    set(handles.text_name,'string','3. Name Regressor')
    set(handles.edit_name,'visible','on')
    set(handles.edit_name,'string','HRF_evt1')
    
    set(handles.text_parameterdescription,'visible','on')
    set(handles.text_parameterdescription,'string','4. HRF parameter Glover, NeuroImage, 9:416-429 gamma1/max(gamma1)-DIP*gamma2/max(gamma2) scaled so that its total integral is 1. ')
    set(handles.text_parameter_1,'visible','on')
    set(handles.text_parameter_1,'string','Time to the peak of the first gamma density')
    set(handles.edit_parameter_1,'visible','on')
    set(handles.edit_parameter_1,'string','5.4')
     
    set(handles.text_parameter_2,'string','FWHM1: approximate FWHM of the first gamma density')
    set(handles.text_parameter_2,'visible','on')
    set(handles.edit_parameter_2,'visible','on')
    set(handles.edit_parameter_2,'string','5.2')
  
    set(handles.text_parameter_3,'string','time to the peak of the second gamma density')
    set(handles.text_parameter_3,'visible','on')
    set(handles.edit_parameter_3,'visible','on')
    set(handles.edit_parameter_3,'string','10.8')
      
    set(handles.text_parameter_4,'string','FWHM2: approximate FWHM of the second gamma density');
    set(handles.text_parameter_4,'visible','on')
    set(handles.edit_parameter_4,'visible','on')
    set(handles.edit_parameter_4,'string','7.35')
     
    set(handles.text_parameter_5,'string','DIP: coefficient of the second gamma density');
    set(handles.text_parameter_5,'visible','on')
    set(handles.edit_parameter_5,'visible','on')
    set(handles.edit_parameter_5,'string','0.35')
    
elseif  strcmp(listoption{idval},'Create HRF convolution .dat AUX signal ')
    set(handles.text_eventfile,'visible','on')
    set(handles.text_eventfile,'string','1. Enter auxiliary to convolve ')
    set(handles.edit_eventfile,'visible','on')
    set(handles.btn_Browseeventfile,'visible','on')
    
    set(handles.text_name,'visible','on')    
    set(handles.text_name,'string','2. Name Regressor')
    set(handles.edit_name,'visible','on')
    set(handles.edit_name,'string','HRF_evt1')
    
    set(handles.text_parameterdescription,'visible','on')
    set(handles.text_parameterdescription,'string','HRF parameter Glover, NeuroImage, 9:416-429 gamma1/max(gamma1)-DIP*gamma2/max(gamma2) scaled so that its total integral is 1. ')
    set(handles.text_parameter_1,'visible','on')
    set(handles.text_parameter_1,'string','Time to the peak of the first gamma density')
    set(handles.edit_parameter_1,'visible','on')
    set(handles.edit_parameter_1,'string','5.4')
     
    set(handles.text_parameter_2,'string','FWHM1: approximate FWHM of the first gamma density')
    set(handles.text_parameter_2,'visible','on')
    set(handles.edit_parameter_2,'visible','on')
    set(handles.edit_parameter_2,'string','5.2')
  
    set(handles.text_parameter_3,'string','time to the peak of the second gamma density')
    set(handles.text_parameter_3,'visible','on')
    set(handles.edit_parameter_3,'visible','on')
    set(handles.edit_parameter_3,'string','10.8')
      
    set(handles.text_parameter_4,'string','FWHM2: approximate FWHM of the second gamma density');
    set(handles.text_parameter_4,'visible','on')
    set(handles.edit_parameter_4,'visible','on')
    set(handles.edit_parameter_4,'string','7.35')
     
    set(handles.text_parameter_5,'string','DIP: coefficient of the second gamma density');
    set(handles.text_parameter_5,'visible','on')
    set(handles.edit_parameter_5,'visible','on')
    set(handles.edit_parameter_5,'string','0.35')
    
    set(handles.text_parameter_6,'string','Enter 1 to binarize the input, enter 0 to keep analog');
    set(handles.text_parameter_6,'visible','on')
    set(handles.edit_parameter_6,'visible','on')
    set(handles.edit_parameter_6,'string','0')
elseif  strcmp(listoption{idval},'Global regressor NIRS.mat average exclude rejected channel') %HRF from event xls
    set(handles.text_rawdata,'visible','on')
    set(handles.text_rawdata,'string','Enter NIRS.mat at the step before normalisation') 
    set(handles.edit_rawdata,'visible','on')
    set(handles.btn_Browserawdata,'visible','on')  
   
    set(handles.text_eventfile,'visible','off')
    set(handles.text_eventfile,'string','Output directory')
    set(handles.edit_eventfile,'visible','off')
    set(handles.btn_Browseeventfile,'visible','off')
    
    set(handles.text_name,'visible','on')    
    set(handles.text_name,'string','Name')
    set(handles.edit_name,'visible','on')
    set(handles.edit_name,'string','Global')
elseif strcmp(listoption{idval},'Short channel regressor NIRS.mat (channel is measurement in the data)')
    set(handles.text_rawdata,'visible','on')
    set(handles.text_rawdata,'string','Enter NIRS.mat at the step before normalisation') 
    set(handles.edit_rawdata,'visible','on')
    set(handles.btn_Browserawdata,'visible','on')  
    
    set(handles.text_eventfile,'visible','off')
    set(handles.text_eventfile,'string','Output directory')
    set(handles.edit_eventfile,'visible','off')
    set(handles.btn_Browseeventfile,'visible','off')
    
    set(handles.text_name,'visible','on')    
    set(handles.text_name,'string','Name ISS')
    set(handles.edit_name,'visible','on')
    set(handles.edit_name,'string','A a1b2')
elseif strcmp(listoption{idval},'Audio .wav (sum mel coefficients)')
    set(handles.text_rawdata,'visible','on')
    set(handles.text_rawdata,'string','Enter .nir or .dat of simulatenous NIRS or EEG recording')
    set(handles.edit_rawdata,'visible','on')
    set(handles.btn_Browserawdata,'visible','on') 
    
    set(handles.text_eventfile,'visible','on')
    set(handles.text_eventfile,'string','Sound data .wav')
    set(handles.edit_eventfile,'visible','on')
    set(handles.btn_Browseeventfile,'visible','on')
    
    set(handles.text_name,'visible','on')    
    set(handles.text_name,'string','Onset delay (s) (WAV event start - nir event start)')
    set(handles.edit_name,'visible','on')
    set(handles.edit_name,'string','0')
    
    set(handles.text_parameter_1,'visible','on')
    set(handles.text_parameter_1,'string','Analysis frame duration (ms)')
    set(handles.edit_parameter_1,'visible','on')
    set(handles.edit_parameter_1,'string','25')
     
    set(handles.text_parameter_2,'string','Analysis frame shift (ms)')
    set(handles.text_parameter_2,'visible','on')
    set(handles.edit_parameter_2,'visible','on')
    set(handles.edit_parameter_2,'string','10')
  
    set(handles.text_parameter_3,'string','Preemphasis coefficient alpha')
    set(handles.text_parameter_3,'visible','on')
    set(handles.edit_parameter_3,'visible','on')
    set(handles.edit_parameter_3,'string','0.97')
      
    set(handles.text_parameter_4,'string','Frequency range to consider');
    set(handles.text_parameter_4,'visible','on')
    set(handles.edit_parameter_4,'visible','on')
    set(handles.edit_parameter_4,'string','300, 3700')
     
    set(handles.text_parameter_5,'string','Number of filterbank channels');
    set(handles.text_parameter_5,'visible','on')
    set(handles.edit_parameter_5,'visible','on')
    set(handles.edit_parameter_5,'string','20')
     
    set(handles.text_parameter_6,'string','Number of cepstral coefficients');
    set(handles.text_parameter_6,'visible','on')
    set(handles.edit_parameter_6,'visible','on')
    set(handles.edit_parameter_6,'string','13')
     
    set(handles.text_parameter_7,'string','Cepstral sine lifter parameter');
    set(handles.text_parameter_7,'visible','on')
    set(handles.edit_parameter_7,'visible','on')
    set(handles.edit_parameter_7,'string','22')
     
    set(handles.text_parameter_8,'string','Hamming window');
    set(handles.text_parameter_8,'visible','on')
    set(handles.edit_parameter_8,'visible','on')
    set(handles.edit_parameter_8,'string','20')
elseif strcmp(listoption{idval},'Concatenate AUX/EEG')
    set(handles.listbox_AUXCONCATENATE,'visible','on')
    set(handles.btn_Browserawdata,'visible','on')
    set(handles.text_rawdata,'visible','on')
    set(handles.text_rawdata,'string',['Open generic data export data to merge, keep the order of the recording' ])
    set(handles.edit_rawdata,'visible','on')
elseif strcmp(listoption{idval},'VIDEO CONVERT and MERGE file CODEC 32 to 64 bit')
    set(handles.listbox_AUXCONCATENATE,'visible','on')
    set(handles.btn_Browserawdata,'visible','on')
    set(handles.text_rawdata,'visible','on')
    set(handles.text_rawdata,'string',['Open .avi LEADMCMPCodec 32bit to be convert in 64bit, use MATLAB 32bit'])
    set(handles.edit_rawdata,'visible','on')
end





function hideall_handles(handles)
     set(handles.text_rawdata,'visible','off')
     set(handles.edit_rawdata,'visible','off')
     set(handles.btn_Browserawdata,'visible','off')     
     set(handles.text_eventfile,'visible','off')
     set(handles.edit_eventfile,'visible','off')
     set(handles.btn_Browseeventfile,'visible','off')     
     set(handles.text_name,'visible','off')
     set(handles.edit_name,'visible','off')
      set(handles.text_parameterdescription,'visible','off')
     set(handles.text_parameter_1,'visible','off')
     set(handles.edit_parameter_1,'visible','off')
     set(handles.text_parameter_2,'visible','off')
     set(handles.edit_parameter_2,'visible','off')
     set(handles.text_parameter_3,'visible','off')
     set(handles.edit_parameter_3,'visible','off')
     set(handles.text_parameter_4,'visible','off')
     set(handles.edit_parameter_4,'visible','off')
     set(handles.text_parameter_5,'visible','off')
     set(handles.edit_parameter_5,'visible','off')
     set(handles.text_parameter_6,'visible','off')
     set(handles.edit_parameter_6,'visible','off')
     set(handles.text_parameter_7,'visible','off')
     set(handles.edit_parameter_7,'visible','off')
     set(handles.text_parameter_8,'visible','off')
     set(handles.edit_parameter_8,'visible','off')
     set(handles.text_parameterdescription,'visible','off')
     set(handles.listbox_AUXCONCATENATE,'visible','off')

% --- Executes during object creation, after setting all properties.
function popupmenuoption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_parameter_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parameter_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parameter_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_parameter_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_parameter_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_parameter_5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parameter_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parameter_5 as text
%        str2double(get(hObject,'String')) returns contents of edit_parameter_5 as a double


% --- Executes during object creation, after setting all properties.
function edit_parameter_5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_parameter_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_parameter_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_parameter_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_rawdata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rawdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_parameter_6_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parameter_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parameter_6 as text
%        str2double(get(hObject,'String')) returns contents of edit_parameter_6 as a double


% --- Executes during object creation, after setting all properties.
function edit_parameter_6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_parameter_7_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parameter_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parameter_7 as text
%        str2double(get(hObject,'String')) returns contents of edit_parameter_7 as a double


% --- Executes during object creation, after setting all properties.
function edit_parameter_7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_parameter_8_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parameter_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parameter_8 as text
%        str2double(get(hObject,'String')) returns contents of edit_parameter_8 as a double


% --- Executes during object creation, after setting all properties.
function edit_parameter_8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function X_cache=fmridesign(frametimes,slicetimes,events,S, ...
   hrf_parameters,shift)

%FMRIDESIGN makes a set of design matrices for fmrilm.
%
% X_CACHE = FMRIDESIGN( FRAME_TIMES [, SLICE_TIMES [, EVENTS , [S , 
%                       [, HRF_PARAMETERS [, SHIFT]]]]]] )
% 
% FRAME_TIMES is a row vector of frame acquisition times in seconds. 
% With just the frametimes, it gives the hemodynamic response function.
% 
% SLICE_TIMES is a row vector of relative slice acquisition times,
% i.e. absolute acquisition time of a slice is FRAME_TIMES + SLICE_TIMES.
% Default is 0.
% 
% EVENTS is a matrix whose rows are events and whose columns are:
% 1. id - an integer from 1:(number of events) to identify event type;
% 2. times - start of event, synchronised with frame and slice times;
% 3. durations (optional - default is 0) - duration of event;
% 4. heights (optional - default is 1) - height of response for event.
% For each event type, the response is a box function starting at the event 
% times, with the specified durations and heights, convolved with the 
% hemodynamic response function (see below). If the duration is zero, the 
% response is the hemodynamic response function whose integral is 
% the specified height - useful for `instantaneous' stimuli such as visual 
% stimuli. The response is then subsampled at the appropriate frame and slice
% times to create a design matrix for each slice, whose columns correspond
% to the event id number. EVENT_TIMES=[] will ignore event times and just 
% use the stimulus design matrix S (see next). Default is [1 0].
% 
% S: Events can also be supplied by a stimulus design matrix, 
% whose rows are the frames, and column are the event types. Events 
% are created for each column, beginning at the frame time for each row
% of S, with a duration equal to the time to the next frame, and a height
% equal to the value of S for that row and column. Note that a
% constant term is not usually required, since it is removed by the
% polynomial trend terms provided N_POLY>=0. Default is [].
% 
% HRF_PARAMETERS is a matrix whose rows are 5 parameters for the 
% hemodynamic response function, one row for each event type and column
% of S (if there is just one row, this is repeated as necessary). 
% The hrf is modeled as the difference of two gamma density functions 
% (Glover, NeuroImage, 9:416-429). The components of HRF_PARAMETERS are:
% 1. PEAK1: time to the peak of the first gamma density;
% 2. FWHM1: approximate FWHM of the first gamma density;
% 3. PEAK2: time to the peak of the second gamma density;
% 4. FWHM2: approximate FWHM of the second gamma density;
% 5. DIP: coefficient of the second gamma density;
%    Final hrf is:   gamma1/max(gamma1)-DIP*gamma2/max(gamma2)
%    scaled so that its total integral is 1. 
% If PEAK1=0 then there is no smoothing of that event type with the hrf.
% If PEAK1>0 but FWHM1=0 then the design is simply lagged by PEAK1.
% Default is: [5.4 5.2 10.8 7.35 0.35] chosen by Glover (1999) for 
% an auditory stimulus. 
% If HRF_PARAMETERS is a structure, then HRF_PARAMETERS.T is a matrix
% whose rows are the times in seconds of a user-supplied HRF, one 
% row for each event type and column of S (if there is just one row, 
% this is repeated as necessary).  Times must start at 0 but need not 
% be equally spaced; spacing of 0.02s is recommended. HRF_PARAMETERS.H 
% are the corresponding values of the HRF at those times. 
%
% SHIFT is a matrix whose rows are the min and max shift for the hrf in 
% seconds, one row for each event type and column of S (if there is just 
% one row, this is repeated as necessary). Default is [-4.5 4.5]*FWHM1/5.2.
% 
% X_CACHE.TR: TR, average time between frames (secs). 
% 
% X_CACHE.X: A cache of the design matrices stored as a 4D array. 
% Dim 1: frames; Dim 2: response variables; Dim 3: 4 values, corresponding 
% to the stimuli convolved with: hrf, derivative of hrf, first and second 
% spectral basis functions over the range in SHIFT; Dim 4: slices. 
%
% X_CACHE.W: A 3D array of coefficients of the basis functions in X_CACHE.X.
% Dim 1: frames; Dim 2: response variables; Dim 3: 5 values: 
% coefficients of the hrf and its derivative, coefficients of the first and 
% second spectral basis functions, shift values.

%############################################################################
% COPYRIGHT:   Copyright 2002 K.J. Worsley 
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              worsley@math.mcgill.ca, liao@math.mcgill.ca
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that the above copyright
%              notice appear in all copies. The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

% Defaults:

if nargin < 2
   slicetimes=0;
end
if nargin < 3
   events=[1 0];
end
if nargin < 4
   S=[];
end
if nargin < 5
   hrf_parameters=[5.4 5.2 10.8 7.35 0.35];
end
if nargin < 6
   if isstruct(hrf_parameters)
      shift=[-4.5 4.5];
   else
     shift=[-4.5 4.5]*max(max(hrf_parameters(:,2))/5.2,1);
  end
end

n=length(frametimes);
numslices=length(slicetimes);

% Keep time points that are not excluded:

if ~isempty(events)
   numevents=size(events,1);
   eventid=events(:,1);
   numeventypes=max(eventid);
   eventime=events(:,2);
   if size(events,2)>=3
      duration=events(:,3);
   else
      duration=zeros(numevents,1);
   end
   if size(events,2)>=4
      height=events(:,4);
   else
      height=ones(numevents,1);
   end
   mineventime=min(eventime);
   maxeventime=max(eventime+duration);
else
   numeventypes=0;
   mineventime=Inf;
   maxeventime=-Inf;
end

if ~isempty(S)
   numcolS=size(S,2);
else
   numcolS=0;
end

% Set up response matrix:

dt=0.02;
startime=min(mineventime,min(frametimes)+min([slicetimes 0]));
finishtime=max(maxeventime,max(frametimes)+max([slicetimes 0]));
numtimes=ceil((finishtime-startime)/dt)+1;
numresponses=numeventypes+numcolS;
response=zeros(numtimes,numresponses);

if ~isempty(events)
   height=height./(1+(duration==0)*(dt-1));
   for k=1:numevents
      type=eventid(k);
      n1=ceil((eventime(k)-startime)/dt)+1;
      n2=ceil((eventime(k)+duration(k)-startime)/dt)+(duration(k)==0);
      if n2>=n1
         response(n1:n2,type)=response(n1:n2,type)+height(k)*ones(n2-n1+1,1);
      end
   end
end

if ~isempty(S)
   for j=1:numcolS
      for i=find(S(:,j)')
         n1=ceil((frametimes(i)-startime)/dt)+1;
         if i<n
            n2=ceil((frametimes(i+1)-startime)/dt);
         else
            n2=numtimes;
         end
         if n2>=n1 
            response(n1:n2,numeventypes+j)= ...
               response(n1:n2,numeventypes+j)+S(i,j)*ones(n2-n1+1,1);
         end
      end
   end
end

if isstruct(hrf_parameters)
   if size(hrf_parameters.T,1)==1
      hrf_parameters.T=repmat(hrf_parameters.T,numresponses,1);
   end
   if size(hrf_parameters.H,1)==1
      hrf_parameters.H=repmat(hrf_parameters.H,numresponses,1);
   end
else
   if size(hrf_parameters,1)==1
      hrf_parameters=repmat(hrf_parameters,numresponses,1);
   end
end
if size(shift,1)==1
   shift=repmat(shift,numresponses,1);
end

eventmatrix=zeros(numtimes,numresponses,4);
nd=41;
X_cache.W=zeros(nd,numresponses,5);

for k=1:numresponses
   Delta1=shift(k,1);
   Delta2=shift(k,2);
   if isstruct(hrf_parameters)
      numlags=ceil((max(hrf_parameters.T(k,:))+Delta2-Delta1)/dt)+1;
   else
      peak1=hrf_parameters(k,1);
      fwhm1=hrf_parameters(k,2);
      peak2=hrf_parameters(k,3);
      fwhm2=hrf_parameters(k,4);
      dip=hrf_parameters(k,5);
      numlags=ceil((max(peak1+3*fwhm1,peak2+3*fwhm2)+Delta2-Delta1)/dt)+1;
   end
   numlags=min(numlags,numtimes);
   time=(0:(numlags-1))'*dt;

   % Taylor:
   if isstruct(hrf_parameters)
      hrf=interp1(hrf_parameters.T(k,:),hrf_parameters.H(k,:),time,'spline',0);
      d_hrf=-gradient(hrf,dt);
   else
      tinv=(time>0)./(time+(time<=0));
      if peak1>0 & fwhm1>0
         alpha1=peak1^2/fwhm1^2*8*log(2);
         beta1=fwhm1^2/peak1/8/log(2);
         gamma1=(time/peak1).^alpha1.*exp(-(time-peak1)./beta1);
         d_gamma1=-(alpha1*tinv-1/beta1).*gamma1;
      else 
         gamma1=min(abs(time-peak1))==abs(time-peak1);
         d_gamma1=zeros(numlags,1);
      end
      if peak2>0 & fwhm2>0
         alpha2=peak2^2/fwhm2^2*8*log(2);
         beta2=fwhm2^2/peak2/8/log(2);
         gamma2=(time/peak2).^alpha2.*exp(-(time-peak2)./beta2);
         d_gamma2=-(alpha2*tinv-1/beta2).*gamma2;
      else 
         gamma2=min(abs(time-peak2))==abs(time-peak2);
         d_gamma2=zeros(numlags,1);
      end
      hrf=gamma1-dip*gamma2;
      d_hrf=d_gamma1-dip*d_gamma2;
   end
   HS=[hrf d_hrf]/sum(hrf);
   temp=conv2(response(:,k),HS);
   eventmatrix(:,k,1:2)=temp(1:numtimes,:);
   
   % Shifted hrfs:
   H=zeros(numlags,nd);
   delta=((1:nd)-1)/(nd-1)*(Delta2-Delta1)+Delta1;
   for id=1:nd
      if isstruct(hrf_parameters)
         t=time+Delta1-delta(id);
         hrf=interp1(hrf_parameters.T(k,:),hrf_parameters.H(k,:),t,'spline',0);
      else
         t=(time+Delta1-delta(id)).*((time+Delta1)>delta(id));
         if peak1>0 & fwhm1>0
            gamma1=(t/peak1).^alpha1.*exp(-(t-peak1)./beta1);
         else 
            gamma1=min(abs(t-peak1))==abs(t-peak1);
         end
         if peak2>0 & fwhm2>0
            gamma2=(t/peak2).^alpha2.*exp(-(t-peak2)./beta2);
         else 
            gamma2=min(abs(t-peak2))==abs(t-peak2);
         end
         hrf=gamma1-dip*gamma2;
      end
      H(:,id)=hrf/sum(hrf);
   end
   
   % Taylor coefs:
   origin=-round(Delta1/dt);
   HS0=[zeros(origin,2); HS(1:(numlags-origin),:)];
   WS=pinv(HS0)*H;
   X_cache.W(:,k,1:2)=WS';
   prcnt_var_taylor=sum(sum(H.*(HS0*WS)))/sum(sum(H.*H))*100;

   % svd:
   [U,SS,V]=svd(H,0);
   prcnt_var_spectral=(SS(1,1)^2+SS(2,2)^2)/sum(diag(SS).^2)*100;
   sumU=sum(U(:,1));
   US=U(:,1:2)/sumU;
   WS=V(:,1:2)*SS(1:2,1:2)*sumU;
   if delta*WS(:,2)<0
      US(:,2)=-US(:,2);
      WS(:,2)=-WS(:,2);
   end
   temp=conv2(response(:,k),US);
   eventmatrix(:,k,3:4)=temp((1:numtimes)-round(Delta1/dt),:);
   X_cache.W(:,k,3:4)=WS;
   X_cache.W(:,k,5)=delta';
   
   if ~all(WS(:,1)>0)
      fprintf(['Warning: use only for magnitudes, not delays \n first coef not positive for stimulus ' num2str(k)]);
   end
   cubic_coef=pinv([delta' delta'.^3])*(WS(:,2)./WS(:,1));
   if prod(cubic_coef)<0
      fprintf(['\nWarning: use only for magnitudes, not delays \n svd ratio not invertible for stimulus ' num2str(k)]);
   end
end 

X_cache.X=zeros(n,numresponses,4,numslices);

for slice = 1:numslices
   subtime=ceil((frametimes+slicetimes(slice)-startime)/dt)+1;
   X_cache.X(:,:,:,slice)=eventmatrix(subtime,:,:);
end

X_cache.TR=(max(frametimes)-min(frametimes))/(length(frametimes)-1);

return


% --- Executes on selection change in listbox_AUXCONCATENATE.
function listbox_AUXCONCATENATE_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_AUXCONCATENATE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_AUXCONCATENATE contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_AUXCONCATENATE


% --- Executes during object creation, after setting all properties.
function listbox_AUXCONCATENATE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_AUXCONCATENATE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pathout_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathout as text
%        str2double(get(hObject,'String')) returns contents of edit_pathout as a double


% --- Executes during object creation, after setting all properties.
function edit_pathout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_getdirout.
function btn_getdirout_Callback(hObject, eventdata, handles)
% hObject    handle to btn_getdirout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

directory_name = uigetdir();
set(handles.edit_pathout,'string',directory_name) ;
