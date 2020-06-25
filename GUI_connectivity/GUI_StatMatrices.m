function varargout = GUI_StatMatrices(varargin)
% GUI_STATMATRICES M-file for GUI_StatMatrices.fig
%      GUI_STATMATRICES, by itself, creates a new GUI_STATMATRICES or raises the existing
%      singleton*.
%
%      H = GUI_STATMATRICES returns the handle to a new GUI_STATMATRICES or the handle to
%      the existing singleton*.
%
%      GUI_STATMATRICES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_STATMATRICES.M with the given input arguments.
%
%      GUI_STATMATRICES('Property','Value',...) creates a new GUI_STATMATRICES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_StatMatrices_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_StatMatrices_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_StatMatrices

% Last Modified by GUIDE v2.5 25-Sep-2018 12:42:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_StatMatrices_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_StatMatrices_OutputFcn, ...
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


% --- Executes just before GUI_StatMatrices is made visible.
function GUI_StatMatrices_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_StatMatrices (see VARARGIN)

% Choose default command line output for GUI_StatMatrices
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_StatMatrices wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_StatMatrices_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_browseXLS.
function btn_browseXLS_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browseXLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uigetfile('.xls');
set(handles.edit_listXLS,'string',[path,file]);
[~,~,ext] =fileparts(fullfile([path,file]));
if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
    [num, txt, raw] = xlsread([path,file]);
elseif strcmp(ext,'.txt')   
    [num, txt, raw] = readtxtfile_asxlsread([path,file]);
end

set(handles.edit_zonemodel,'string',fullfile(raw{2,1},raw{2,3}));

function edit_listXLS_Callback(hObject, eventdata, handles)
% hObject    handle to edit_listXLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_listXLS as text
%        str2double(get(hObject,'String')) returns contents of edit_listXLS as a double


% --- Executes during object creation, after setting all properties.
function edit_listXLS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_listXLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_zscore.
function btn_zscore_Callback(hObject, eventdata, handles)
% hObject    handle to btn_zscore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlslistfile = get(handles.edit_listXLS,'string')
[~,~,ext] =fileparts([xlslistfile]);
if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
    [raw, txt, info]=xlsread([xlslistfile]);
elseif strcmp(ext,'.txt')   
    [raw, txt, info] = readtxtfile_asxlsread([xlslistfile]);
end

groupeall = [];
for isubject=2:size(info,1)
    id = isubject-1;
    MAT = load(fullfile(info{isubject,1}, info{isubject,2}));
    DATA{id}.ZoneList = MAT.ZoneList;
    if isfield(MAT, 'meancorrPearsonFisher')
        DATA{id}.MAT = MAT.meancorrPearsonFisher;
    end
    if isfield(MAT, 'meancorr')
        matcorr = MAT.meancorr;
        % matcorr =  1/2*log((1+nanmean(matcorr,3))./(1-nanmean(matcorr,3)));
        DATA{id}.MAT = matcorr
    end
    
    DATA{id}.name = info{isubject,2};
    DATA{id}.MATtrial =  MAT.matcorr;
    DATA{id}.GR = info{isubject,4};
    DATA{id}.System = 'ISS'
    %if isubject==10
    load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');
    %end
    names = fieldnames(zone);
    for iname = 1:numel(names)
        eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname}]);
    end
    list_subject{id} =DATA{id}.name;
    groupeall = [groupeall; info{isubject,4}];
end
idsubject = 1:numel(groupeall)
for isubject = 1:numel(groupeall)
    MATall(:,:,isubject) = DATA{idsubject(isubject)}.MAT;
    groupid(isubject)= DATA{idsubject(isubject)}.GR
end
nbnode = size( MATall,2);
MATall = reshape(MATall,nbnode*nbnode,size( MATall,3));
g1 = find(groupid==1)
x = MATall;


%GROUPE 1 zscore
mu = nanmean(MATall(:,g1),2)
sigma= nanstd(MATall(:,g1),1,2)
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);
z = reshape(z, nbnode, nbnode,numel(groupid))
%Save Zscore
for isubject=2:size(info,1)
    id = isubject-1;
    filenameout = fullfile(info{isubject,1}, ['zscore',info{isubject,2}])
    groupid
    ZoneList =  DATA{id}.ZoneList
    matcorr = z(:,:, id)
    meancorr = z(:,:, id)
    save(filenameout,'ZoneList','matcorr','meancorr');
    info{isubject,2} = ['zscore',info{isubject,2}]
end

[pathstr, name, ext] = fileparts(xlslistfile);
[file,path]= uiputfile([pathstr, filesep,name,'_zscore',ext]);

if ismac
    % Code to run on Mac platform
    writetxtfile([path,file],info);
elseif isunix
    % Code to run on Linux platform
    xlswrite([path,file],info)
elseif ispc
    % Code to run on Windows platform
    xlswrite([path,file],info)
else
    disp('Platform not supported')
end


% --- Executes on button press in btn_Permutationttest.
function btn_Permutationttest_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Permutationttest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
xlslistfile = get(handles.edit_listXLS,'string');
[~,~,ext] =fileparts([xlslistfile]);
if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
    [raw, txt, info]=xlsread([xlslistfile]);
elseif strcmp(ext,'.txt')   
    [raw, txt, info] = readtxtfile_asxlsread([xlslistfile]);
end

testcardiac =0
if testcardiac
   if isfile('C:\data\Data_NIRS\BebeResting\connectivityMAT_COH_cardiac\Listsujet_Bydetector.xls')
       [rawcard, txtcard, infocard]=xlsread('C:\data\Data_NIRS\BebeResting\connectivityMAT_COH_cardiac\Listsujet_Bydetector.xls');
   elseif isfile('C:\data\Data_NIRS\BebeResting\connectivityMAT_COH_cardiac\Listsujet_Bydetector.txt')
       [rawcard, txtcard, infocard] = readtxtfile_asxlsread('C:\data\Data_NIRS\BebeResting\connectivityMAT_COH_cardiac\Listsujet_Bydetector.txt');
   end
    
end
groupeall = []; 
%Load the data
for isubject=2:size(info,1)
  
    id = isubject-1;
    MAT = load(fullfile(info{isubject,1},[ info{isubject,2},'.mat']));
    DATA{id}.ZoneList = MAT.ZoneList;
    if isfield(MAT, 'meancorrPearsonFisher')
        DATA{id}.MAT = MAT.meancorrPearsonFisher;
    end
    if isfield(MAT, 'meancorr')
        matcorr = MAT.meancorr;
        DATA{id}.MAT = matcorr;
        % matcorr =  1/2*log((1+nanmean(matcorr,3))./(1-nanmean(matcorr,3)));
        if get(handles.radio_viewdistribution,'value')
            figure;subplot(2,1,1);hist(matcorr(:));xlim([0 1]);title([info{isubject,2}])
            if testcardiac
                MATCartiac = load(fullfile(infocard{isubject,1}, infocard{isubject,2}));
                DATA{id}.MAT = matcorr .*double(MATCartiac.meancorr>0.2);
                subplot(2,1,2);hist(MATCartiac.meancorr(:));xlim([0 1])
                title([info{isubject,2}])
                %  figure; imagesc(matcorr);caxis([0 0.1])
            end
        end
    end
    
    DATA{id}.name = info{isubject,2};
    DATA{id}.MATtrial =  MAT.matcorr;
    DATA{id}.GR = info{isubject,4};
    DATA{id}.System = 'ISS';
    load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');
    names = fieldnames(zone);
    for iname = 1:numel(names)
        eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname},';']);
    end
    list_subject{id} =DATA{id}.name;
    groupeall = [groupeall; info{isubject,4}];
    clear MAT
end
if get(handles.popup_optionpermutation,'value')==1 %channel mode
    idsubject = 1:numel(groupeall);
    for isubject = 1:numel(groupeall)
        MATall(isubject,:,:)=DATA{isubject}.MAT;
        groupid(isubject)= DATA{idsubject(isubject)}.GR;
    end
    ZONEid = [info{end,3}];
    ZoneList =  DATA{end}.ZoneList;
elseif  get(handles.popup_optionpermutation,'value')==2       
    MATall =zeros(numel(DATA),numel(DATA{id}.zone.label),numel(DATA{id}.zone.label));
    for isubject = 1:numel(groupeall)
        List = DATA{isubject}.ZoneList;        
        for izone = 1:numel(DATA{isubject}.zone.label)
            ML = DATA{isubject}.zone.ml;
            DATA{isubject}.zone.plotLst;
            idlisti = [];
            idliststr = []; 
            chzone = DATA{isubject}.zone.plotLst{izone};
            for ichzone = 1:numel(chzone);
                ich = chzone(ichzone);
                if strcmp(DATA{isubject}.System,'ISS')
                    strDet = SDDet2strboxy_ISS(ML(ich,2));
                    strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                    idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                end
                idliststr =[idliststr,{[strDet, ' ',strSrs ]}];
                idlisti = [idlisti, idch];
            end
            for jzone = 1:numel(DATA{isubject}.zone.label)
                idlistj = [];
                chzone = DATA{isubject}.zone.plotLst{jzone};
                for ichzone = 1:numel(chzone);
                    ich = chzone(ichzone);
                    if strcmp(DATA{isubject}.System,'ISS')
                        strDet = SDDet2strboxy_ISS(ML(ich,2));
                        strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                        idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                    end
                    idlistj = [idlistj, idch];
                end 
                matROI = DATA{isubject}.MAT(idlisti,idlistj);
                id = find(matROI==0);
                if isempty(id)
                    matROI(id)=nan;
                end
                MATall(isubject,izone,jzone) = nanmean(matROI(:));
                if izone==jzone
                    matnbnanbyizone(isubject,izone)=numel(find(sum(double(isnan(matROI))) ==size(matROI,1)));
                    matnbtotchbyizone(isubject,izone) = size(matROI,1);
                end
            end
            
            
            
        end
        groupid(isubject)= DATA{isubject}.GR;
    end
    zoneuse=DATA{isubject}.zone;
    ZoneList = [];
    plottmp=[];
    plotLst = [];
    for izoneList = 1:size(MATall,2)
        MLfake(izoneList,1) = izoneList;%source
        MLfake(izoneList,2) = 1; %detecteur
        MLfake(izoneList,3) = 1;
        MLfake(izoneList,4) = 1;
        strDet = SDDet2strboxy_ISS(MLfake(izoneList,2));
        strSrs = SDPairs2strboxy_ISS(MLfake(izoneList,1));
        ZoneLabel{izoneList,1}=zoneuse.label{izoneList};
        ZoneList{izoneList,1} = [strDet,' ', strSrs];
        plottmp{izoneList} = [izoneList,1];
        plotLst{izoneList} = [izoneList];
    end
    %save zone list associate
    zone.plot = plottmp;
    zone.plotLst = plotLst;
    zone.label = ZoneLabel;
    zone.color = zone.color;
    zone.ml = MLfake;
    zone.chMAT = plotLst;
    save(fullfile(info{isubject,1},['avg', info{isubject,3}]),'zone','-mat');
    ZONEid = ['avg', info{isubject,3}];
end 
thr_p = str2num(get(handles.edit_maskthreshold,'string'));
%POUR TEST 
% MATall = rand(40,160,160) 
% groupid= [ones(20,1,1); ones(20,1,1)*2]
% idnan=find(MATall<0.25);
% if  ~isempty(idnan)
%     pb1(idnan)=nan;
% end
if get(handles.popup_optionpermutation,'value')==1 
    modelabel = 'ch';
elseif get(handles.popup_optionpermutation,'value')==2 
    modelabel = 'zoneavg';
end
if get(handles.popup_STAT2groupe, 'val') ==1 %PERMUTATION    
    ESTAD = 4;
    INDEP = 1;
    NPERM = 500;
    g1 = find(groupid==1);
    cb1 = MATall(g1,:,:);
    g2 = find(groupid==2);
    pb1 = MATall(g2,:,:);   
    if (sum(isnan(pb1(:)))+sum(isnan(cb1(:))))>0 %if presence NAN do the stat without them...
        for i=1:size(pb1,2) % row
            for j=1:size(cb1,3) %col
                tmpcb1 = cb1(:,i,j);
                tmppb1 = pb1(:,i,j);
                idnan = find(isnan(tmpcb1));
                if ~isempty(idnan)
                    tmpcb1(idnan)=[];
                end
                idnan = find(isnan(tmppb1))
                if ~isempty(idnan)
                    tmppb1(idnan)=[];
                end
                
                [FSupSup,FSupDeriv,FSupTime,pij,tij] = TestPermut2Grupos(ESTAD,INDEP,tmpcb1,tmppb1,NPERM);
                ncb1(i,j) = numel(tmpcb1);
                npb1(i,j) = numel(tmppb1);
                FUniv(i,j) = pij;
                Toij(i,j) = tij;
            end
        end
    else
        [FSupSup,FSupDeriv,FSupTime,FUniv,Toij] = TestPermut2Grupos(ESTAD,INDEP,pb1,cb1,NPERM);
        ncb1 = ones(size(cb1,2),size(cb1,3))*numel(size(cb1,1));  
        npb1 = ones(size(pb1,2),size(pb1,3))*numel(size(pb1,1));
    end
    MeanG1 = squeeze(nanmean(cb1));
    MeanG2 = squeeze(nanmean(pb1));
    
    %WRITE IN A NEW FILE
    infonew = [info(1,:);repmat(info(2,:),6,1)];    
    isubject=2;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation tstat']);
    matcorr = Toij;
    meancorr = Toij;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel,'permutation tstat'];
    infonew{isubject,3} = ZONEid;    
    isubject=3;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation 1-Pval']);
    matcorr = 1-FUniv;
    meancorr = 1-FUniv;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel,'permutation 1-Pval'];
    infonew{isubject,3} = ZONEid;    
    isubject=4;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation G1-G2']);
    matcorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
    meancorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel,'permutation G1-G2'];
    infonew{isubject,3} = ZONEid;    
    isubject=5;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation G2-G1']);
    matcorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
    meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel,'permutation G2-G1'];
    infonew{isubject,3} = ZONEid;
    isubject=6;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation G1-G2 p' sprintf('%02.0f',thr_p*100)]);
    matcorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<thr_p);
    meancorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<thr_p);
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel,'permutation G1-G2 p',sprintf('%02.0f',thr_p*100)];
    infonew{isubject,3} = ZONEid;
    isubject=7;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation G2-G1 p',sprintf('%02.0f',thr_p*100)]);
    matcorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<thr_p);
    meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<thr_p);
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel,'permutation G2-G1 p',sprintf('%02.0f',thr_p*100)];
    infonew{isubject,3} = ZONEid;    
    isubject=8;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation mean G1']);
    matcorr = real(squeeze(nanmean(cb1,1)));
    meancorr = real(squeeze(nanmean(cb1,1)));
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel,'permutation mean G1'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    isubject=9;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation mean G2']);
    matcorr = real(squeeze(nanmean(pb1,1)));
    meancorr = real(squeeze(nanmean(pb1,1)));
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel,'permutation mean G2'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;        
    isubject=10;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation N G1']);
    matcorr = ncb1;
    meancorr = ncb1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel,'permutation N G1'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;    
    isubject=11;
    filenameout = fullfile(info{2,1}, [modelabel,'permutation N G2']);
    matcorr = npb1;
    meancorr = npb1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel,'permutation N G2'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;    
    [pathstr, name, ext] = fileparts(xlslistfile);
    [file,path]= uiputfile([pathstr, filesep,name,modelabel ,'_permutation',ext]); 
if ismac
    % Code to run on Mac platform
    writetxtfile([path,file],infonew);
elseif isunix
    % Code to run on Linux platform
     xlswrite([path,file],infonew); 
elseif ispc
    % Code to run on Windows platform
    xlswrite([path,file],infonew); 
else
    disp('Platform not supported')
end
    
end
%TEST
%idx = 1
%idy = 2
%save('MATRICEIRSC.mat','MATall','groupid','ZoneLabel','xlslistfile','-mat')


if get(handles.popup_STAT2groupe, 'val') ==2 %Unpaired student ttest    
     g1 = find(groupid==1);
    cb1 = MATall(g1,:,:); 
    g2 = find(groupid==2);
    pb1 = MATall(g2,:,:);     
    for i=1:size(MATall,2)
        for j=1:size(MATall,2)
                tmpcb1 = cb1(:,i,j);
                tmppb1 = pb1(:,i,j);
                idnan = find(isnan(tmpcb1));
                if ~isempty(idnan)
                    tmpcb1(idnan)=[];
                end
                idnan = find(isnan(tmppb1));
                if ~isempty(idnan)
                    tmppb1(idnan)=[];
                end
                try
                    STATS = testt(tmpcb1,tmppb1);
                    FUniv(i,j) = STATS.tpvalue;
                    Toij(i,j) = STATS.tvalue;
                    MeanG1(i,j) = nanmean(tmpcb1);
                    MeanG2(i,j) = nanmean(tmppb1);
                    ncb1(i,j) = numel(tmpcb1);
                    npb1(i,j) = numel(tmppb1);
                catch
                    disp('error testt')
                    FUniv(i,j) = 1;
                    Toij(i,j) = 0;
                    MeanG1(i,j) = 0;
                    MeanG2(i,j) = 0;
                    ncb1(i,j) = 0;
                    npb1(i,j) =0;
                end
            %figure;hist(val)
                  
        end
    end
    [FDRmat,pcritic] = fdr_bh(FUniv, thr_p);
    %figure;cdfplot(FUniv(:))
        %WRITE IN A NEW FILE
    infonew = [info(1,:);repmat(info(2,:),6,1)];    
    isubject=2;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired ttest tval']);
    matcorr = Toij;
    meancorr = Toij;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel, 'unpaired ttest tval'];
    infonew{isubject,3} = ZONEid;
    isubject=3;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired 1-Pval']);
    matcorr = 1-FUniv;
    meancorr = 1-FUniv;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel, 'unpaired 1-Pval'];
    infonew{isubject,3} = ZONEid;
%     [pID, p_masked] = fdr(FUniv, 0.05);
%     figure;imagesc(p_masked);
    isubject=4;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired G1-G2']);
    matcorr = MeanG1-MeanG2;
    meancorr =  MeanG1-MeanG2;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel, 'unpaired G1-G2'];
    infonew{isubject,3} = ZONEid;
    
    isubject=5;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired G2-G1']);
    matcorr = MeanG2-MeanG1;
    meancorr = MeanG2-MeanG1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel, 'unpaired G2-G1'];
    infonew{isubject,3} = ZONEid;
    
    isubject=6;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired G1-G2 p' sprintf('%02.0f',thr_p*100)]);
    matcorr = (MeanG1-MeanG2).*double(FUniv<thr_p);
    meancorr = (MeanG1-MeanG2).*double(FUniv<thr_p);
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = [modelabel, 'unpaired G1-G2 p',sprintf('%02.0f',thr_p*100)];
    infonew{isubject,3} = ZONEid;
    
    isubject=7;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired G2-G1 p',sprintf('%02.0f',thr_p*100)]);
    matcorr =  (MeanG2-MeanG1).*double(FUniv<thr_p);
    meancorr = (MeanG2-MeanG1).*double(FUniv<thr_p);
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired G2-G1 p',sprintf('%02.0f',thr_p*100)];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    
    isubject=8;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired G1-G2 qFDR' sprintf('%02.0f',thr_p*100)]);
    matcorr = (MeanG1-MeanG2).*double(FDRmat);
    meancorr = (MeanG1-MeanG2).*double(FDRmat);
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired G1-G2 qFDR',sprintf('%02.0f',thr_p*100)];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    isubject=9;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired G2-G1 qFDR',sprintf('%02.0f',thr_p*100)]);
    matcorr =  (MeanG2-MeanG1).*double(FDRmat);
    meancorr = (MeanG2-MeanG1).*double(FDRmat);
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired G2-G1 qFDR',sprintf('%02.0f',thr_p*100)];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    
    isubject=10;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired mean G1']);
    matcorr = MeanG1;
    meancorr = MeanG1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired mean G1'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    isubject=11;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired mean G2']);
    matcorr = MeanG2;
    meancorr = MeanG2;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired mean G2'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1; 
    isubject=12;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired N G1']);
    matcorr = ncb1;
    meancorr = ncb1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired N G1'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    isubject=13;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired N G2']);
    matcorr = npb1;
    meancorr = npb1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired N G2'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;    
    
    isubject=14;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired N G1 POURCENTAGE']);
    matcorr = ncb1./numel(g1)*100;
    meancorr = ncb1./numel(g1)*100;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired N G1 POURCENTAGE'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    isubject=15;
    filenameout = fullfile(info{2,1}, [modelabel, 'unpaired N G2 POURCENTAGE']);
    matcorr = npb1./numel(g2)*100;
    meancorr = npb1./numel(g2)*100;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = [modelabel, 'unpaired N G2 POURCENTAGE'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;    
    [pathstr, name, ext] = fileparts(xlslistfile);
    [file,path]= uiputfile([pathstr, filesep, name,modelabel,'_unpairedtest',ext]);
    if ismac        
        % Code to run on Mac platform
        writetxtfile([path,file],infonew);       
    elseif isunix
        % Code to run on Linux platform
        xlswrite([path,file],infonew);   
    elseif ispc
        % Code to run on Windows platform
        xlswrite([path,file],infonew);    
    else
        disp('Platform not supported')
    end
    
    
    
    
end

if get(handles.popup_STAT2groupe, 'val') ==3 %Wilcoxon Rank sum    
    g1 = find(groupid==1);
    cb1 = MATall(g1,:,:);
    g2 = find(groupid==2);
    pb1 = MATall(g2,:,:);    
    
   
    for i=1:size(MATall,2)
        for j=1:size(MATall,2)
                tmpcb1 = cb1(:,i,j);
                tmppb1 = pb1(:,i,j);
                idnan = find(isnan(tmpcb1));
                if ~isempty(idnan)
                    tmpcb1(idnan)=[];
                end
                idnan = find(isnan(tmppb1));
                if ~isempty(idnan)
                    tmppb1(idnan)=[];
                end
            [pval,h,STATS] = ranksum(tmpcb1,tmppb1);
            %figure;hist(val)
            FUniv(i,j) = pval;
            Toij(i,j) = STATS.ranksum;
            MeanG1(i,j) = nanmean(tmpcb1);
            MeanG2(i,j) = nanmean(tmppb1);
            ncb1(i,j) = numel(tmpcb1);
            npb1(i,j) = numel(tmppb1);           
        end
    end
        %WRITE IN A NEW FILE
    infonew = [info(1,:);repmat(info(2,:),6,1)];    
    isubject=2;
    filenameout = fullfile(info{2,1}, ['Wilcoxon ranksum']);
    matcorr = Toij;
    meancorr = Toij;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = ['Wilcoxon ttest tval'];
    infonew{isubject,3} = ZONEid;
    
    isubject=3;
    filenameout = fullfile(info{2,1}, ['1-Pval']);
    matcorr = 1-FUniv;
    meancorr = 1-FUniv;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = ['1-Pval'];
    infonew{isubject,3} = ZONEid;
    
    isubject=4;
    filenameout = fullfile(info{2,1}, ['G1-G2']);
    matcorr = MeanG1-MeanG2;
    meancorr =  MeanG1-MeanG2;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = ['G1-G2'];
    infonew{isubject,3} = ZONEid;
    
    isubject=5;
    filenameout = fullfile(info{2,1}, ['G2-G1']);
    matcorr = MeanG2-MeanG1;
    meancorr = MeanG2-MeanG1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = ['G2-G1'];
    infonew{isubject,3} = ZONEid;
    
    isubject=6;
    filenameout = fullfile(info{2,1}, ['G1-G2 p' sprintf('%02.0f',thr_p*100)]);
    matcorr = (MeanG1-MeanG2).*double(FUniv<thr_p);
    meancorr = (MeanG1-MeanG2).*double(FUniv<thr_p);
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = ['G1-G2 p',sprintf('%02.0f',thr_p*100)];
    infonew{isubject,3} = ZONEid;
    
    isubject=7;
    filenameout = fullfile(info{2,1}, ['G2-G1 p',sprintf('%02.0f',thr_p*100)]);
    matcorr =  (MeanG2-MeanG1).*double(FUniv<thr_p);
    meancorr = (MeanG2-MeanG1).*double(FUniv<thr_p);
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,2} = ['G2-G1 p',sprintf('%02.0f',thr_p*100)];
    infonew{isubject,3} = ZONEid;
    
    isubject=8;
    filenameout = fullfile(info{2,1}, ['mean G1']);
    matcorr = MeanG1;
    meancorr = MeanG1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = ['mean G1'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    isubject=9;
    filenameout = fullfile(info{2,1}, ['mean G2']);
    matcorr = MeanG2;
    meancorr = MeanG2;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = ['mean G2'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1; 
    isubject=10;
    filenameout = fullfile(info{2,1}, ['N G1']);
    matcorr = ncb1;
    meancorr = ncb1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = ['N G1'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;
    isubject=11;
    filenameout = fullfile(info{2,1}, ['N G2']);
    matcorr = npb1;
    meancorr = npb1;
    save(filenameout,'ZoneList','matcorr','meancorr');
    infonew{isubject,1} = info{2,1};
    infonew{isubject,2} = ['N G2'];
    infonew{isubject,3} = ZONEid;
    infonew{isubject,4} = 1;    
    [pathstr, name, ext] = fileparts(xlslistfile);
    [file,path]= uiputfile([pathstr, filesep,name,modelabel,'_Wilcoxon',ext]);
    
    
    if ismac        
        % Code to run on Mac platform
        writetxtfile([path,file],infonew);       
    elseif isunix
        % Code to run on Linux platform
        xlswrite([path,file],infonew);   
    elseif ispc
        % Code to run on Windows platform
        xlswrite([path,file],infonew);    
    else
        disp('Platform not supported')
    end
    
end


%NBCHANNEL BY ZONE
try
    cellnbsujet = txt;
    linedescription = [DATA{isubject}.zone.label; DATA{isubject}.zone.label];
    linedescription = linedescription(:)';
    val = zeros(size(matnbtotchbyizone,1),size(linedescription,2));
    val(:,1:2:end)=matnbnanbyizone;
    val(:,2:2:end)=matnbtotchbyizone;
    
    cellnbsujet = [cellnbsujet,[linedescription;num2cell(val) ]];
    [file,path]= uiputfile([cellnbsujet{2,1},'nbchannebyzone.xls'])
    if ismac        
        % Code to run on Mac platform
        writetxtfile([path,file],cellnbsujet);       
    elseif isunix
        % Code to run on Linux platform
        xlswrite([path,file],cellnbsujet);   
    elseif ispc
        % Code to run on Windows platform
        xlswrite([path,file],cellnbsujet);    
    else
        disp('Platform not supported')
    end
    statst =testt(matnbnanbyizone(:,1), matnbnanbyizone(:,2))
catch
end





% --- Executes on button press in btn_NBS.
function btn_NBS_Callback(hObject, eventdata, handles)
% hObject    handle to btn_NBS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlslistfile = get(handles.edit_listXLS,'string')
[~,~,ext] =fileparts(xlslistfile);
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        [raw, txt, info]=xlsread([xlslistfile]);
    elseif strcmp(ext,'.txt')   
        [raw, txt, info] = readtxtfile_asxlsread([xlslistfile]);
    end
groupeall = [];
for isubject=2:size(info,1)
    id = isubject-1;
    MAT = load(fullfile(info{isubject,1}, info{isubject,2}));
    DATA{id}.ZoneList = MAT.ZoneList;
    if isfield(MAT, 'meancorrPearsonFisher')
        DATA{id}.MAT = MAT.meancorrPearsonFisher;
    end
    if isfield(MAT, 'meancorr')
        matcorr = MAT.meancorr;
        matcorr =  1/2*log((1+nanmean(matcorr,3))./(1-nanmean(matcorr,3)));
        DATA{id}.MAT = matcorr
    end
    
    DATA{id}.name = info{isubject,2};
    DATA{id}.MATtrial =  MAT.matcorr;
    DATA{id}.GR = info{isubject,4};
    DATA{id}.System = 'ISS'
    %if isubject==10
    load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');
    %end
    names = fieldnames(zone);
    for iname = 1:numel(names)
        eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname},';']);
    end
    list_subject{id} =DATA{id}.name;
    groupeall = [groupeall; info{isubject,4}];
end
idsubject = 1:numel(groupeall)


for isubject = 1:numel(groupeall)
    MATall(:,:,isubject) = DATA{idsubject(isubject)}.MAT;
    groupid(isubject)= DATA{idsubject(isubject)}.GR
end
mkdir([info{isubject,1},'/NBS'])
Mat = MATall
save(fullfile([info{isubject,1},filesep,'NBS'],'DATA.mat'),'Mat')
design= double([ groupeall==1,groupeall==2])
save(fullfile([info{isubject,1},filesep,'NBS'],'design.mat'),'design')
contrast= [ 1,-1]
save(fullfile([info{isubject,1},filesep,'NBS'],'G1higherG2contrast.mat'),'contrast','-mat')

listelectrode = DATA{1}.ZoneList
fid = fopen(fullfile([info{isubject,1},filesep,'NBS'],'nodeLabels.txt'),'w')
for i=1:numel(listelectrode)
    fprintf(fid,'%s',listelectrode{i})
    if i<numel(listelectrode)
        fprintf(fid,'\r')
    end
end
fclose(fid)

%OUTOUT POUR NBSToolbox not realisitic position... to be built if needed
fid = fopen(fullfile([info{isubject,1},filesep,'NBS'],'COG.txt'),'w')
for i=1:numel(listelectrode)
    fprintf(fid,'%2.3f %2.3f %2.3f\r',0,0,0 )
end
fclose(fid)




% --- Executes on button press in btn_ZoneCH.
function btn_ZoneCH_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ZoneCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xlslistfile = get(handles.edit_listXLS,'string')
[~,~,ext] =fileparts(fileevent);
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        [raw, txt, info]=xlsread([xlslistfile]);
    elseif strcmp(ext,'.txt')   
        [raw, txt, info] = readtxtfile_asxlsread([xlslistfile]);
    end

groupeall = [];
for isubject=2:size(info,1)
    id = isubject-1;
    MAT = load(fullfile(info{isubject,1}, info{isubject,2}));
    DATA{id}.ZoneList = MAT.ZoneList;
    if isfield(MAT, 'meancorrPearsonFisher')
        DATA{id}.MAT = MAT.meancorrPearsonFisher;
    end
    if isfield(MAT, 'meancorr')
        matcorr = MAT.meancorr;
        matcorr =  1/2*log((1+nanmean(matcorr,3))./(1-nanmean(matcorr,3)));
        DATA{id}.MAT = matcorr
    end
    
    DATA{id}.name = info{isubject,2};
    DATA{id}.MATtrial =  MAT.matcorr;
    DATA{id}.GR = info{isubject,4};
    DATA{id}.System = 'ISS'
    %if isubject==10
    load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');
    %end
    names = fieldnames(zone);
    for iname = 1:numel(names)
        eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname}]);
    end
    list_subject{id} =DATA{id}.name;
    groupeall = [groupeall; info{isubject,4}];
    groupid(id)=  DATA{id}.GR;
    
end
%%%
idsubject = 1:numel(DATA)
%zone seed stats
idlabelall= DATA{1}.zone.label; %zone premier sujet du groupe
MATAVGall = zeros(numel(idsubject),numel(idlabelall),numel(idlabelall));

for isubject=1:numel(idsubject)
    ML=DATA{idsubject(isubject)}.zone.ml;
    List=strvcat(DATA{idsubject(isubject)}.ZoneList);
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
                strDet = SDDet2strboxy_ISS(ML(ich,2));
                strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                idch = strmatch([strDet, ' ',strSrs ],List,'exact');
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
                strDet = SDDet2strboxy_ISS(ML(ich,2));
                strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                idch = strmatch([strDet, ' ',strSrs ],List,'exact');
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
    MATAVGall(isubject,:,:) = MATAVG;
end

ESTAD = 4;
INDEP = 1;
NPERM = 500;
g1 = find(groupid==1);
cb1 =  MATAVGall(g1,:,:);
g2 = find(groupid==2);
pb1 =  MATAVGall(g2,:,:)
[FSupSup,FSupDeriv,FSupTime,FUniv,Toij] = TestPermut2Grupos(ESTAD,INDEP,pb1,cb1,NPERM);


figure
imagesc(Toij.*FUniv<0.05)

ZoneList = [];
plot= [];
plotLst = [];
for izoneList = 1:size(MATAVGall,2)
    MLfake(izoneList,1) = izoneList;%source
    MLfake(izoneList,2) = 1; %detecteur
    MLfake(izoneList,3) = 1;
    MLfake(izoneList,4) = 1;
    strDet = SDDet2strboxy_ISS(MLfake(izoneList,2));
    strSrs = SDPairs2strboxy_ISS(MLfake(izoneList,1));
    ZoneList{izoneList,1}=[strDet,' ', strSrs]
    plottmp{izoneList} = [izoneList,1];
    plotLst{izoneList} = [izoneList];
end

%save zone list associate
zone.plot = plottmp;
zone.plotLst = plotLst;
zone.label = idlabelall;
zone.color = DATA{idsubject(1)}.zone.color;
zone.ml = MLfake;
zone.chMAT = plotLst;
save(fullfile(info{isubject,1},'zonecombine.zone'),'zone','-mat')
%  zone.chMAT =
%save .mat associate
% ZoneList, matcorr, meancorr
for isubject=2:size(info,1)
    id = isubject-1;
    filenameout = fullfile(info{isubject,1}, ['zone',info{isubject,2}]);
    groupid;
    ZoneList =  DATA{id}.ZoneList;
    matcorr = squeeze(MATAVGall(id,:,:));
    meancorr =  squeeze(MATAVGall(id,:,:));
    save(filenameout,'ZoneList','matcorr','meancorr');
    info{isubject,2} = ['zone',info{isubject,2}]
    info{isubject,3} = 'zonecombine.zone'
end
info(isubject+1,:)  = info(isubject,:)
info{isubject+1,2} = ['pval']
ZoneList =  DATA{id}.ZoneList;
filenameout = fullfile(info{isubject,1}, ['pval'])
matcorr = FUniv;
meancorr = FUniv;
save(filenameout,'ZoneList','matcorr','meancorr');


[file,path]= uiputfile(xlslistfile)

if ismac
    % Code to run on Mac platform
    writetxtfile([path,file],info);
elseif isunix
    % Code to run on Linux platform
    xlswrite([path,file],info);
elseif ispc
    % Code to run on Windows platform
    xlswrite([path,file],info);
else
    disp('Platform not supported')
end




function edit_maskthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maskthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maskthreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_maskthreshold as a double


% --- Executes during object creation, after setting all properties.
function edit_maskthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maskthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_permutationLevel.
function popup_permutationLevel_Callback(hObject, eventdata, handles)
% hObject    handle to popup_permutationLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_permutationLevel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_permutationLevel


% --- Executes during object creation, after setting all properties.
function popup_permutationLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_permutationLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_optionpermutation.
function popup_optionpermutation_Callback(hObject, eventdata, handles)
% hObject    handle to popup_optionpermutation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_optionpermutation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_optionpermutation
if get(handles.popup_optionpermutation,'val')==1
    set(handles.edit_zonemodel,'enable','off')
    set(handles.radio_nbvalidchannel,'enable','off')
elseif get(handles.popup_optionpermutation,'val')==2
    set(handles.edit_zonemodel,'enable','on')
    set(handles.radio_nbvalidchannel,'enable','on')
end


% --- Executes during object creation, after setting all properties.
function popup_optionpermutation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_optionpermutation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_viewdistribution.
function radio_viewdistribution_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewdistribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_viewdistribution


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_zonemodel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zonemodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zonemodel as text
%        str2double(get(hObject,'String')) returns contents of edit_zonemodel as a double


% --- Executes during object creation, after setting all properties.
function edit_zonemodel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zonemodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_nbvalidchannel.
function radio_nbvalidchannel_Callback(hObject, eventdata, handles)
% hObject    handle to radio_nbvalidchannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_nbvalidchannel


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on selection change in popup_STAT2groupe.
function popup_STAT2groupe_Callback(hObject, eventdata, handles)
% hObject    handle to popup_STAT2groupe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_STAT2groupe contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_STAT2groupe


% --- Executes during object creation, after setting all properties.
function popup_STAT2groupe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_STAT2groupe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over popup_STAT2groupe.
function popup_STAT2groupe_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to popup_STAT2groupe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.popup_STAT2groupe,'value')==2
    %FDR
end 

% fdr_bh() - Executes the Benjamini & Hochberg (1995) and the Benjamini &
%            Yekutieli (2001) procedure for controlling the false discovery 
%            rate (FDR) of a family of hypothesis tests. FDR is the expected
%            proportion of rejected hypotheses that are mistakenly rejected 
%            (i.e., the null hypothesis is actually true for those tests). 
%            FDR is a somewhat less conservative/more powerful method for 
%            correcting for multiple comparisons than procedures like Bonferroni
%            correction that provide strong control of the family-wise
%            error rate (i.e., the probability that one or more null
%            hypotheses are mistakenly rejected).
%
%            This function also returns the false coverage-statement rate 
%            (FCR)-adjusted selected confidence interval coverage (i.e.,
%            the coverage needed to construct multiple comparison corrected
%            confidence intervals that correspond to the FDR-adjusted p-values).
%
%
% Usage:
%  >> [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., Gaussian variables that are positively correlated or
%             independent).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure (e.g.,
%             Gaussian variables with any covariance matrix) is used. 'dep'
%             is always appropriate to use but is less powerful than 'pdep.'
%             {default: 'pdep'}
%   report  - ['yes' or 'no'] If 'yes', a brief summary of FDR results are
%             output to the MATLAB command line {default: 'no'}
%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the 
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All uncorrected p-values less than or equal to crit_p are 
%             significant (i.e., their null hypotheses are rejected).  If 
%             no p-values are significant, crit_p=0.
%   adj_ci_cvrg - The FCR-adjusted BH- or BY-selected 
%             confidence interval coverage. For any p-values that 
%             are significant after FDR adjustment, this gives you the
%             proportion of coverage (e.g., 0.99) you should use when generating
%             confidence intervals for those parameters. In other words,
%             this allows you to correct your confidence intervals for
%             multiple comparisons. You can NOT obtain confidence intervals 
%             for non-significant p-values. The adjusted confidence intervals
%             guarantee that the expected FCR is less than or equal to q
%             if using the appropriate FDR control algorithm for the  
%             dependency structure of your data (Benjamini & Yekutieli, 2005).
%             FCR (i.e., false coverage-statement rate) is the proportion 
%             of confidence intervals you construct
%             that miss the true value of the parameter. adj_ci=NaN if no
%             p-values are significant after adjustment.
%   adj_p   - All adjusted p-values less than or equal to q are significant
%             (i.e., their null hypotheses are rejected). Note, adjusted 
%             p-values can be greater than 1.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
%   Benjamini, Y., & Yekutieli, D. (2005). False discovery rate?adjusted 
%     multiple confidence intervals for selected parameters. Journal of the 
%     American Statistical Association, 100(469), 71?81. doi:10.1198/016214504000001907
%
%
% Example:
%  nullVars=randn(12,15);
%  [~, p_null]=ttest(nullVars); %15 tests where the null hypothesis
%  %is true
%  effectVars=randn(12,5)+1;
%  [~, p_effect]=ttest(effectVars); %5 tests where the null
%  %hypothesis is false
%  [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_null p_effect],.05,'pdep','yes');
%  data=[nullVars effectVars];
%  fcr_adj_cis=NaN*zeros(2,20); %initialize confidence interval bounds to NaN
%  if ~isnan(adj_ci_cvrg),
%     sigIds=find(h);
%     fcr_adj_cis(:,sigIds)=tCIs(data(:,sigIds),adj_ci_cvrg); % tCIs.m is available on the
%     %Mathworks File Exchagne
%  end
%
%
% For a review of false discovery rate control and other contemporary
% techniques for correcting for multiple comparisons see:
%
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis 
% of event-related brain potentials/fields I: A critical tutorial review. 
% Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x 
% http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint1.pdf
%
%
% For a review of FCR-adjusted confidence intervals (CIs) and other techniques 
% for adjusting CIs for multiple comparisons see:
%
%   Groppe, D.M. (in press) Combating the scientific decline effect with 
% confidence (intervals). Psychophysiology.
% http://biorxiv.org/content/biorxiv/early/2015/12/10/034074.full.pdf
%
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 5/7/2010-Added FDR adjusted p-values
% 5/14/2013- D.H.J. Poot, Erasmus MC, improved run-time complexity
% 10/2015- Now returns FCR adjusted confidence intervals

function [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report)
% Copyright (c) 2015, David Groppe
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvals>1,1)),
        error('Some p-values are greater than 1.');
    end
end

if nargin<2,
    q=.05;
end

if nargin<3,
    method='pdep';
end

if nargin<4,
    report='no';
end

s=size(pvals);
if (length(s)>2) || s(1)>1,
    [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    [p_sorted, sort_ids]=sort(pvals);
end
[dummy, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests

if strcmpi(method,'pdep'),
    %BH procedure for independence or positive dependence
    thresh=(1:m)*q/m;
    wtd_p=m*p_sorted./(1:m);
    
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom=m*sum(1./(1:m));
    thresh=(1:m)*q/denom;
    wtd_p=denom*p_sorted./[1:m];
    %Note, it can produce adjusted p-values greater than 1!
    %compute adjusted p-values
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end

if nargout>3,
    %compute adjusted p-values; This can be a bit computationally intensive
    adj_p=zeros(1,m)*NaN;
    [wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
    nextfill = 1;
    for k = 1 : m
        if wtd_p_sindex(k)>=nextfill
            adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
            nextfill = wtd_p_sindex(k)+1;
            if nextfill>m
                break;
            end;
        end;
    end;
    adj_p=reshape(adj_p(unsort_ids),s);
end

rej=p_sorted<=thresh;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id),
    crit_p=0;
    h=pvals*0;
    adj_ci_cvrg=NaN;
else
    crit_p=p_sorted(max_id);
    h=pvals<=crit_p;
    adj_ci_cvrg=1-thresh(max_id);
end

if strcmpi(report,'yes'),
    n_sig=sum(p_sorted<=crit_p);
    if n_sig==1,
        fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n',m,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n',m,n_sig,q);
    end
    if strcmpi(method,'pdep'),
        fprintf('FDR/FCR procedure used is guaranteed valid for independent or positively dependent tests.\n');
    else
        fprintf('FDR/FCR procedure used is guaranteed valid for independent or dependent tests.\n');
    end
end
