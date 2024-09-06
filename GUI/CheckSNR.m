function varargout = CheckSNR(varargin)
% CHECKSNR M-file for CheckSNR.fig
%      CHECKSNR, by itself, creates a new CHECKSNR or raises the existing
%      singleton*.
%
%      H = CHECKSNR returns the handle to a new CHECKSNR or the handle to
%      the existing singleton*.
%
%      CHECKSNR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHECKSNR.M with the given input arguments.
%
%      CHECKSNR('Property','Value',...) creates a new CHECKSNR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CheckSNR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CheckSNR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheckSNR

% Last Modified by GUIDE v2.5 18-Dec-2018 16:36:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CheckSNR_OpeningFcn, ...
                   'gui_OutputFcn',  @CheckSNR_OutputFcn, ...
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


% --- Executes just before CheckSNR is made visible.
function CheckSNR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CheckSNR (see VARARGIN)

% Choose default command line output for CheckSNR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CheckSNR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CheckSNR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in btn_SNRstd.
function btn_SNRstd_Callback(hObject, eventdata, handles)
% hObject    handle to btn_SNRstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
config{1}.NIRSini = get(handles.edit_Ini,'string')
config{1}.NIRSPCA = get(handles.edit_PCA,'string');
config{1}.NIRSparafac =get(handles.edit_PARAFAC,'string')
config{1}.NIRSICA =  get(handles.edit_ICA,'string');
[pathhome,file,ext]= fileparts(config{1}.NIRSPCA)
config{1}.CorrectionApply = [pathhome,'\CorrectionApply'];
[pathhome,file,ext]= fileparts(config{1}.NIRSini)
showfigure  =get(handles.radio_stdFigure,'value');
id = 1;
[num, txt, raw] =xlsread(get(handles.edit_baselinefile,'string'))
config{1}.baseline = num;
for iconfig =1
    tstart = config{iconfig}.baseline(:,1);
    tstop  = config{iconfig}.baseline(:,2);
    load(config{iconfig}.CorrectionApply);% PARCORR  structure toute les
    
    
   for ievent = 1:numel(PARCORR)
    
    bloc = PARCORR(ievent).file ;
    chlst =[];
    %open DATAini (non corriger) DATApar(parafac) DATAica(ICA correctin) 
    
    load(config{iconfig}.NIRSini);
    NC = NIRS.Cf.H.C.N;
    xname = NIRS.Dt.fir.pp(end).p{bloc};
    %REPLACE BEGINING PATH by acual NIRS.mat localisation
    [xnamepath, xnamefile, xnameext]= fileparts(xname);
    [xpath,xfile,xext] = fileparts(config{iconfig}.NIRSini);
    xname= fullfile(xpath,[xnamefile, xnameext]);
    DATAini = fopen_NIR(xname,NC)';
    tHRF = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(DATAini,1)*1/NIRS.Cf.dev.fs;


    % open PARAFAC
    load(config{iconfig}.NIRSparafac);
    NC = NIRS.Cf.H.C.N;
    xname = NIRS.Dt.fir.pp(end).p{bloc};
    %REPLACE BEGINING PATH by acual NIRS.mat localisation
    [xnamepath, xnamefile, xnameext]= fileparts(xname);
    [xpath,xfile,xext] = fileparts(config{iconfig}.NIRSparafac);
    xname= fullfile(xpath,[xnamefile, xnameext]);
    DATApar = fopen_NIR(xname,NC)';

    load(config{iconfig}.NIRSICA);
    NC = NIRS.Cf.H.C.N;
    xname = NIRS.Dt.fir.pp(end).p{bloc};
    [xnamepath, xnamefile, xnameext]= fileparts(xname);
    [xpath,xfile,xext] = fileparts(config{iconfig}.NIRSICA);
    xname= fullfile(xpath,[xnamefile, xnameext]);
    DATAica = fopen_NIR(xname,NC)';
    
    load(config{iconfig}.NIRSPCA);
    NC = NIRS.Cf.H.C.N;
    xname = NIRS.Dt.fir.pp(end).p{bloc};
    [xnamepath, xnamefile, xnameext]= fileparts(xname);
    [xpath,xfile,xext] = fileparts(config{iconfig}.NIRSPCA);
    xname= fullfile(xpath,[xnamefile, xnameext]);
    DATApca = fopen_NIR(xname,NC)';
    
%period artefact
indt = PARCORR(ievent).indt;
lstch = PARCORR(ievent).listgood;
 HbRlist = find(lstch > NC/2);

 if ~isempty(HbRlist)
  if sum(lstch(HbRlist) - NC/2 == lstch(1:end/2))==NC/2
     msgbox('problem HbO HbR problem')
  end
  lstch(HbRlist) = [];

 end
if showfigure
figure
hold on
plot(tHRF(indt),DATAini(indt,lstch),'k','displayname','ini')
plot(tHRF(indt), DATApar(indt,lstch),'m','displayname','PARAFAC')
plot(tHRF(indt),DATAica(indt,lstch),'g','displayname','ICA')
plot(tHRF(indt), DATApca(indt,lstch),'c','displayname','pca')
figure
hold on
plot(tHRF(:),DATAini(:,lstch),'k','displayname','ini')
plot(tHRF(:), DATApar(:,lstch),'m','displayname','PARAFAC')
plot(tHRF(:),DATAica(:,lstch),'g','displayname','ICA')
plot(tHRF(:), DATApca(:,lstch),'c','displayname','pca')
end
chHbO = lstch;
chHbR = lstch + NC/2;
Xm = PARCORR(ievent).Xm;
    label = {'sujet','bloc','CH','WAVELENGT','sample',' ',' ', 'INI','PCA'};
idRHO = 1
for ich = 1:numel(lstch)
    timeco = indt;
    timeok= find(tHRF>tstart(bloc)& tHRF<tstop(bloc));
    ch = chHbO(ich);
    SNR(id,1) = iconfig;      %SUJET  
    SNR(id,2) = bloc ;        %BLOC 
    SNR(id,3) = lstch(ich); %CH
    SNR(id,4) =  1;          %1 HbO 2 HbR
    SNR(id,5) =  numel(timeco);   
    
 
    SNR(id,8) = log10(nanvar(DATAini(timeok,ch))/nanvar(DATAini(timeco,ch)));     %INI
    SNR(id,9) = log10(nanvar(DATAini(timeok,ch))/nanvar(DATApca(timeco,ch)));    %PCA
    SNR(id,10) = log10(nanvar(DATAini(timeok,ch))/nanvar(DATApar(timeco,ch)));     %PARAFAC
    SNR(id,11) = log10(nanvar(DATAini(timeok,ch))/nanvar(DATAica(timeco,ch)));    %ICA 
    id = id+ 1
    
    ch1 = chHbO(ich);
    ch2 = chHbR(ich);

    [RHO(idRHO,1),PVAL] = corr(DATAini(timeco,ch1), DATAini(timeco,ch2));
    [RHO(idRHO,2),PVAL] = corr(DATApca(timeco,ch1), DATApca(timeco,ch2));
    [RHO(idRHO,3),PVAL] = corr(DATApar(timeco,ch1), DATApar(timeco,ch2));
    [RHO(idRHO,4),PVAL] = corr(DATAica(timeco,ch1), DATAica(timeco,ch2));
    idRHO = idRHO+1;
    %    [RHO(id,5),PVAL] = corr(DATAini(timeok,ch1),DATAini(timeok,ch2));
    %HBR    
    ch = chHbR(ich);
    SNR(id,1) = iconfig;      %SUJET  
    SNR(id,2) = bloc ;        %BLOC
    SNR(id,3) =  lstch(ich); %CH
    SNR(id,4) =  2;          %1 HbO 2 HbR
    SNR(id,5) =  numel(timeco);       
    SNR(id,8) = log10(nanvar(DATAini(timeok,ch))/nanvar(DATAini(timeco,ch)));     %INI
    SNR(id,9) = log10(nanvar(DATAini(timeok,ch))/nanvar(DATApca(timeco,ch)));     %PCA
    SNR(id,10) = log10(nanvar(DATAini(timeok,ch))/nanvar(DATApar(timeco,ch)));    %PARAFAC
    SNR(id,11) = log10(nanvar(DATAini(timeok,ch))/nanvar(DATAica(timeco,ch)));    %ICA 
    
    SNR(id,13) = RHO(id,1);
    SNR(id,14) = RHO(id,2);
    SNR(id,15) = RHO(id,3);
    SNR(id,16) = RHO(id,4);
    
    id = id+ 1;
    if 0 %numel(timeco)>numel(timeok)
        timecocut = timeco(1:numel(timeok));
        Fs = NIRS.Cf.dev.fs;
            f = Fs/2*linspace(0,1,NFFT/2+1);

        y = DATAini(timeok,ch);
        ycorr = DATAini(timecocut,ch);
        ypar = DATApar(timecocut,ch);
        ypca = DATApca(timecocut,ch);
        yica =   DATAica(timecocut,ch);
        L = numel(timeok);
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        %Y = fft(DATAini(timeok,:)-ones(size(DATAini(timeok,:),1),1)*DATAini(1,:),NFFT)/L;
        Y = fft(y,NFFT)/L;
        Ycorr = fft(ycorr,NFFT)/L;
        Ypar = fft(ypar,NFFT)/L;
        Ypca = fft(ypca,NFFT)/L;     
         Yica = fft(yica,NFFT)/L;  
        fstart = 5;
        PWD =  abs(Y(fstart:NFFT/2+1,:));
        PWDcorr = abs(Ycorr(fstart:NFFT/2+1,:));
        PWDpar =  abs(Ypar(fstart:NFFT/2+1,:));
        PWDpca =  abs(Ypca(fstart:NFFT/2+1,:));
        PWDica =  abs(Yica(fstart:NFFT/2+1,:));

        fPWD = f(fstart:end);
%         figure;hold on
%         plot(fPWD,PWD,'x')
%         plot(fPWD,PWDcorr,'xr')
%         plot(fPWD,PWDpar,'xg')
%         plot(fPWD,PWDpca,'xk')
        idfft = idfft+1;
        
          SNRPWD( idfft ,1) =  mean(PWD)/mean(PWDcorr);
          SNRPWD( idfft ,2) =  mean(PWD)/mean(PWDpar) ;         
          SNRPWD( idfft ,3) =  mean(PWD)/mean(PWDica);
          SNRPWD( idfft ,4) =  mean(PWD)/mean(PWDpca);
         
    end
        end
    end
end 

idkeep = find(~isnan(sum(SNR(:,8:11),2))&  ~isinf(sum(SNR(:,8:11),2)));
X = SNR(idkeep,8:11);
[p,tbl,stats] = anova1(X)
set(gca,'XTickLabel',{'SNRini','PCA','PARAFAC', 'ICA'},'fontsize',20)
set(gca,'XTick',[1,2,3,4],'fontsize',20)
[c,m] = multcompare(stats)
set(gca,'yTickLabel',{'SNRICA','Parafac','PCA','RAW'});
set(gca,'yTick',[1,2,3,4]);


X = RHO(2:2:end,:);
[p,tbl,stats] = anova1(X)
set(gca,'XTickLabel',{'R2ini','PCA','PARAFAC', 'ICA'},'fontsize',20)
set(gca,'XTick',[1,2,3,4],'fontsize',20)
[c,m] = multcompare(stats)
set(gca,'yTickLabel',{'R2ICA','Parafac','PCA','RAW'});
set(gca,'yTick',[1,2,3,4]);

[file,path]= uiputfile([pathhome,'\SNRSTD.xls'])
label = {'sujet','bloc','CH','WAVELENGT','sample',' ','SNR var', 'INI','PCA','PARAFAC','ICA','R2wavelenght','INI','RhoPCA','RhoPARAFAC','RhoICA'};
tmp = [label; num2cell(SNR)];
xlswrite([path,file],tmp);


% --- Executes on button press in btn_Wavelet.
function btn_Wavelet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Wavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
config{1}.NIRSini = get(handles.edit_Ini,'string')
config{1}.NIRSparafac =get(handles.edit_PARAFAC,'string')
config{1}.NIRSPCA = get(handles.edit_PCA,'string');
[pathhome,file,ext]= fileparts(config{1}.NIRSPCA)
load([pathhome,'\CorrectionApply'])
[pathhome,file,ext]= fileparts(config{1}.NIRSini)

showfigure =get(handles.radio_wavfigure,'value')
trdb = str2num(get(handles.edit_tresholddb,'string'))

idsnr = 1 ;
icorr =  str2num(get(handles.edit_event,'string'));
ifile = PARCORR(icorr).file;

    load(config{1}.NIRSini)
    NC = NIRS.Cf.H.C.N;
    xname = NIRS.Dt.fir.pp(end).p{ifile};
    DATAini = fopen_NIR(xname,NC)';
    
    load(config{1}.NIRSparafac)
    NC = NIRS.Cf.H.C.N;
    xname = NIRS.Dt.fir.pp(end).p{ifile};
    DATAPARAFAC= fopen_NIR(xname,NC)';
    
    load(config{1}.NIRSPCA)
    NC = NIRS.Cf.H.C.N;
    xname = NIRS.Dt.fir.pp(end).p{ifile};
    DATAPCA= fopen_NIR(xname,NC)';
    listgood = PARCORR(icorr).listgood;
    pretime = str2num(get(handles.edit_pretime,'string'))
    posttime= str2num(get(handles.edit_posttime,'string'))
    baseline =  str2num(get(handles.edit_baseline,'string'))
    Fs = NIRS.Cf.dev.fs;% 19.5312
    freqVec = str2num(get(handles.edit_FR,'string'));
    width = str2num(get(handles.edit_morlet,'string'));
    tvec = 1/Fs:1/Fs:(1/Fs)*size(DATAini,1);
  indstart = sum(tvec<pretime);    
  indstop = sum(tvec<posttime);  
  indprebas = sum(tvec<baseline);
  num = str2num(get(handles.edit_specificch,'string'))
  if isempty(num)
      nmax = numel(listgood);
  else
      nmax = num;
  end
  

 lstch = listgood;
 HbRlist = find(lstch > NC/2);
 if ~isempty(HbRlist)
  if sum(lstch(HbRlist) - NC/2 == lstch(1:end/2))==NC/2
     msgbox('problem HbO HbR problem')
  end
  lstch(HbRlist) = [];
 end
  chHbO = lstch;
  chHbR = lstch + NC/2;

for ilistgood = 1:numel(listgood)
    ch = listgood(ilistgood);
     
     
    indt = PARCORR(icorr).indt;
    indt = indt(1)-indstart:indt(end)+indstop;
    indtnoise = PARCORR(icorr).indt- indt(1);
    indtok = indtnoise-indprebas;

  
tvecsub = tvec(indt);
Tfrini = morletcwtEd(DATAini(indt,ch),freqVec,Fs,width);
TfrPARAFAC = morletcwtEd(DATAPARAFAC(indt,ch),freqVec,Fs,width);
TfrPCA = morletcwtEd(DATAPCA(indt,ch),freqVec,Fs,width);
   PwavSNR(idsnr,1) = 10*log10(mean(mean(abs(Tfrini(indtok,:)).^2))./ mean(mean(abs(Tfrini(indtnoise,:)).^2)));    %INI
   PwavSNR(idsnr,2) = 10*log10(mean(mean(abs(Tfrini(indtok,:)).^2))./ mean(mean(abs(TfrPCA(indtnoise,:)).^2)));     %PCA
   PwavSNR(idsnr,3) = 10*log10(mean(mean(abs(Tfrini(indtok,:)).^2))./mean(mean(abs(TfrPARAFAC(indtnoise,:)).^2))); %PARAFAC
   Dini = DATAini(indt,ch);  
   DPCA = DATAPCA(indt,ch);
   DPARAFAC = DATAPARAFAC(indt,ch);
   PwavSNR(idsnr,5) = 10*log10(nanvar(Dini(indtok))/nanvar(Dini(indtnoise)));
   PwavSNR(idsnr,6) = 10*log10(nanvar(Dini(indtok))/nanvar(DPCA(indtnoise)));
   PwavSNR(idsnr,7) = 10*log10(nanvar(Dini(indtok))/nanvar(DPARAFAC(indtnoise)));
   if ilistgood<=numel(chHbO)
       ch1 = chHbO(ilistgood); 
       DiniHbO = DATAini(indt,ch1);  
       DPCAHbO = DATAPCA(indt,ch1);
       DPARAFACHbO = DATAPARAFAC(indt,ch1);
       ch2 = chHbR(ilistgood);
       DiniHbR = DATAini(indt,ch2);  
       DPCAHbR = DATAPCA(indt,ch2);
       DPARAFACHbR = DATAPARAFAC(indt,ch2);
       [PwavSNR(idsnr,9),PVAL] = corr(DiniHbO(indtok), DiniHbR(indtok));
       [PwavSNR(idsnr,10),PVAL] = corr(DiniHbO(indtnoise), DiniHbR(indtnoise));
       [PwavSNR(idsnr,11),PVAL] = corr(DPCAHbO(indtnoise), DPCAHbR(indtnoise));
       [PwavSNR(idsnr,12),PVAL] = corr(DPARAFACHbO(indtnoise), DPARAFACHbR(indtnoise));
   end
  
   
   
  if showfigure
      if  PwavSNR(idsnr,1)<trdb
    figure
    axesplot = subplot(5,1,[1,2]);hold on
    plot(tvecsub,DATAPARAFAC(indt,ch),'g','displayname','PARAFAC')
    plot(tvecsub,DATAPCA(indt,ch),'k','displayname','PCA')
    plot(tvecsub,DATAini(indt,ch),'displayname','ini')
    xlim([tvecsub(1), tvecsub(end)])
    plot([tvecsub(indtok(1)),tvecsub(indtok(end))],[0.7 0.7],'k','displayname','BASELINE SIGNAL')
    plot([tvecsub(indtnoise(1)),tvecsub(indtnoise(end))],[0.7 0.7],'r','displayname','NOISE SIGNAL')
    ylim([0.6,1.6])
    subplot(5,1,3);hold on
    S =  log(abs(Tfrini).^2); % Spectogram amplitude or power
    Measure = S;fvec = freqVec;  
    imagesc(tvecsub,freqVec,Measure(:,:,1)');
    axis xy;
    xlim([tvecsub(1), tvecsub(end)])
    title('ini')
    cmin = [min(Measure(:))];
    cmax = [max(Measure(:))]
    caxis([cmin,cmax])
  
    subplot(5,1,4)
    S =  log(abs(TfrPCA).^2); % Spectogram amplitude or power
    Measure = S;fvec = freqVec;
    imagesc(tvecsub,freqVec,Measure(:,:,1)');
    axis xy;
    caxis([cmin,cmax])
    xlim([tvecsub(1), tvecsub(end)])
    title('PCA')
      subplot(5,1,5)
    S =  log(abs(TfrPARAFAC).^2); % Spectogram amplitude or power
    Measure = S;fvec = freqVec; 
    imagesc(tvecsub,freqVec,Measure(:,:,1)');
    axis xy;
    caxis([cmin,cmax])
    xlim([tvecsub(1), tvecsub(end)])
    title('PARAFAC')
%     subplot(5,1,5)
%     Measure = abs(TfrPCA).^2 - abs(TfrPARAFAC).^2
%     imagesc(tvecsub,freqVec,Measure(:,:,1)');
%     axis xy;
%     xlim([tvecsub(1), tvecsub(end)])
%    caxis([-0.01, 0.01])
    axes(axesplot)
    title(['CH',num2str(ch),' SNR ini ', sprintf('%02.1f',(PwavSNR(idsnr,1))), ' PCA ', sprintf('%02.1f',(PwavSNR(idsnr,2))),' PARAFAC ', sprintf('%02.1f',(PwavSNR(idsnr,3))),'dB'])

    if 0
    figure   
    subplot(2,2,1)
    Measure =log(abs(Tfrini).^2)- log(abs(TfrPARAFAC).^2); 
    imagesc(tvecsub,freqVec,Measure(:,:,1)');
    cmin = [min(Measure(:))];
    cmax = [max(Measure(:))]
    caxis([cmin,cmax])
    axis xy;
    title(' ini - PARAFAC')
    subplot(2,2,3)
    plot(tvecsub,sum(Measure(:,freqint),2))
    title(['sum difference ini and PARAFAC freq',num2str(freqVec(freqint(1))),'to', num2str(freqVec(freqint(end))), 'mean ', num2str(sum(sum(Measure(:,freqint),2))) ] )

    xlim([tvecsub(1), tvecsub(end)])
    subplot(2,2,2)
    Measure =log(abs(Tfrini).^2)- log(abs(TfrPCA).^2); 
    imagesc(tvecsub,freqVec,Measure(:,:,1)');
    caxis([cmin,cmax])
    axis xy; 
    title(' ini - PCA')
    xlim([tvecsub(1), tvecsub(end)])
    subplot(2,2,4)
    plot(tvecsub,sum(Measure(:,freqint),2))
    title(['sum difference ini and PCA freq',num2str(freqVec(freqint(1))),'to', num2str(freqVec(freqint(end))), 'sum ', num2str(sum(sum(Measure(:,freqint),2)))] )
    end
      end
   end
 
idsnr= idsnr +1;
 

end
idbad = find(PwavSNR(:,1)<trdb)
% temp = 1:7
% [p,table,stats] =anova1(PwavSNR(idbad,temp))
% set(gca,'XTickLabel',{'ini','PCA','PARAFAC'},'fontsize',20)

temp = 1:3
[p,table,stats] =anova1(PwavSNR(idbad,temp))
set(gca,'XTickLabel',{'ini','PCA','PARAFAC'},'fontsize',20)
title('SNR wavelet')
h=figure
multcompare(stats)
set(gca,'yTickLabel',{'Parafac','PCA','Ini',});
set(gca,'yTick',[1,2,3]);
figure
temp = 5:7
[p,table,stats] =anova1(PwavSNR(idbad,temp))
set(gca,'XTickLabel',{'ini','PCA','PARAFAC'},'fontsize',20)
title('SNR STD')
figure
multcompare(stats)
set(gca,'yTickLabel',{'Parafac','PCA','Ini',});
set(gca,'yTick',[1,2,3]);

h=figure
temp = 9:12
[p,table,stats] =anova1(PwavSNR(idbad,temp))
set(gca,'XTickLabel',{'RAW','no corr','PCA','PARAFAC'},'fontsize',20)
title('CORR')
h=figure
multcompare(stats)
set(gca,'yTickLabel',{'Parafac','PCA','Ini','Raw'});
set(gca,'yTick',[1,2,3,4]);



[file,pathout]= uiputfile([pathhome,'SNRWAV',num2str(icorr),'.xls'])
tmp = [{'snrwavINI','PCA','PARAFAC', 'SNRvar', 'INI','PCA','PARAFAC','corr', 'raw','INI','PCA','PARAFAC'}; num2cell(PwavSNR)]
xlswrite([pathout,file],tmp);

function edit_Ini_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Ini (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Ini as text
%        str2double(get(hObject,'String')) returns contents of edit_Ini as a double


% --- Executes during object creation, after setting all properties.
function edit_Ini_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ini (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_opensubject.
function btn_opensubject_Callback(hObject, eventdata, handles)
% hObject    handle to btn_opensubject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name,path]=uigetfile('NIRS.mat')
set(handles.edit_Ini,'string',[path,name])
set(handles.edit_PCA,'string',[path, '\PCA\NIRS.mat'])
set(handles.edit_PARAFAC,'string',[path,'\PARAFAC\NIRS.mat'])

function edit_PCA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PCA as text
%        str2double(get(hObject,'String')) returns contents of edit_PCA as a double


% --- Executes during object creation, after setting all properties.
function edit_PCA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_PCA.
function btn_PCA_Callback(hObject, eventdata, handles)
% hObject    handle to btn_PCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_PARAFAC_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PARAFAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PARAFAC as text
%        str2double(get(hObject,'String')) returns contents of edit_PARAFAC as a double


% --- Executes during object creation, after setting all properties.
function edit_PARAFAC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PARAFAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_PARAFAC.
function btn_PARAFAC_Callback(hObject, eventdata, handles)
% hObject    handle to btn_PARAFAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_pretime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pretime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pretime as text
%        str2double(get(hObject,'String')) returns contents of edit_pretime as a double


% --- Executes during object creation, after setting all properties.
function edit_pretime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pretime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_event_Callback(hObject, eventdata, handles)
% hObject    handle to edit_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_event as text
%        str2double(get(hObject,'String')) returns contents of edit_event as a double


% --- Executes during object creation, after setting all properties.
function edit_event_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_posttime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_posttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_posttime as text
%        str2double(get(hObject,'String')) returns contents of edit_posttime as a double


% --- Executes during object creation, after setting all properties.
function edit_posttime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_posttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tresholddb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tresholddb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tresholddb as text
%        str2double(get(hObject,'String')) returns contents of edit_tresholddb as a double


% --- Executes during object creation, after setting all properties.
function edit_tresholddb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tresholddb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_baselinefile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_baselinefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_baselinefile as text
%        str2double(get(hObject,'String')) returns contents of edit_baselinefile as a double


% --- Executes during object creation, after setting all properties.
function edit_baselinefile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_baselinefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ICA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ICA as text
%        str2double(get(hObject,'String')) returns contents of edit_ICA as a double


% --- Executes during object creation, after setting all properties.
function edit_ICA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radio_wavfigure.
function radio_wavfigure_Callback(hObject, eventdata, handles)
% hObject    handle to radio_wavfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_wavfigure


% --- Executes on button press in btn_loadconfig.
function btn_loadconfig_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]= uigetfile('.txt');
fid = fopen([path,file]);
tline = fgetl(fid);
set(handles.edit_Ini,'string',tline);
tline = fgetl(fid);
set(handles.edit_PCA,'string',tline);
tline = fgetl(fid);
set(handles.edit_PARAFAC,'string',tline);
tline = fgetl(fid);
set(handles.edit_ICA,'string',tline);
tline = fgetl(fid);
set(handles.edit_baselinefile,'string',tline);
fclose(fid);
   


% --- Executes on button press in radio_stdFigure.
function radio_stdFigure_Callback(hObject, eventdata, handles)
% hObject    handle to radio_stdFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_stdFigure



function edit_morlet_Callback(hObject, eventdata, handles)
% hObject    handle to edit_morlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_morlet as text
%        str2double(get(hObject,'String')) returns contents of edit_morlet as a double


% --- Executes during object creation, after setting all properties.
function edit_morlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_morlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_FR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FR as text
%        str2double(get(hObject,'String')) returns contents of edit_FR as a double


% --- Executes during object creation, after setting all properties.
function edit_FR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to edit_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_baseline as text
%        str2double(get(hObject,'String')) returns contents of edit_baseline as a double


% --- Executes during object creation, after setting all properties.
function edit_baseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_specificch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_specificch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_specificch as text
%        str2double(get(hObject,'String')) returns contents of edit_specificch as a double


% --- Executes during object creation, after setting all properties.
function edit_specificch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_specificch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
