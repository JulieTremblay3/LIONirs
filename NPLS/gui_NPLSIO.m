function varargout = gui_NPLSIO(varargin)
% GUI_NPLSIO M-file for gui_NPLSIO.fig
%      GUI_NPLSIO, by itself, creates a new GUI_NPLSIO or raises the existing
%      singleton*.
%
%      H = GUI_NPLSIO returns the handle to a new GUI_NPLSIO or the handle to
%      the existing singleton*.
%
%      GUI_NPLSIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_NPLSIO.M with the given input arguments.
%
%      GUI_NPLSIO('Property','Value',...) creates a new GUI_NPLSIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_NPLSIO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_NPLSIO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_NPLSIO

% Last Modified by GUIDE v2.5 28-Nov-2022 12:30:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_NPLSIO_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_NPLSIO_OutputFcn, ...
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
 
% --- Executes just before gui_NPLSIO is made visible.
function gui_NPLSIO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_NPLSIO (see VARARGIN)

% Choose default command line output for gui_NPLSIO


%CALCUL INTERACTIF PARAFAC SELECTIONNNER POUR LA FENETRE DE TEMPS DANS LE
%GUI
handles.output = hObject;
tstart = varargin{1};
tstop = varargin{2};
NIRS = varargin{3};
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
d = PMI{currentsub}.data(cf).HRF.AvgC;

indt = tstart(end):tstop(end);
if tstart(end)==tstop(end)
    disp('Please verify time start and time stop.')
    msgbox('Please verify time start and time stop.')
end
intensnorm = d(indt,:);
  
%Detrent DATA segment for centrering 
X = 1:1:size(intensnorm,1);
Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
Mb2 =  intensnorm(1,:)'; %offset    
A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
spar = intensnorm - A;

%Ajuste le nombre de canaux au canaux selectionner pour la nirs
listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
MeasListActplotLst = zeros(size(A,2),1);
MeasListActplotLst(PMI{1}.plotLst,:)=1; 
if min(PMI{1}.plotLst-size(A,2)/2) >0
    MeasListActplotLst(PMI{1}.plotLst-size(A,2)/2,:)=1;
end
MeasListActplotLst((size(A,2)/2+1):end,1)=0; %forcer la fin a zeros


listgood = find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));



t = PMI{currentsub}.data(cf).HRF.tHRF(indt);
spar = cat(3,spar(:,1:end/2),spar(:,end/2+1:end));
spartmp = spar(:,listgood,:);
spar = spar(:,listgood,:);
% PMI{1}.tmpPARAFAC.checksumd = sum(intensnorm(:)); 
% PMI{1}.tmpPARAFAC.selected =1;
% PMI{currentsub}.tmpPARAFAC.spar = spartmp;                    %Data initial
% PMI{currentsub}.tmpPARAFAC.listgood = listgood ; 
% PMI{currentsub}.tmpPARAFAC.indt = [tstart(end),tstop(end)];%Time indice
try 
zone = PMI{currentsub}.zone;
catch
   disp('Please enter a zone it''s needed to ') 
end
%set(guiHOMER,'UserData',PMI);

nbcomp = str2num(get(handles.edit_nbComponent,'string'));
label = [];
for i=1:nbcomp
    label = [label,{['F',num2str(i)]}];
end
set(handles.popupmenu_Factor,'string',label);


%% Prepare Hemodynamic Response Function
FsN = NIRS.Cf.dev.fs; %1000/SpN; % Sampling Frequency for the NIRS data.
 Th = 20; % length in seconds of the hrf (for the current parameters it is 0 at 33 seconds)
 delay = -2.3; % additional delay assumed
Nthrf = round(Th*FsN); frametimes = (0:(Nthrf-1))./FsN;% time points for the Th seconds with the same resolution as nirs.
 H = fmridesign(frametimes,0,[1 0],[],[5.4+delay 5.2 10.8+delay 7.35 0.35 0]);
 H = squeeze(H(:,1,1));
disp(['Hemodynamic response function :', num2str(frametimes(end)), 's'])

ifile = 1;
%LOAD THE EEG
offset = NIRS.Dt.EEG.pp(end).sync_timesec{ifile};
tstartEEG = PMI{currentsub}.data(cf).HRF.tHRF(tstart)+offset;
tstopEEG = PMI{currentsub}.data(cf).HRF.tHRF(tstop)+offset;
filelist = NIRS.Dt.EEG.pp(end).p{ifile};
[dirNPLS,file,ext] = fileparts(filelist);
[dataseg,info,label,ind_dur_ch] = fopen_EEG(filelist, tstartEEG, tstopEEG );
%add_mirror_padding to avoid edge effect and cut the 
if 1 
    dataseg = [fliplr(dataseg);dataseg;fliplr(dataseg)];
end

%Ajuste le nombre de canaux au canaux selectionner pour la nirs
        try
            [nameele,x,y,z] = readelefile('C:\data\NPLS\fcNIRSepilepsie\ElectrodeDoubleBanana.ele');
            [nameele,x,y,z] = readelefile('C:\data\NPLS\fcNIRSepilepsie\Electrode1020.ele');
        catch
            disp(['Failed could not open the ele file : ', job.zonecorrelation{1} ,' for channel selection'])
            return
        end
       channeluse = [];
       infilevhdr = fullfile(dirNPLS,[file '.vhdr']);
       infovhdr = read_vhdr_brainvision(infilevhdr);
       for ilabel=1:numel(nameele)
           for ifind = 1:numel(infovhdr.label)
               if strcmp(infovhdr.label{ifind }, nameele{ilabel});
                 channeluse = [channeluse,ifind ];  
               end
           end   
       end
dataseg = dataseg(:,channeluse);
chLabels = info.name_ele(1,channeluse);
Fs = 1/(info.SamplingInterval*1e-6);   

% Downsample to 256 to perform wavelet transform 
Fsnew = 256;
dataseg1 =  resample_old( dataseg, Fsnew, Fs);
tEEGini = 1/Fs :1/Fs :size(dataseg,1)*1/Fs ;
dataseg = dataseg1; clear dataseg1
Fs = Fsnew;
freqVec = 1:1:40;
width = 5 ;
%figure;plot(dataseg(:,:))
tic; TFRseg = morletcwtEd(dataseg,freqVec,Fs,width); tori=toc; % Complex Wavelets Coefficients
%S =  abs(TFRseg);   %nonormalized
S =  abs(TFRseg)./mean(mean(abs(TFRseg),3),1).*ones(size(TFRseg));%normalized

tEEG = 1/Fs:1/Fs:size(dataseg,1)*1/Fs;
        


 %% LOAD the model HRF
 if isfield( PMI{currentsub}.data(cf),'AUX')  
    for ich=1:numel(PMI{currentsub}.data(cf).AUX.view)
        if PMI{currentsub}.data(cf).AUX.view{ich}
              fs = PMI{currentsub}.data(cf).AUX.fs{ich};
            daux = PMI{currentsub}.data(cf).AUX.data{ich};
            Label = PMI{currentsub}.data(cf).AUX.label{ich};
            taux=1/fs:1/fs:(1/fs*numel(daux));
            dModel = daux(tstart:tstop);
            %plot(taux(tstart:tstop),daux(tstart:tstop),'linewidth',2, 'displayname',Label )
        end
    end  
end
 %% Downsampling EEG Spectrogram
 Ntseg=size(S,1);Nf=size(S,2);Nc=size(S,3);
 Fsnew = FsN; %19.5312; %125; % Hz
 [~,indL] = min(mod(FsN*(1:200),1)); % look for the closest integer multiplo
 Sc = reshape(S,Ntseg,Nf*Nc);
 Sc = resample_old(Sc,round(Fsnew*indL),Fs*indL);
 %
 ntime = size(Sc,1)
 Sc(Sc(:)<=0)=eps*min(abs(Sc(:)));  % resampling can lead to non-positive values.

 Ntds = size(Sc,1);
 if 1	%
 Sc = reshape(Sc,size(Sc,1),Nf,Nc);
 tEEGdsc = ( 1/Fsnew:1/Fsnew:size(Sc,1)*1/Fsnew) - (t(end)-t(1));% divided by round(Fsnew*indL)/indL??

     1
     dNIRS = spar;
     disp(['Run parafac fNIRS on ', num2str(t(end)-t(1))  's ', 'datatime:',  num2str(t(1)), 'to', num2str(t(end)), ' CH=', num2str(size(dNIRS,2)), ' WV=2']);
     dEEG = Sc;
     disp(['Run parafac EEG on ', num2str(tEEGdsc(end)), 's Wavelet:', num2str(freqVec(1)), 'to', num2str(freqVec(end)), ' Hz Channel: ', num2str(size(Sc,3))]);

    %PARAFAC EEG WITHOUT CONVOLUTION 


    for Ncomp = 2:4
        opt = zeros(1,6); opt(3)=0; % only plotting
        const = [0 2 0]; % constraints, 1 for orthogonality, 2 for nonnegativity, 3 for unimodality
        id =  sum(tEEGdsc<=0);
        idstop = id+ size(dNIRS,1)-1;
        [Factors,it,err,corcondia] = parafac(Sc,Ncomp,opt,const);       
        [AEEG,BEEG,CEEG]=fac2let(Factors);
        FactorComp{Ncomp}.EEG =  Factors;  
        HRF_TimeComponent = conv2(H.X(:,1,1)/10,AEEG);       
        %figure;plot(t,HRF_TimeComponent(id:idstop,:))
        [rEEGHRF_MODEL(Ncomp,1:Ncomp),pvalEEGHRF_MODEL(Ncomp,1:Ncomp)]  =  corr(HRF_TimeComponent(id:idstop,:),dModel );
        rEEGHRF_MODEL(Ncomp,6)  = corcondia;
%         if concordia < 0.80
%             rEEGHRF_MODEL(Ncomp,1:4)=0
%         end
       
        const = [0 0 0]; % constraints, 1 for orthogonality, 2 for nonnegativity, 3 for unimodality
        [Factors,it,err,corcondia] = parafac(dNIRS,Ncomp,opt,const);
        FactorComp{Ncomp}.NIRS =  Factors;
        [ANIRS,BNIRS,CNIRS]=fac2let(Factors);
        [rdNIRS_MODEL(Ncomp,1:Ncomp),pvalNIRS_MODEL(Ncomp,1:Ncomp)] = corr(ANIRS,dModel )        
        rdNIRS_MODEL(Ncomp,6) = corcondia;
    end

    %chose highest CORR for fNIRS and model 
    figure ;
    [val, NcompNIRSmax]=max(abs(rdNIRS_MODEL(:,1:5)));
    [val, icompNIRSmax]= max(val);
    NcompNIRSmax =  NcompNIRSmax( icompNIRSmax);
     [Anirs,Bnirs,Cnirs]=fac2let(FactorComp{  NcompNIRSmax}.NIRS);
    %FLIP Time and channel if corr neg
    if rdNIRS_MODEL(NcompNIRSmax,icompNIRSmax)<0
        Anirs = -Anirs;
        Bnirs = -Bnirs;
        rdNIRS_MODEL = -rdNIRS_MODEL;
    end
    Bchannel  = nan(numel(NIRS.Cf.H.C.ok)/2,size(Bnirs,2));
    Bchannel(listgood,:)=Bnirs;
    chplot = [];
    idtick  = [];
    try
    for izone=1:numel(zone.plotLst)
        chplot = [chplot, zone.plotLst{izone}];
        idtick = [idtick, numel(chplot)];
    end
    catch
       chplot = 1:size(Bchannel,1) ;
    end
    subplot(5,3,1); plot(t,Anirs(:,icompNIRSmax));xlabel('Time (s)');ylabel('fNIRS'); %ts_eeg(indseg)
    hold on;plot(t, dModel);
    title('BEST CORRELATION NIRS Model')
    subplot(5,3,3); plot(Bchannel(chplot,icompNIRSmax)); axis tight; xlabel('Channels');ylabel('Normalized');
   try; set(gca,'Xtick',idtick ,'XTickLabel',zone.label);catch;end;
    title(['CoorNIRS model=', num2str(rdNIRS_MODEL(  NcompNIRSmax,icompNIRSmax)),'p=', num2str(pvalNIRS_MODEL(  NcompNIRSmax,icompNIRSmax))]);
    subplot(5,3,2); plot(Cnirs(:,icompNIRSmax)); axis tight; xlabel('Wavelenght');ylabel('Normalized');
   [~,fileNIRS,~] =fileparts(NIRS.Dt.fir.pp(1).p{1});
   title(fileNIRS);
    %Save component  for topo or stat 
    A = nan(numel(NIRS.Cf.H.C.ok)/2,1);
    A(listgood)=Bnirs(:,icompNIRSmax);

      fileout = fullfile(dirNPLS,['CompNIRS_Time',...
        sprintf('%2.1f',PMI{currentsub}.data(cf).HRF.tHRF(tstart)),...
        'to',sprintf('%2.1f',PMI{currentsub}.data(cf).HRF.tHRF(tstop)),...
        '_corrHRF', sprintf('%0.2f',rdNIRS_MODEL(  NcompNIRSmax,icompNIRSmax)),'.mat']);
    save(fileout,'A' );
    disp(['Save: ',fileout ]);
    

    
    %[Xfactors,Yfactors,Core,B,ypred,ssx,ssy] = npls(Sc(id:idstop,:,:),dNIRS(:,:,:),Ncomp);
    %CORR 
    HRFNIRS_Corrwith_MODEL = Anirs(:,icompNIRSmax);
    for Ncomp  = 2:4
          [AEEG,BEEG,CEEG]=fac2let( FactorComp{Ncomp}.EEG);
          HRF_TimeComponent = conv2(H.X(:,1,1)/10,AEEG);  
         [rEEGHRF_NIRS(Ncomp,1:Ncomp),pvalEEGHRF_NIRS(Ncomp,1:Ncomp)]  =  corr(HRF_TimeComponent(id:idstop,:),HRFNIRS_Corrwith_MODEL );
         %rEEGHRF_Model(Ncomp,1:Ncomp)  =  corr(HRF_TimeComponent(id:idstop,:),HRFNIRS_Corrwith_MODEL ); rEEGHRF_MODEL
    end
     [val, NcompmaxEEG] = max(abs(rEEGHRF_NIRS));
      [val, icompmaxEEG]= max(val);
      NcompmaxEEG =       NcompmaxEEG ( icompmaxEEG);
     [AEEG,BEEG,CEEG]=fac2let(FactorComp{NcompmaxEEG}.EEG);
    subplot(5,3,4); plot(t,detrend(AEEG(id:idstop,icompmaxEEG)));xlabel('Time (s)');ylabel('EEG'); %ts_eeg(indseg)
    hold on;plot(t,dModel)
    title('BEST CORRELATION NIRS EEG')
    HRF_TimeComponent = conv2(H.X(:,1,1)/10,AEEG);
     hold on;plot(t, detrend(HRF_TimeComponent(id:idstop,icompmaxEEG)))
    subplot(5,3,5); plot(freqVec,BEEG(:,icompmaxEEG)); axis tight; xlabel('Frequence');ylabel('Normalized');
     title(file);
    subplot(5,3,6); plot(CEEG(:,icompmaxEEG)); axis tight; ylabel('Normalized');
    set(gca,'Xtick',1:numel(chLabels),'XTickLabel',chLabels); xlabel('Channels');ylabel('Normalized');

    title(['EEG Coor NIRS =', num2str(rEEGHRF_NIRS(  NcompmaxEEG,icompmaxEEG)),' p=',num2str(pvalEEGHRF_NIRS(  NcompmaxEEG,icompmaxEEG)),'EEG Coor model =', num2str(rEEGHRF_MODEL(  NcompmaxEEG,icompmaxEEG)), 'p=',num2str(pvalEEGHRF_MODEL(  NcompmaxEEG,icompmaxEEG))]);
    tmp=abs(rEEGHRF_NIRS(NcompmaxEEG,:));
    tmp(icompmaxEEG) = 0;
     [val, icompmaxEEG] = max(tmp);
     subplot(5,3,7); plot(t,detrend(AEEG(id:idstop,icompmaxEEG)));xlabel('Time (s)');ylabel('EEG'); %ts_eeg(indseg)
    hold on;plot(t,dModel)
    title('BEST CORRELATION NIRS EEG')
    HRF_TimeComponent = conv2(H.X(:,1,1)/10,AEEG);
     hold on;plot( t,detrend(HRF_TimeComponent(id:idstop,icompmaxEEG)),'displayname','EEG with HRF convolution')
    subplot(5,3,8); plot(freqVec,BEEG(:,icompmaxEEG)); axis tight; xlabel('Frequence');ylabel('Normalized');
    subplot(5,3,9); plot(CEEG(:,icompmaxEEG)); axis tight; ylabel('Normalized');
    set(gca,'Xtick',1:numel(chLabels),'XTickLabel',chLabels); 
    title(['EEG Coor NIRS =', num2str(rEEGHRF_NIRS(  NcompmaxEEG,icompmaxEEG)),' ','EEG Coor model =', num2str(rEEGHRF_MODEL(  NcompmaxEEG,icompmaxEEG))]);

    
    valtest = rEEGHRF_MODEL(NcompmaxEEG,:);
    
    subplot(5,3,13); plot(detrend(AEEG(id:idstop,:)));xlabel('Time (s)');ylabel('EEG'); %ts_eeg(indseg)
    %hold on;plot(dModel)
    title('BEST CORRELATION NIRS EEG')
    HRF_TimeComponent = conv2(H.X(:,1,1)/10,AEEG);
    hold on;plot( detrend(HRF_TimeComponent(id:idstop,:)),'displayname','EEG with HRF convolution');
   
    subplot(5,3,14);hold on 
     for i=2:numel(NcompmaxEEG)
        plot(freqVec,BEEG(:,i),'displayname', ['Corr:', valtest(i)] );
     end
     axis tight; xlabel('Frequence');ylabel('Normalized');
    subplot(5,3,15); plot(CEEG(:,:)); axis tight; ylabel('Normalized');

    FactorNIRS = FactorComp{NcompNIRSmax}.NIRS;
    FactorEEG = FactorComp{NcompmaxEEG}.EEG;
    
    rEEGHRF_NIRS
    rEEGHRF_MODEL
    
    disp('Save Factor' )
    
    disp('done: gui_NPLSIO' )
  %  
 elseif 1 %old test 
    %% Convolution with the Hemodynamic Response Function  whole spectrum
    Sc = conv2(H.X(:,1,1),Sc);
    Sc = Sc(1:ntime,:);


  
  Fsnew = 256;
  FsN
Sc1 =  resample_old( Sc,   FsN*10000,Fsnew*10000);

Sc2 = reshape(Sc1,size(Sc1,1),Nf,Nc);
  
tconv = 1/Fsnew:1/Fsnew:size(Sc1,1)*1/Fsnew;
 id_cut_conv = numel(tconv) - numel(H.X(:,1,1));
Ntds = size(Sc1,1);
tdsc = tEEGdsc(1)+(0:Ntds-1)./Fsnew; % divided by round(Fsnew*indL)/indL??
id_cut_start = 1; %round(10/(1/Fsnew))
%PARAFAC EEG WITHOUT CONVOLUTION 
opt = zeros(1,6); opt(3)=0;  figure;% only plotting
const = [0 0 0]; % constraints, 1 for orthogonality, 2 for nonnegativity, 3 for unimodality
%figure; % This is for plotting concordia
Ncomp = 3
%Sc EEG sans conv
%Sc1 EEG avec conv
dNIRS = spar;
%[Xfactors,Yfactors,Core,B,ypred,ssx,ssy] = npls(dNIRS(1:id_cut_conv,:,:),dModel(1:id_cut_conv),Ncomp);


[Xfactors,Yfactors,Core,B,ypred,ssx,ssy] = npls(Sc(:,:,:),dNIRS(:,:,:),Ncomp);
%[Xfactors,Yfactors,Core,B,ypred,ssx,ssy] = npls(HRF_TimeComponent(1:id_cut_conv,:),dNIRS,Ncomp);
%imagesc(abs(corrcoef([Anirs(:,:),HRF_TimeComponent(1:id_cut_conv,:)]))) 
tHRF = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:id_cut_conv*1/NIRS.Cf.dev.fs;
[AEEG,BEEG,CEEG]=fac2let(Xfactors);
figure; % This is for plotting concordia
subplot(5,3,1); plot(tHRF(id_cut_start:id_cut_conv),AEEG);xlabel('Time (s)');ylabel('EEG spectrum with HRF conv'); %ts_eeg(indseg)
subplot(5,3,2); plot(freqVec, BEEG); axis tight; xlabel('Frequence (Hz)');ylabel('Normalized');
subplot(5,3,3); plot(CEEG); 
[Anirs,Bnirs,Cnirs]=fac2let(Yfactors);
subplot(5,3,4); plot(tHRF(id_cut_start:id_cut_conv),Anirs);xlabel('Time (s)');ylabel('fNIRS'); %ts_eeg(indseg)
subplot(5,3,5); plot( Bnirs); axis tight; xlabel('Channels');ylabel('Normalized');
subplot(5,3,6); plot(Cnirs); axis tight; xlabel('Wavelenght');ylabel('Normalized');
try
for icomp = 1:Ncomp
   rEEG(icomp) = corr(dModel(1:size(AEEG,1),:),AEEG(:,icomp)) 
   rNIRS(icomp) =  corr(dModel(1:size(AEEG,1),:),Anirs(:,icomp)) 
end
catch    
end


dNIRS = spar;
%PARAFAC EEG WITHOUT CONVOLUTION 
opt = zeros(1,6); opt(3)=0; % only plotting
const = [0 2 0]; % constraints, 1 for orthogonality, 2 for nonnegativity, 3 for unimodality
%figure; % This is for plotting concordia
Ncomp = 2
[FactorsEEG,it,err,corcondia] = parafac(Sc(id_cut_start:id_cut_conv,:,:),Ncomp,opt,const);
const = [0 0 0]; % constraints, 1 for orthogonality, 2 for nonnegativity, 3 for unimodality

[FactorsNIRS,it,err,corcondia] = parafac(dNIRS(id_cut_start:id_cut_conv,:,:),Ncomp,opt,const);
%[Xfactors,Yfactors,Core,B,ypred,ssx,ssy] = npls(Sc(1:id_cut_conv,:,:),dNIRS(1:id_cut_conv,:,:),Ncomp);
%[Xfactors,Yfactors,Core,B,ypred,ssx,ssy] = npls(HRF_TimeComponent(1:id_cut_conv,:),dNIRS,Ncomp);
%imagesc(abs(corrcoef([Anirs(:,:),HRF_TimeComponent(1:id_cut_conv,:)]))) 
tHRF = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:id_cut_conv*1/NIRS.Cf.dev.fs;
[AEEG,BEEG,CEEG]=fac2let(FactorsEEG);
subplot(5,3,7); plot(tHRF(id_cut_start:id_cut_conv),AEEG);xlabel('Time (s)');ylabel('EEG spectrum without conv'); %ts_eeg(indseg)
subplot(5,3,8); plot(freqVec, BEEG); axis tight; xlabel('Frequence (Hz)');ylabel('Normalized');
subplot(5,3,9); plot(CEEG);  set(gca,'Xtick',1:numel(chLabels),'XTickLabel',chLabels); xlabel('Channels');ylabel('Normalized');

HRF_TimeComponent = conv2(H.X(:,1,1)/10,AEEG);
subplot(5,3,10); plot(tHRF(id_cut_start:id_cut_conv),HRF_TimeComponent(id_cut_start:id_cut_conv,:,:));xlabel('Time (s)');ylabel('EEG spectrum without conv'); %ts_eeg(indseg)
try
subplot(5,3,11); plot(taux(tstart:tstop),daux(tstart:tstop));xlabel('Time (s)');ylabel('EEG spectrum without conv'); %ts_eeg(indseg)
hold on

catch
end
[Anirs,Bnirs,Cnirs]=fac2let(FactorsNIRS);
subplot(5,3,13); plot(tHRF(id_cut_start:id_cut_conv),Anirs);xlabel('Time (s)');ylabel('fNIRS'); %ts_eeg(indseg)
subplot(5,3,14); plot( Bnirs); axis tight; xlabel('Channels');ylabel('Normalized');
subplot(5,3,15); plot(Cnirs); axis tight; xlabel('Wavelenght');ylabel('Normalized');

%EEG x Model
try
EEG_Mod = corr(HRF_TimeComponent(id_cut_start:id_cut_conv,:),daux(id_cut_start:id_cut_conv,:))
NIRS_Mod =corr(Anirs(id_cut_start:id_cut_conv,:),daux(id_cut_start:id_cut_conv,:))
catch
end
NIRS_EEG = corr(Anirs(id_cut_start:id_cut_conv,:),HRF_TimeComponent(id_cut_start:id_cut_conv,:))
figure;

subplot(2,3,1);hold on
[val,idnirs] = max(NIRS_Mod)
[val,ideeg] =  max(EEG_Mod)
NIRS_Mod(idnirs) 
EEG_Mod(ideeg) 
subplot(2,3,1);
plot(daux(id_cut_start:id_cut_conv,:),'displayname','modelHRF')
plot(HRF_TimeComponent(id_cut_start:id_cut_conv,ideeg),'displayname','convEEG')
title(['EEG vs model =', num2str(corr(HRF_TimeComponent(id_cut_start:id_cut_conv,ideeg),daux(id_cut_start:id_cut_conv,:))) ])
subplot(2,3,2);hold on
plot(Anirs(id_cut_start:id_cut_conv,idnirs),'displayname','NIRS')
 plot(daux(id_cut_start:id_cut_conv,:),'displayname','modelHRF')
 title(['NIRS vs model =', num2str(corr(Anirs(id_cut_start:id_cut_conv,idnirs),daux(id_cut_start:id_cut_conv,:))) ])
subplot(2,3,3);hold on
plot(Anirs(id_cut_start:id_cut_conv,idnirs),'displayname','NIRS')
plot(HRF_TimeComponent(id_cut_start:id_cut_conv,ideeg),'displayname','convEEG')
 title(['NIRS vs EEG =', num2str(corr(Anirs(id_cut_start:id_cut_conv,idnirs),HRF_TimeComponent(id_cut_start:id_cut_conv,ideeg))) ])

 
%NIRS x Model
subplot(2,3,4); plot(tHRF(id_cut_start:id_cut_conv),AEEG(:,ideeg));xlabel('Time (s)');ylabel('EEG spectrum without conv'); %ts_eeg(indseg)
subplot(2,3,5); plot(freqVec, BEEG(:,ideeg)); axis tight; xlabel('Frequence (Hz)');ylabel('Normalized');
subplot(2,3,6); plot(CEEG(:,ideeg));  set(gca,'Xtick',1:numel(chLabels),'XTickLabel',chLabels); xlabel('Channels');ylabel('Normalized');

end 
%EEG x NIRS

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_NPLSIO wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_NPLSIO_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in popup_parafactID.
function popup_parafactID_Callback(hObject, eventdata, handles)
% hObject    handle to popup_parafactID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_parafactID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_parafactID


% --- Executes during object creation, after setting all properties.
function popup_parafactID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_parafactID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Parafac.
function btn_Parafac_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Parafac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
spar = PMI{currentsub}.tmpPARAFAC.spar;
Nc = str2num(get(handles.edit_nbComponent,'string'));
opt(1) = 1e-6; 
opt(2) = 0;    %initialisation 10 all methods fit using DTLD/GRAM
opt(3) = 0;     %plotting
opt(5) = 0;     %how often to show fit.
opt(6) = 10;     %Max num iterations
const = [0 0 0]; % constraints 0- nothing; 1-orthog; 2- nonneg; 3- unimod
Oldload{2} = []; % rand(21,3);

if get(handles.radio_equallambda,'value'); 
    fixMode = [0 0 1];
   % weights{3} = [1,0,,0,0]
   weights{1} = [];
   weights{2} = [];
   weights{3} = [1 1;
       1,0]
else
    fixMode = [0 0 0];
    weights = []; %Mean for example put one at good time point and zero at noisy one
end

% figure 
% plot(spar(:,:,1))
tic
[Factors,it,err,corcondia] = parafac(spar(:,:,:),Nc,opt,const,Oldload,fixMode,weights);
toc
PMI{currentsub}.tmpPARAFAC.Factors = Factors;

set(guiHOMER,'UserData',PMI);
guidata(hObject, handles);

set(handles.text_ConcordiaVal,'string', [num2str(round(corcondia)),' ',num2str(round(err))] )
set(handles.popupmenu_Factor,'value',1)
resetview(handles,0)

function edit_nbComponent_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nbComponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nbComponent as text
%        str2double(get(hObject,'String')) returns contents of edit_nbComponent as a double
nbcomp = str2num(get(handles.edit_nbComponent,'string'));
label = [];
for i=1:nbcomp
    label = [label,{['F',num2str(i)]}];
end

set(handles.popupmenu_Factor,'string',label)
 %for the topo
btn_Parafac_Callback(hObject, eventdata, handles)


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

function resetview(handles,newfigure)
%All normal 

axes(handles.axes_timecourse);
cla
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');


    FacTime = PMI{currentsub}.tmpPARAFAC.Factors{1};
    FacSpatial = PMI{currentsub}.tmpPARAFAC.Factors{2};
    FacHemo = PMI{currentsub}.tmpPARAFAC.Factors{3};
    Nc = str2num(get(handles.edit_nbComponent,'string'));
    colorlist = [0 0 1,
        1 0 0,
        0 1 0,
        0 0 0,
        1 0.5 0];
    cf = PMI{currentsub}.currentFile;
    if newfigure==0
        axes(handles.axes_timecourse);
        cla
        hold on
    else
        figure
        subplot(2,2,1)
        hold on
    end
    t = PMI{currentsub}.data(cf).HRF.tHRF(PMI{currentsub}.tmpPARAFAC.indt(1):PMI{currentsub}.tmpPARAFAC.indt(2));

    for i=1:Nc
         if PMI{1}.tmpPARAFAC.selected==i
            plot(t, FacTime(:,i),'color',colorlist(i,:),'linewidth',4 )
         else
            plot(t, FacTime(:,i),'color',colorlist(i,:) )
         end
    end

    if newfigure==0
        axes(handles.axes_Topo);
        cla
        hold on
    else
        subplot(2,2,2)
        hold on
    end
    for i=1:Nc
        if PMI{1}.tmpPARAFAC.selected==i
            plot(FacSpatial(:,i),'color',colorlist(i,:),'linewidth',4 )
        else
            plot(FacSpatial(:,i),'color',colorlist(i,:) )
        end
    end
    if newfigure==0
        axes(handles.axes_hemo);
        cla
        hold on
     else
        subplot(2,2,3)
        hold on
    end
    for i=1:Nc
        if PMI{1}.tmpPARAFAC.selected==i
            plot(FacHemo(:,i),'color',colorlist(i,:),'linewidth',4 )
        else
            plot(FacHemo(:,i),'color',colorlist(i,:) )        
        end
    end
    


% --- Executes on selection change in popupmenu_Factor.
function popupmenu_Factor_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Factor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Factor
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
PMI{1}.tmpPARAFAC.selected =get(handles.popupmenu_Factor,'value');
set(guiHOMER,'UserData',PMI);
resetview(handles,0)
% guiHelmet = getappdata(0,'guiHelmet');
% if ~isempty(guiHelmet)  
%         Helmethandles = guihandles(guiHelmet);
%         fhresetview = getappdata(guiHelmet,'fhresetview');
%         fhresetview(Helmethandles);
% end
% --- Executes during object creation, after setting all properties.
function popupmenu_Factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function CreateNewFigure_Callback(hObject, eventdata, handles)
% hObject    handle to CreateNewFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_3dhelmet.
function btn_3dhelmet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_3dhelmet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

option.conc = 1;
option.SPMgui = 1;
IO_HelmetMTG_Display(handles,option)


% --- Executes on button press in radio_equallambda.
function radio_equallambda_Callback(hObject, eventdata, handles)
% hObject    handle to radio_equallambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_equallambda


% --------------------------------------------------------------------
function context_Newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to context_Newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetview(handles,1)


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function  [y, h] = resample_old( x, p, q, N, bta )
%RESAMPLE  Change the sampling rate of a signal.
%   Y = RESAMPLE(X,P,Q) resamples the sequence in vector X at P/Q times
%   the original sample rate using a polyphase implementation.  Y is P/Q 
%   times the length of X (or the ceiling of this if P/Q is not an integer).  
%   P and Q must be positive integers.
%
%   RESAMPLE applies an anti-aliasing (lowpass) FIR filter to X during the 
%   resampling process, and compensates for the filter's delay.  The filter 
%   is designed using FIRLS.  RESAMPLE provides an easy-to-use alternative
%   to UPFIRDN, relieving the user of the need to supply a filter or
%   compensate for the signal delay introduced by filtering.
%
%   In its filtering process, RESAMPLE assumes the samples at times before
%   and after the given samples in X are equal to zero. Thus large
%   deviations from zero at the end points of the sequence X can cause
%   inaccuracies in Y at its end points.
%
%   Y = RESAMPLE(X,P,Q,N) uses a weighted sum of 2*N*max(1,Q/P) samples of X 
%   to compute each sample of Y.  The length of the FIR filter RESAMPLE applies
%   is proportional to N; by increasing N you will get better accuracy at the 
%   expense of a longer computation time.  If you don't specify N, RESAMPLE uses
%   N = 10 by default.  If you let N = 0, RESAMPLE performs a nearest
%   neighbor interpolation; that is, the output Y(n) is X(round((n-1)*Q/P)+1)
%   ( Y(n) = 0 if round((n-1)*Q/P)+1 > length(X) ).
%
%   Y = RESAMPLE(X,P,Q,N,BTA) uses BTA as the BETA design parameter for the 
%   Kaiser window used to design the filter.  RESAMPLE uses BTA = 5 if
%   you don't specify a value.
%
%   Y = RESAMPLE(X,P,Q,B) uses B to filter X (after upsampling) if B is a 
%   vector of filter coefficients.  RESAMPLE assumes B has odd length and
%   linear phase when compensating for the filter's delay; for even length 
%   filters, the delay is overcompensated by 1/2 sample.  For non-linear 
%   phase filters consider using UPFIRDN.
%
%   [Y,B] = RESAMPLE(X,P,Q,...) returns in B the coefficients of the filter
%   applied to X during the resampling process (after upsampling).
%
%   If X is a matrix, RESAMPLE resamples the columns of X.
%
%   See also UPFIRDN, INTERP, DECIMATE, FIRLS, KAISER, INTFILT,
%   MFILT/FIRSRC in the Filter Design Toolbox.

%   NOTE-1: digital anti-alias filter is desiged via windowing

%   Author(s): James McClellan, 6-11-93
%              Modified to use upfirdn, T. Krauss, 2-27-96
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.9.4.7 $  $Date: 2007/12/14 15:05:59 $

if nargin < 5,  bta = 5;  end   %--- design parameter for Kaiser window LPF
if nargin < 4,   N = 10;   end
if abs(round(p))~=p || p==0, error(generatemsgid('MustBePosInteger'),'P must be a positive integer.'), end
if abs(round(q))~=q || q==0, error(generatemsgid('MustBePosInteger'),'Q must be a positive integer.'), end

[p,q] = rat( p/q, 1e-12 );  %--- reduce to lowest terms 
   % (usually exact, sometimes not; loses at most 1 second every 10^12 seconds)
if (p==1) && (q==1)
    y = x; 
    h = 1;
    return
end
pqmax = max(p,q);
if length(N)>1      % use input filter
   L = length(N);
   h = N;
else                % design filter
   if( N>0 )
      fc = 1/2/pqmax;
      L = 2*N*pqmax + 1;
      h = p*firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,bta)' ;
      % h = p*fir1( L-1, 2*fc, kaiser(L,bta)) ;
   else
      L = p;
      h = ones(1,p);
   end
end

Lhalf = (L-1)/2;
isvect = any(size(x)==1);
if isvect
    Lx = length(x);
else
    Lx = size(x, 1);
end

% Need to delay output so that downsampling by q hits center tap of filter.
nz = floor(q-mod(Lhalf,q));
z = zeros(1,nz);
h = [z h(:).'];  % ensure that h is a row vector.
Lhalf = Lhalf + nz;

% Number of samples removed from beginning of output sequence 
% to compensate for delay of linear phase filter:
delay = floor(ceil(Lhalf)/q);

% Need to zero-pad so output length is exactly ceil(Lx*p/q).
nz1 = 0;
while ceil( ((Lx-1)*p+length(h)+nz1 )/q ) - delay < ceil(Lx*p/q)
    nz1 = nz1+1;
end
h = [h zeros(1,nz1)];

% ----  HERE'S THE CALL TO UPFIRDN  ----------------------------
y = upfirdn(x,h,p,q);

% Get rid of trailing and leading data so input and output signals line up
% temporally:
Ly = ceil(Lx*p/q);  % output length
% Ly = floor((Lx-1)*p/q+1);  <-- alternately, to prevent "running-off" the
%                                data (extrapolation)
if isvect
    y(1:delay) = [];
    y(Ly+1:end) = [];
else
    y(1:delay,:) = [];
    y(Ly+1:end,:) = [];
end

h([1:nz (end-nz1+1):end]) = [];  % get rid of leading and trailing zeros 
                                 % in case filter is output

function [h,a]=firls(N,F,M,W,ftype)
% FIRLS Linear-phase FIR filter design using least-squares error minimization.
%   B=FIRLS(N,F,A) returns a length N+1 linear phase (real, symmetric
%   coefficients) FIR filter which has the best approximation to the
%   desired frequency response described by F and A in the least squares
%   sense. F is a vector of frequency band edges in pairs, in ascending
%   order between 0 and 1. 1 corresponds to the Nyquist frequency or half
%   the sampling frequency. A is a real vector the same size as F
%   which specifies the desired amplitude of the frequency response of the
%   resultant filter B. The desired response is the line connecting the
%   points (F(k),A(k)) and (F(k+1),A(k+1)) for odd k; FIRLS treats the
%   bands between F(k+1) and F(k+2) for odd k as "transition bands" or
%   "don't care" regions. Thus the desired amplitude is piecewise linear
%   with transition bands.  The integrated squared error is minimized.
%
%   For filters with a gain other than zero at Fs/2, e.g., highpass
%   and bandstop filters, N must be even.  Otherwise, N will be
%   incremented by one. Alternatively, you can use a trailing 'h' flag to
%   design a type 4 linear phase filter and avoid incrementing N.
%
%   B=FIRLS(N,F,A,W) uses the weights in W to weight the error. W has one
%   entry per band (so it is half the length of F and A) which tells
%   FIRLS how much emphasis to put on minimizing the integral squared error
%   in each band relative to the other bands.
%
%   B=FIRLS(N,F,A,'Hilbert') and B=FIRLS(N,F,A,W,'Hilbert') design filters
%   that have odd symmetry, that is, B(k) = -B(N+2-k) for k = 1, ..., N+1.
%   A special case is a Hilbert transformer which has an approx. amplitude
%   of 1 across the entire band, e.g. B=FIRLS(30,[.1 .9],[1 1],'Hilbert').
%
%   B=FIRLS(N,F,A,'differentiator') and B=FIRLS(N,F,A,W,'differentiator')
%   also design filters with odd symmetry, but with a special weighting
%   scheme for non-zero amplitude bands. The weight is assumed to be equal
%   to the inverse of frequency, squared, times the weight W. Thus the
%   filter has a much better fit at low frequency than at high frequency.
%   This designs FIR differentiators.
%
%   % Example of a length 31 lowpass filter.
%   h=firls(30,[0 .1 .2 .5]*2,[1 1 0 0]);
%   fvtool(h);
%
%   % Example of a length 45 lowpass differentiator.
%   h=firls(44,[0 .3 .4 1],[0 .2 0 0],'differentiator');
%   fvtool(h);
%
%   % Example of a length 26 type 4 highpass filter.
%   h=firls(25,[0 .4 .5 1],[0 0 1 1],'h');
%   fvtool(h);
%
%   See also FIRPM, FIR1, FIR2, FREQZ and FILTER.

%       Author(s): T. Krauss
%   History: 10-18-91, original version
%            3-30-93, updated
%            9-1-95, optimize adjacent band case
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.7 $  $Date: 2009/05/23 08:13:29 $

% check number of arguments, set up defaults.
error(nargchk(3,5,nargin,'struct'));

if (max(F)>1) || (min(F)<0)
    error(generatemsgid('InvalidRange'),'Frequencies in F must be in range [0,1].')
end
if (rem(length(F),2)~=0)
    error(generatemsgid('InvalidDimensions'),'F must have even length.');
end
if (length(F) ~= length(M))
    error(generatemsgid('InvalidDimensions'),'F and A must be equal lengths.');
end
if (nargin==3),
    W = ones(length(F)/2,1);
    ftype = '';
end
if (nargin==4),
    if ischar(W),
        ftype = W; W = ones(length(F)/2,1);
    else
        ftype = '';
    end
end
if (nargin==5),
    if isempty(W),
        W = ones(length(F)/2,1);
    end
end
if isempty(ftype)
    ftype = 0;  differ = 0;
else
    ftype = lower(ftype);
    if strcmpi(ftype,'h') || strcmpi(ftype,'hilbert')
        ftype = 1;  differ = 0;
    elseif strcmpi(ftype,'d') || strcmpi(ftype,'differentiator')
        ftype = 1;  differ = 1;
    else
        error(generatemsgid('InvalidEnum'),'Requires symmetry to be ''Hilbert'' or ''differentiator''.')
    end
end

% Check for valid filter length
[N,msg1,msg2] = firchk(N,F(end),M,ftype);
if ~isempty(msg1)
    error(generatemsgid('InvalidFilterOrder'),msg1);
end

if ~isempty(msg2),
    msg2 = sprintf([msg2,'\r',...
        '\nAlternatively, you can pass a trailing ''h'' argument,\r',...
        'as in firls(N,F,A,W,''h''), to design a type 4 linear phase filter.']);
    warning(generatemsgid('OrderIncreasedByOne'),msg2); 
end



N = N+1;                   % filter length
F=F(:)/2;  M=M(:);  W=sqrt(W(:));  % make these guys columns
dF = diff(F);

if (length(F) ~= length(W)*2)
    error(generatemsgid('InvalidDimensions'),'There should be one weight per band.');
end;
if any(dF<0),
    error(generatemsgid('InvalidFreqVec'),'Frequencies in F must be nondecreasing.')
end

% Fix for 67187
if all(dF(2:2:length(dF)-1)==0) && length(dF) > 1,
    fullband = 1;
else
    fullband = 0;
end
if all((W-W(1))==0)
    constant_weights = 1;
else
    constant_weights = 0;
end

L=(N-1)/2;

Nodd = rem(N,2);

if (ftype == 0),  % Type I and Type II linear phase FIR
    % basis vectors are cos(2*pi*m*f) (see m below)
    if ~Nodd
        m=(0:L)+.5;   % type II
    else
        m=(0:L);      % type I
    end
    k=m';
    need_matrix = (~fullband) || (~constant_weights);
    if need_matrix
        I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
        I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
        G=zeros(size(I1));
    end

    if Nodd
        k=k(2:length(k));
        b0=0;       %  first entry must be handled separately (where k(1)=0)
    end;
    b=zeros(size(k));
    for s=1:2:length(F),
        m=(M(s+1)-M(s))/(F(s+1)-F(s));    %  slope
        b1=M(s)-m*F(s);                   %  y-intercept
        if Nodd
            b0 = b0 + (b1*(F(s+1)-F(s)) + m/2*(F(s+1)*F(s+1)-F(s)*F(s)))...
                * abs(W((s+1)/2)^2) ;
        end
        b = b+(m/(4*pi*pi)*(cos(2*pi*k*F(s+1))-cos(2*pi*k*F(s)))./(k.*k))...
            * abs(W((s+1)/2)^2);
        b = b + (F(s+1)*(m*F(s+1)+b1)*sinc(2*k*F(s+1)) ...
            - F(s)*(m*F(s)+b1)*sinc(2*k*F(s))) ...
            * abs(W((s+1)/2)^2);
        if need_matrix
            G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))+sinc(2*I2*F(s+1))) ...
                - .5*F(s)*(sinc(2*I1*F(s))+sinc(2*I2*F(s))) ) ...
                * abs(W((s+1)/2)^2);
        end
    end;
    if Nodd
        b=[b0; b];
    end;

    if need_matrix
        a=G\b;
    else
        a=(W(1)^2)*4*b;
        if Nodd
            a(1) = a(1)/2;
        end
    end
    if Nodd
        h=[a(L+1:-1:2)/2; a(1); a(2:L+1)/2].';
    else
        h=.5*[flipud(a); a].';
    end;
elseif (ftype == 1),  % Type III and Type IV linear phase FIR
    %  basis vectors are sin(2*pi*m*f) (see m below)
    if (differ),      % weight non-zero bands with 1/f^2
        do_weight = ( abs(M(1:2:length(M))) +  abs(M(2:2:length(M))) ) > 0;
    else
        do_weight = zeros(size(F));
    end

    if Nodd
        m=(1:L);      % type III
    else
        m=(0:L)+.5;   % type IV
    end;
    k=m';
    b=zeros(size(k));

    need_matrix = (~fullband) || (any(do_weight)) || (~constant_weights);
    if need_matrix
        I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
        I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
        G=zeros(size(I1));
    end

    i = sqrt(-1);
    for s=1:2:length(F),
        if (do_weight((s+1)/2)),      % weight bands with 1/f^2
            if F(s) == 0, F(s) = 1e-5; end     % avoid singularities
            m=(M(s+1)-M(s))/(F(s+1)-F(s));
            b1=M(s)-m*F(s);
            snint1 = sineint(2*pi*k*F(s+1)) - sineint(2*pi*k*F(s));
            %snint1 = (-1/2/i)*(expint(i*2*pi*k*F(s+1)) ...
            %    -expint(-i*2*pi*k*F(s+1)) -expint(i*2*pi*k*F(s)) ...
            %    +expint(-i*2*pi*k*F(s)) );
            % csint1 = cosint(2*pi*k*F(s+1)) - cosint(2*pi*k*F(s)) ;
            csint1 = (-1/2)*(expint(i*2*pi*k*F(s+1))+expint(-i*2*pi*k*F(s+1))...
                -expint(i*2*pi*k*F(s))  -expint(-i*2*pi*k*F(s)) );
            b=b + ( m*snint1 ...
                + b1*2*pi*k.*( -sinc(2*k*F(s+1)) + sinc(2*k*F(s)) + csint1 ))...
                * abs(W((s+1)/2)^2);
            snint1 = sineint(2*pi*F(s+1)*(-I2));
            snint2 = sineint(2*pi*F(s+1)*I1);
            snint3 = sineint(2*pi*F(s)*(-I2));
            snint4 = sineint(2*pi*F(s)*I1);
            G = G - ( ( -1/2*( cos(2*pi*F(s+1)*(-I2))/F(s+1)  ...
                - 2*snint1*pi.*I2 ...
                - cos(2*pi*F(s+1)*I1)/F(s+1) ...
                - 2*snint2*pi.*I1 )) ...
                - ( -1/2*( cos(2*pi*F(s)*(-I2))/F(s)  ...
                - 2*snint3*pi.*I2 ...
                - cos(2*pi*F(s)*I1)/F(s) ...
                - 2*snint4*pi.*I1) ) ) ...
                * abs(W((s+1)/2)^2);
        else      % use usual weights
            m=(M(s+1)-M(s))/(F(s+1)-F(s));
            b1=M(s)-m*F(s);
            b=b+(m/(4*pi*pi)*(sin(2*pi*k*F(s+1))-sin(2*pi*k*F(s)))./(k.*k))...
                * abs(W((s+1)/2)^2) ;
            b = b + (((m*F(s)+b1)*cos(2*pi*k*F(s)) - ...
                (m*F(s+1)+b1)*cos(2*pi*k*F(s+1)))./(2*pi*k)) ...
                * abs(W((s+1)/2)^2) ;
            if need_matrix
                G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))-sinc(2*I2*F(s+1))) ...
                    - .5*F(s)*(sinc(2*I1*F(s))-sinc(2*I2*F(s)))) * ...
                    abs(W((s+1)/2)^2);
            end
        end;
    end

    if need_matrix
        a=G\b;
    else
        a=-4*b*(W(1)^2);
    end
    if Nodd
        h=.5*[flipud(a); 0; -a].';
    else
        h=.5*[flipud(a); -a].';
    end
    if differ, h=-h; end
end

if nargout > 1
    a = 1;
end

%----------------------------------------------------------------------------
function y = sineint(x)
% SINEINT (a.k.a. SININT)   Numerical Sine Integral
%   Used by FIRLS in the Signal Processing Toolbox.
%   Untested for complex or imaginary inputs.
%
%   See also SININT in the Symbolic Toolbox.

%   Was Revision: 1.5, Date: 1996/03/15 20:55:51

i1 = find(real(x)<0);   % this equation is not valid if x is in the
% left-hand plane of the complex plane.
% use relation Si(-z) = -Si(z) in this case (Eq 5.2.19, Abramowitz
%  & Stegun).
x(i1) = -x(i1);
y = zeros(size(x));
ind = find(x);
% equation 5.2.21 Abramowitz & Stegun
%  y(ind) = (1/(2*i))*(expint(i*x(ind)) - expint(-i*x(ind))) + pi/2;
y(ind) = imag(expint(i*x(ind))) + pi/2;
y(i1) = -y(i1);

function [n,msg1,msg2] = firchk(n,Fend,a,exception)
%FIRCHK   Check if specified filter order is valid.
%   FIRCHK(N,Fend,A) checks if the specified order N is valid given the
%   final frequency point Fend and the desired magnitude response vector A.
%   Type 2 linear phase FIR filters (symmetric, odd order) must have a
%   desired magnitude response vector that ends in zero if Fend = 1.  This
%   is because type 2 filters necessarily have a zero at w = pi.
%
%   If the order is not valid, a warning is given and the order
%   of the filter is incremented by one.
%
%   If A is a scalar (as when called from fircls1), A = 0 is
%   interpreted as lowpass and A = 1 is interpreted as highpass.
%
%   FIRCHK(N,Fend,A,EXCEPTION) will not warn or increase the order
%   if EXCEPTION = 1.  Examples of EXCEPTIONS are type 4 filters
%   (such as differentiators or hilbert transformers) or non-linear
%   phase filters (such as minimum and maximum phase filters).

%   Author : R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.5 $  $Date: 2007/12/14 15:15:06 $

error(nargchk(3,4,nargin,'struct'));

if nargin == 3,
    exception = false;
end

msg1 = '';
msg2 = '';
oddord = false; % Flag, initially we assume even order

if isempty(n) || length(n) > 1 || ~isnumeric(n) || ~isreal(n) || n~=round(n) || n<=0,
    msg1 = 'Filter order must be a real, positive integer.';
    return
end

if rem(n,2) == 1,
    oddord = true; % Overwrite flag
end
 
if (a(end) ~= 0) && Fend == 1 && oddord && ~exception,
    str = ['Odd order symmetric FIR filters must have a gain of zero \n'...
     'at the Nyquist frequency. The order is being increased by one.'];
    msg2 = sprintf(str);
    n = n+1;
end

function y=sinc(x)
%SINC Sin(pi*x)/(pi*x) function.
%   SINC(X) returns a matrix whose elements are the sinc of the elements 
%   of X, i.e.
%        y = sin(pi*x)/(pi*x)    if x ~= 0
%          = 1                   if x == 0
%   where x is an element of the input matrix and y is the resultant
%   output element.
%
%   % Example of a sinc function for a linearly spaced vector:
%   t = linspace(-5,5);
%   y = sinc(t);
%   plot(t,y);
%   xlabel('Time (sec)');ylabel('Amplitude'); title('Sinc Function')
%
%   See also SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS,
%   and TRIPULS.

%   Author(s): T. Krauss, 1-14-93
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.1 $  $Date: 2004/08/10 02:11:27 $

i=find(x==0);                                                              
x(i)= 1;      % From LS: don't need this is /0 warning is off                           
y = sin(pi*x)./(pi*x);                                                     
y(i) = 1;   

function w = kaiser(n_est,bta)
%KAISER Kaiser window.
%   W = KAISER(N) returns an N-point Kaiser window in the column vector W.
% 
%   W = KAISER(N,BTA) returns the BETA-valued N-point Kaiser window.
%       If omitted, BTA is set to 0.500.
%
%   See also CHEBWIN, GAUSSWIN, TUKEYWIN, WINDOW.

%   Author(s): L. Shure, 3-4-87
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.17.4.4 $  $Date: 2007/12/14 15:05:16 $

error(nargchk(1,2,nargin,'struct'));

% Default value for the BETA parameter.
if nargin < 2 || isempty(bta), 
    bta = 0.500;
end

[nn,w,trivialwin] = check_order(n_est);
if trivialwin, return, end;

nw = round(nn);
bes = abs(besseli(0,bta));
odd = rem(nw,2);
xind = (nw-1)^2;
n = fix((nw+1)/2);
xi = (0:n-1) + .5*(1-odd);
xi = 4*xi.^2;
w = besseli(0,bta*sqrt(1-xi/xind))/bes;
w = abs([w(n:-1:odd+1) w])';

    
% [EOF] kaiser.m

function [n_out, w, trivalwin] = check_order(n_in)
%CHECK_ORDER Checks the order passed to the window functions.
% [N,W,TRIVALWIN] = CHECK_ORDER(N_ESTIMATE) will round N_ESTIMATE to the
% nearest integer if it is not already an integer. In special cases (N is
% [], 0, or 1), TRIVALWIN will be set to flag that W has been modified.

%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.6.4.2 $  $Date: 2009/05/23 08:16:17 $

w = [];
trivalwin = 0;

if ~(isnumeric(n_in) & isfinite(n_in)),
    error(generatemsgid('InvalidOrder'),'The order N must be finite.');
end

% Special case of negative orders:
if n_in < 0,
   error(generatemsgid('InvalidOrder'),'Order cannot be less than zero.');
end

% Check if order is already an integer or empty
% If not, round to nearest integer.
if isempty(n_in) | n_in == floor(n_in),
   n_out = n_in;
else
   n_out = round(n_in);
   warning(generatemsgid('InvalidOrder'),'Rounding order to nearest integer.');
end

% Special cases:
if isempty(n_out) | n_out == 0,
   w = zeros(0,1);               % Empty matrix: 0-by-1
   trivalwin = 1; 
elseif n_out == 1,
   w = 1;
   trivalwin = 1;   
end

function Y = upfirdn(x,h,varargin)
%UPFIRDN  Upsample, apply a specified FIR filter, and downsample a signal.
%   UPFIRDN(X,H,P,Q) is a cascade of three systems applied to input signal X:
%         (1) Upsampling by P (zero insertion).  P defaults to 1 if not 
%             specified.
%         (2) FIR filtering with the filter specified by the impulse response 
%             given in H.
%         (3) Downsampling by Q (throwing away samples).  Q defaults to 1 if not 
%             specified.
%   UPFIRDN uses an efficient polyphase implementation.
%
%   Usually X and H are vectors, and the output is a (signal) vector. 
%   UPFIRDN permits matrix arguments under the following rules:
%   If X is a matrix and H is a vector, each column of X is filtered through H.
%   If X is a vector and H is a matrix, each column of H is used to filter X.
%   If X and H are both matrices with the same number of columns, then the i-th
%      column of H is used to filter the i-th column of X.
%
%   Specifically, these rules are carried out as follows.  Note that the length
%   of the output is Ly = ceil( ((Lx-1)*P + Lh)/Q ) where Lx = length(X) and 
%   Lh = length(H). 
%
%      Input Signal X    Input Filter H    Output Signal Y   Notes
%      -----------------------------------------------------------------
%   1) length Lx vector  length Lh vector  length Ly vector  Usual case.
%   2) Lx-by-Nx matrix   length Lh vector  Ly-by-Nx matrix   Each column of X
%                                                            is filtered by H.
%   3) length Lx vector  Lh-by-Nh matrix   Ly-by-Nh matrix   Each column of H is
%                                                            used to filter X.
%   4) Lx-by-N matrix    Lh-by-N matrix    Ly-by-N matrix    i-th column of H is
%                                                            used to filter i-th
%                                                            column of X.
%
%   For an easy-to-use alternative to UPFIRDN, which does not require you to 
%   supply a filter or compensate for the signal delay introduced by filtering,
%   use RESAMPLE.
%
%   EXAMPLE: Sample-rate conversion by a factor of 147/160. It is used to
%            % downconvert from 48kHz to 44.1kHz.
%            L = 147; M = 160;                   % Interpolation/decimation factors.
%            Lp = 24;                            % Filter length of each phase
%            N = Lp*L-1;                         % Filter Order
%            h = fir1(N,1/M,kaiser(N+1,7.8562));
%            h = L*h; % Passband gain = L
%            Fs = 48e3;                          % Original sampling frequency: 48kHz
%            n = 0:10239;                        % 10240 samples, 0.213 seconds long
%            x  = sin(2*pi*1e3/Fs*n);            % Original signal, sinusoid at 1kHz
%            y = upfirdn(x,h,L,M);               % 9408 samples, still 0.213 seconds
%
%            % Overlay original (48kHz) with resampled signal (44.1kHz) in red.
%            stem(n(1:49)/Fs,x(1:49)); hold on 
%            stem(n(1:45)/(Fs*L/M),y(12:56),'r','filled'); 
%            xlabel('Time (sec)');ylabel('Signal value');
%    
%   See also RESAMPLE, INTERP, DECIMATE, FIR1, INTFILT, MFILT/FIRSRC in the
%   Filter Design Toolbox.
  
%   Author(s): Paul Pacheco
%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.6.4.5 $  $Date: 2008/09/13 07:14:26 $

%   This M-file validates the inputs, sets defaults, and then calls the C MEX-file.

% Validate number of I/O args.
error(nargchk(2,4,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

% Force to be a column if input is a vector
[mx,nx] = size(x);
if find([mx nx]==1),
  x = x(:);  % columnize it.
end
[Lx,nChans] = size(x);

% Force to be a column if filter is a vector
if find(size(h)==1),
    h = h(:);  % columnize it.
end
[Lh,hCols] = size(h);

% Validate input args and define defaults.
[p,q,errid,errmsg] = validateinput(x,h,varargin);
if ~isempty(errmsg), error(errid,errmsg); end

% Call the MEX-file
Y = upfirdnmex(x,h,p,q,Lx,Lh,hCols,nChans);

% Convert output to be a row vector (if x was a row and H is NOT a matrix)
if (mx==1) && (hCols == 1)
    Y = Y(:).';
end


%----------------------------------------------------------------------
function [p,q,errid,errmsg] = validateinput(x,h,opts)

% Default values
p = 1;
q = 1;
errid = '';
errmsg = '';

% Validate 1st two input args: signal and filter.
if isempty(x) || issparse(x) || ~isa(x,'double'),
    errid = generatemsgid('invalidInput');
    errmsg = 'The input signal X must be a double-precision vector.';
    return;
end
if isempty(h) || issparse(h) || ~isa(h,'double'),
    errid = generatemsgid('invalidFilter');
    errmsg = 'The filter H must be a double-precision vector.';
    return;
end

% The following check is for case 4 (as seen on the reference page), i.e., 
% x and h are matrices, check that they both have the same number of
% columns. 
nChans = size(x, 2);
hCols  = size(h, 2);
if (nChans > 1) && (hCols > 1) && (hCols ~= nChans),
    errid = generatemsgid('xNhSizemismatch');
    errmsg = 'Signal X and filter H must have the same number of columns.';
    return;
end

% Validate optional input args: upsample and downsample factors.
nopts = length(opts);
if (nopts >= 1),
    p = opts{1};
    if isempty(p) || ~isa(p,'double') || p<1 || ~isequal(round(p),p),
        errid = generatemsgid('invalidP');
        errmsg = 'The upsample factor P must be a positive, double-precision, integer.';
        return;

    elseif (nopts == 2),
        q = opts{2};
        if isempty(q) || ~isa(q,'double') || q<1 || ~isequal(round(q),q),
            errid = generatemsgid('invalidQ');
            errmsg = 'The downsample factor Q must be a positive, double-precision, integer.';
            return;
        end
    end
    if p*q > intmax('int32'),
        errid = generatemsgid('ProdPNQTooLarge');
        errmsg = ['The product of the downsample factor Q and the upsample factor P must be',...
            ' less than 2^31.'];
        return;
    end
end

% [EOF]


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
