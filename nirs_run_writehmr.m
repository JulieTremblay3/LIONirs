function out = nirs_run_writehmr(job)
%JT
fileout = job.FileOutput;
% Matrice NIRS

%Big loop over all subjects
%% Loop over all nirs.mat file defined 
tic
scaling = 1 %1e-6;
for ifile=1:numel(job.NIRSmat)   %each session defined in a datafile
    load(job.NIRSmat{ifile});
    
    
   currentsub = 1;
   [path,name,ext]= fileparts(NIRS.Dt.fir.pp(1).p{1});
  
   DOT{currentsub}.subjectNum = {name(1:end-4)};
   %SD structure 
   DOT{currentsub}.SD.Lambda = NIRS.Cf.dev.wl;
   DOT{currentsub}.SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
   DOT{currentsub}.SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
   DOT{currentsub}.SD.nSrcs  = numel(DOT{currentsub}.SD.SrcPos)/3;
   DOT{currentsub}.SD.nDets  = numel(DOT{currentsub}.SD.DetPos)/3;
   DOT{currentsub}.SD.SrcAmp = ones(DOT{currentsub}.SD.nSrcs ,2);
   DOT{currentsub}.SD.DetAmp = ones(DOT{currentsub}.SD.nDets,2);
   %SD
   DOT{currentsub}.color = lines(size(NIRS.Cf.H.C.id,2));
   
   DOT{currentsub}.currentFile = 1;
%    if ~isfield(DOT{currentsub},'data')
%     DOT{currentsub}.data = [];  
%    end
   cf = ifile;

    %DATA ********************************************************************** 
    
    DOT{currentsub}.data(cf).FILEpathnm = [path,'\'];
        %Load module epoch average info
        if isfield(NIRS.Dt.fir.pp(end),'chok')
            MeasListAct = NIRS.Dt.fir.pp(end).chok;    
        else
             MeasListAct = NIRS.Cf.H.C.ok;
        end

        for imodule = numel(NIRS.Dt.fir.pp):-1:1
                if ~isempty(strfind(NIRS.Dt.fir.pp(imodule).pre,'Epoch averaging'))&isempty(strfind(NIRS.Dt.fir.pp(imodule).pre,'Manual Gui'))
                   module_Epoch = imodule;
                   break
                end
        end
        MeasListAct = NIRS.Dt.fir.pp(module_Epoch).chok;    
        pretime = -abs(NIRS.Dt.fir.pp(module_Epoch).job.choiceave.pretime);
        posttime = NIRS.Dt.fir.pp(module_Epoch).job.choiceave.posttime;
        [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(module_Epoch).p{1});     
        [rem,tok]=strtok(fliplr(pathstr),{'\','/'});
        pathout = fliplr(tok);
        %Fin Load module epoch average info
        
        for imodule = numel(NIRS.Dt.fir.pp):-1:1
        if ~isempty(strfind(NIRS.Dt.fir.pp(imodule).pre,'Filtered'))             
               lpf = NIRS.Dt.fir.pp(imodule).job.lowcutfreq;
               hpf = NIRS.Dt.fir.pp(imodule).job.highcutfreq;
               break
            end
        end
        if ~isfield(DOT{currentsub}.data(cf),'lpf')
            lpf = [0.1];
            hpf = [0];
        end
        %LOAD RAW epoch average trial
        trialname = fullfile(pathstr,[name,'.mat']);
        data = load(trialname,'-mat');
        numStim = size(data.A,3);
       
%         if strcmp(NIRS.Dt.fir.pp(end).pre(1:27),'Epoch averaging Help memory')
        %HELP MEMORY AVG
%             raw = ones(numel(MeasListAct)*2,numel(MeasListAct))
        if strfind(NIRS.Dt.fir.pp(module_Epoch).pre,'Epoch averaging') 
            raw = reshape(data.A,size(data.A,1), size(data.A,2)*size(data.A,3))';
        end
        
        
        idNaN = find(isnan(raw));
        if ~isempty(idNaN)
            noise = zeros(size(raw));
            noise(idNaN) = 4;
            raw(idNaN)=0; %Pour éviter les problèmes de compatibilité dans homer min max... 
        else
            noise = zeros(size(raw));
        end
        stim = zeros(size(raw,1),1);
        
        %Load HRF data
       idHbO = 1:NIRS.Cf.H.C.N/2;
       idHbR = NIRS.Cf.H.C.N/2+1:NIRS.Cf.H.C.N;
       %LOAD  HRF.AvgC average epoch
        xname =fullfile(pathstr,[name,ext]);     
        d1 = fopen_NIR(xname,NIRS.Cf.H.C.N)'; 
        idNaN = find(isnan(d1));       
        if ~isempty(idNaN);d1(idNaN)=0;end       
       
        xname = fullfile(pathstr,[name,'_std',ext]);
        d1std = fopen_NIR(xname,NIRS.Cf.H.C.N)';
        idNaN = find(isnan(d1std));        
        if ~isempty(idNaN);d1std(idNaN)=0;end       
        
        xname = fullfile(pathstr,[name,'_events',ext]);
        d1event = fopen_NIR(xname,NIRS.Cf.H.C.N)';     
        idNaN = find(isnan(d1event)|d1event==0);
        if ~isempty(idNaN);d1event(idNaN)=1;end
    
       %trig artificiel sur epoch average pre et post time      
       tHRF = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(d1,1)*1/NIRS.Cf.dev.fs;
       tHRF= tHRF-abs(pretime); 
       id = find(tHRF <= 0);
       t =  1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(raw,1)*1/NIRS.Cf.dev.fs;
       stim = zeros(numel(t),1);
       trigid =id(end):numel(tHRF):((numel(tHRF)*numStim));
       stim(trigid)=1;

    DOT{currentsub}.data(cf).filenm = fliplr(rem);
    DOT{currentsub}.data(cf).BaselineData = 0;
    DOT{currentsub}.data(cf).UserStim = num2str(t(trigid));
    DOT{currentsub}.data(cf).Toggle_UserStim = 1;
    DOT{currentsub}.data(cf).lpf = lpf;
    DOT{currentsub}.data(cf).nSV = 0;
    DOT{currentsub}.data(cf).hpf = hpf;
    DOT{currentsub}.data(cf).nSV_dOD = [0 0];


   DOT{currentsub}.data(cf).raw = raw.*scaling;
   DOT{currentsub}.data(cf).Aux = zeros(size(raw,1),6);
   DOT{currentsub}.data(cf).t = t;
   DOT{currentsub}.data(cf).MeasList = [NIRS.Cf.H.C.id(2:3,:)',...
       ones(size(NIRS.Cf.H.C.id,2),1),...
    [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
   DOT{currentsub}.data(cf).MeasListini = DOT{currentsub}.data(cf).MeasList;     
   DOT{currentsub}.data(cf).AuxLabels = {'aux channel 1 ','aux channel 2 ','aux channel 3 ','aux channel 4 ','aux channel 5 ','aux channel 6 '};
   DOT{currentsub}.data(cf).phase = [];
   DOT{currentsub}.data(cf).stimRaw = stim;
   DOT{currentsub}.data(cf).nMeas = numel(DOT{currentsub}.data(cf).MeasList)/4;
   DOT{currentsub}.data(cf).Sact = ones(1,DOT{currentsub}.SD.nSrcs);
   DOT{currentsub}.data(cf).Dact = ones(1,DOT{currentsub}.SD.nDets);
   DOT{currentsub}.data(cf).MeasListAct = double(MeasListAct);
   DOT{currentsub}.data(cf).stimOrig = stim;
   DOT{currentsub}.data(cf).numStim = numStim;
   DOT{currentsub}.data(cf).stimOn =stim;
   %SECTION HRF
   
   
   
   DOT{currentsub}.data(cf).HRF.regdata = stim;
   DOT{currentsub}.data(cf).HRF.StimOn = stim;
   DOT{currentsub}.data(cf).HRF.numStim = size(data.A,3);
   DOT{currentsub}.data(cf).HRF.type = 'User Defined';
   DOT{currentsub}.data(cf).HRF.pretime =  pretime;
   DOT{currentsub}.data(cf).HRF.posttime = posttime;
   DOT{currentsub}.data(cf).HRF.lpf = lpf;
   DOT{currentsub}.data(cf).HRF.hpf = hpf;
   DOT{currentsub}.data(cf).HRF.StoredRegData = stim;
   DOT{currentsub}.data(cf).HRF.Type = 'UserDefined';
   DOT{currentsub}.data(cf).HRF.stimOrig = stim;
   DOT{currentsub}.data(cf).HRF.UserStim = num2str(t(trigid));
   DOT{currentsub}.data(cf).HRF.tHRF = tHRF;
   DOT{currentsub}.data(cf).HRF.Conf_High = [];
   DOT{currentsub}.data(cf).HRF.Conf_Low = [];
   DOT{currentsub}.data(cf).HRF.Conf_High_Conc = [];
   DOT{currentsub}.data(cf).HRF.Conf_High_Conc = [];  
    
    DOT{currentsub}.data(cf).HRF.Avg = zeros(size(d1));
    DOT{currentsub}.data(cf).HRF.AvgStd = zeros(size(d1));
    DOT{currentsub}.data(cf).HRF.AvgStdErr = zeros(size(d1));
    DOT{currentsub}.data(cf).HRF.AvgOdd = zeros(size(d1));
    DOT{currentsub}.data(cf).HRF.AvgEven = zeros(size(d1)); 
    
    DOT{currentsub}.data(cf).HRF.AvgC = [];
    DOT{currentsub}.data(cf).HRF.AvgC(:,:,1) = d1(:,idHbO).*scaling;
    DOT{currentsub}.data(cf).HRF.AvgC(:,:,2) = d1(:,idHbR).*scaling;
    DOT{currentsub}.data(cf).HRF.AvgC(:,:,3) = d1(:,idHbO).*scaling+d1(:,idHbR).*scaling;
    %LOAD HRF std average epoch
   
    DOT{currentsub}.data(cf).HRF.AvgCStd(:,:,1) = d1std(:,idHbO).*scaling;
    DOT{currentsub}.data(cf).HRF.AvgCStd(:,:,2) = d1std(:,idHbR).*scaling;
    DOT{currentsub}.data(cf).HRF.AvgCStd(:,:,3) = zeros(size(d1std(:,idHbR)));%HbT ici
    
    
    DOT{currentsub}.data(cf).HRF.AvgCStdErr(:,:,1) = d1std(:,idHbO).*scaling./ d1event(:,idHbO).^0.5;
    DOT{currentsub}.data(cf).HRF.AvgCStdErr(:,:,2) = d1std(:,idHbO).*scaling./ d1event(:,idHbO).^0.5;
    DOT{currentsub}.data(cf).HRF.AvgCStdErr(:,:,3) = zeros(size(d1std(:,idHbR))).*scaling;%HbT ici
    
    %DOT{currentsub}.data(cf).HRF.STDBaseline  
    d1event = fopen_NIR(xname,NIRS.Cf.H.C.N)';
    DOT{currentsub}.data(cf).HRF.navg = nanmean(d1event(:,idHbO)); %en moyenne le nombre d'évènement retenue
    %FIN SECTION HRF    
    DOT{currentsub}.data(cf).nSV_dConc = [0 0];
    DOT{currentsub}.data(cf).dOD_UpToDate = 0;
    DOT{currentsub}.data(cf).svs = 0; %
    DOT{currentsub}.data(cf).Intens_hpf =[];
    DOT{currentsub}.data(cf).norm =[];
    DOT{currentsub}.data(cf).norm_cov =[];
    DOT{currentsub}.data(cf).dOD =[];
    DOT{currentsub}.data(cf).dODc = raw;
    DOT{currentsub}.data(cf).dConc(:,:,1)=raw(:,idHbO);  
    DOT{currentsub}.data(cf).dConc(:,:,2)=raw(:,idHbR);
    DOT{currentsub}.data(cf).dConc(:,:,3)=raw(:,idHbO)+raw(:,idHbR);
    DOT{currentsub}.data(cf).dConcc(:,:,1)=raw(:,idHbO);  
    DOT{currentsub}.data(cf).dConcc(:,:,2)=raw(:,idHbR);
    DOT{currentsub}.data(cf).dConcc(:,:,3)=raw(:,idHbO)+raw(:,idHbR);
    DOT{currentsub}.data(cf).svsConc = 0;
    DOT{currentsub}.data(cf).avgUpToDate = 1;
    DOT{currentsub}.data(cf).MeasListActNoisy = MeasListAct;
    DOT{currentsub}.data(cf).nkp = 0;
    DOT{currentsub}.data(cf).Thresh = 0.2; 
    DOT{currentsub}.data(cf).Tsep = 2; 
    DOT{currentsub}.data(cf).HPF = [];
    
   
    

    
   DOT{currentsub}.data(cf).noise=noise;
    
   DOT{currentsub}.prj_name = job.prjfile{1}; 
   DOT{currentsub}.plot = [1,1];
   DOT{currentsub}.plotLst = [1];
 % DOT{currentsub}.prune
 % DOT{currentsub}.zone
   DOT{currentsub}.imgAvgLen = 5;
   
   DOT{currentsub}.AdvOptions.DispOptions.Fig1ShowModel     = 0;
   DOT{currentsub}.AdvOptions.DispOptions.Fig3ShowStdErr    = 0;
   DOT{currentsub}.AdvOptions.DispOptions.Fig3ShowStdDev    = 1;
   DOT{currentsub}.AdvOptions.DispOptions.Fig3ShowEBars     = 0;
   DOT{currentsub}.AdvOptions.DispOptions.Fig3ShowConf      = 0;
   DOT{currentsub}.AdvOptions.DispOptions.Fig3ShowOddEven   = 0;
   
   DOT{currentsub}.AdvOptions.FiltOptions.PartialVolumeCorr = [50 50];
   DOT{currentsub}.AdvOptions.FiltOptions.UsePartialVolumeCorr = 0;
   DOT{currentsub}.AdvOptions.FiltOptions.DPF = [6,6];
   DOT{currentsub}.AdvOptions.FiltOptions.UsePathLenghtCorr = 0;
   
   DOT{currentsub}.AdvOptions.FiltOptions.UseCustomFilter = 0;
   DOT{currentsub}.AdvOptions.FiltOptions.Wiener = 0;
   DOT{currentsub}.AdvOptions.FiltOptions.FilterType = [1];
   DOT{currentsub}.AdvOptions.FiltOptions.Detrend = 0;
   DOT{currentsub}.AdvOptions.FiltOptions.Detrendsize = 30;
   DOT{currentsub}.AdvOptions.FiltOptions.Viewnoise = 0;
   DOT{currentsub}.AdvOptions.FiltOptions.ppod = 1;
   DOT{currentsub}.AdvOptions.FiltOptions.FilterOrder = 3;
   DOT{currentsub}.AdvOptions.FiltOptions.AbsSpect = 1;
   
   DOT{currentsub}.AdvOptions.AdvOptions.UseIterReg = 0;
   DOT{currentsub}.AdvOptions.AdvOptions.UseDetrend = 0;
   
   DOT{currentsub}.AdvOptions.NormOptions.Normalize = 0;
   DOT{currentsub}.AdvOptions.NormOptions.Custom = 0;
   DOT{currentsub}.AdvOptions.NormOptions.UserNormalize = [];
end %end for  Big loop over subjects
            [tok,rem]=strtok(job.FileOutput,':');
            if isempty(rem) %si le nom du repertoire est juste un nom à ajouter dans le fichier de donné
                     save([pathout,job.FileOutput ,'.hmr'],'DOT','-MAT');
                    fprintf('%s\n',['File created : ',pathout,job.FileOutput ,'.hmr'])

            else %si le nom du repertoire est un répertoire complet
                try
                     label.pathout = [job.FileOutput,filesep];
                     save([job.FileOutput ,'.hmr'],'DOT','-MAT');
                     fprintf('%s\n',['File created : ',job.FileOutput ,'.hmr'])

                catch
                    msgbox(['The name of the output directory : ',job.FileOutput,' is invalid'])
                end
            end

    toc
    out.NIRSmat = job.NIRSmat;