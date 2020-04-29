function out = nirs_run_segment(job)
%Use to segment NIRS event around trigger
%Synchronised according to mutual trigger information EEG, AUX and Video Field
%
%filename prefix 
prefix = 's'; %for "segment"
DelPreviousData  = job.DelPreviousData;
flagsegmentAUX = 0; %write offset time instead off create a copy of the data segmented 0 keep whole original segment 1 segment
flagsegmentEEG = 0;
for filenb = 1:size(job.NIRSmat,1)
    NIRS = [];
    load(job.NIRSmat{filenb,1});
    [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1})
    NIRSold = NIRS;
    %use last step of operation
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    fprintf('%s\n','File processed');
    if isfield(NIRS.Dt,'EEG')
        moduleEEG =  numel(NIRS.Dt.EEG(end).pp);
    end
    if isfield(NIRS.Dt,'AUX')
        moduleaux = numel(NIRS.Dt.AUX(end).pp);
    end
    if isfield(NIRS.Dt,'Video')
        moduleVideo = numel(NIRS.Dt.Video(end).pp);
    end
    if isfield(NIRS.Dt,'Audio')
        moduleAudio = numel(NIRS.Dt.Audio(end).pp);
    end
    ifile = 1; %Utiliser si on remets chaque stim normaliser dans des fichier différent
    
    for f=1:numel(rDtp) %Loop over all files of a NIRS.mat
        d = fopen_NIR(rDtp{f},NC);     
       [dir1,fil1,ext1] = fileparts(rDtp{f});
        vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
        [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
        nsample = size(d,2);
         noise = ind_dur_ch2mat(ind_dur_ch, nsample,NC)';
  
        if isfield(NIRS.Dt,'EEG')
            try
                fileeeg=NIRS.Dt.EEG.pp(moduleEEG).p{f}; %OPEN EEG HERE
                [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]=fopen_EEG(fileeeg);
            catch
            end
            
          if 1
              figure 
              hold on
              actualNIRStrig = NIRS.Dt.fir.aux5{1};
             for itrig=1:size(actualNIRStrig,1)             
                 plot([actualNIRStrig(itrig,2)*1/fs actualNIRStrig(itrig,2)*1/fs],[0,1],'k','displayname',['NIRS S ',num2str(actualNIRStrig(itrig,1))]);
             end
             for itrig=1:size(EEG.marker,1) 
                 plot([EEG.ind_dur_ch(itrig,1)*EEG.infoBV.SamplingInterval/1000000 EEG.ind_dur_ch(itrig,1)*EEG.infoBV.SamplingInterval/1000000],[0,1],'r','displayname',['EEG ',EEG.marker{itrig,2}]);
             end
          end
        end 
        
        try
        if isfield(NIRS.Dt,'AUX')            
            for iaux = 1:numel(NIRS.Dt.AUX)
                fileaux = NIRS.Dt.AUX(iaux).pp(moduleaux).p{f}; %OPEN EEG HERE
                [AUX(iaux).data,AUX(iaux).infoBV,AUX(iaux).marker,AUX(iaux).ind_dur_ch]=fopen_EEG(fileaux);
            end
        end
        catch
            msgbox('AUX could not be segment')
        end
        
        trigger = job.trigger;
        switch job.pretime
            case 'start'
                  pretime = 1;  
            otherwise
                 pretime = str2num(job.pretime); %in seconde
        end
        
       
        switch job.posttime
            case 'end'
                 posttime = size(d,2);  
            otherwise
                 posttime = str2num(job.posttime); %in seconde %round(fs*str2num(job.posttime)); in sample
        end
        
        NIRS.Dt.fir.aux5bloc = NIRS.Dt.fir.aux5{f};
        aux5 = NIRS.Dt.fir.aux5{f};
            indstim = [];
            itrigger = [];
        if trigger == 0 %If 0, consider every trigger
            indstim = aux5(:,2);
            ind = [indstim; size(d,2)];
            itrigger = 1;
        else
        
            for itypestim = 1:numel(trigger)
                idstim= aux5((aux5(:,1) == trigger(itypestim)),2)
                indstim  = [indstim  ;idstim ];
                ind = [1; idstim; size(d,2)];
                itrigger =[itrigger, ones(1,numel(idstim)).*trigger(itypestim)];
            end
            if isfield(NIRS.Dt,'EEG')
                try
                    tmp = EEG.marker{2,2};
                    idstimEEG = [];
                    tmp = EEG.marker{2,2};
                    if iscell(tmp)
                        tmp = tmp{1};
                    end
                    if isempty(tmp)
                        tmp = {'nul'};
                    end
                    for itypestim = 1:numel(trigger)
                        if strcmp(tmp(1),'S') %S suffix nirx but not ISS
                            idtrigEEG = ['S',sprintf('%3.0f',trigger(itypestim))];
                        else
                            idtrigEEG = [sprintf('%3.0f',trigger(itypestim))];
                        end
                        
                        idstimEEG = [idstimEEG;strmatch(idtrigEEG,EEG.marker(:,2))];
                    end
                      if isfield(NIRS.Dt.EEG.pp(moduleEEG),'sync_timesec')
                        tstart= NIRS.Dt.EEG.pp(moduleEEG).sync_timesec{f}
                        tstop = tstart+size(d,2)*1/NIRS.Cf.dev.fs                    
                        valtimestim =       EEG.ind_dur_ch(idstimEEG,1).*EEG.infoBV.SamplingInterval/1000000
                        triginside = find(valtimestim>tstart & valtimestim<tstop)                   
                         tmp = idstimEEG ;
                         idstimEEG = tmp(triginside);
                    end
                    
                    idtrigEEG = ['S',sprintf('%3.0f',trigger)];
                    if numel(indstim)==numel(idstimEEG) %gerer les cas stim agree dimention only
                        indstim_EEG = EEG.ind_dur_ch(idstimEEG,1);
                    else
                        msgbox('unequal EEG trigger identification')
                        indstim_EEG = EEG.ind_dur_ch(idstimEEG(1:5),1);
                    end
                    fsEEG = 1/(EEG.infoBV.SamplingInterval/1000000); %uS en seconde
                 
                    pretimeEEG = round(fsEEG*pretime)-1;
                    posttimeEEG = round(fsEEG*posttime);
                catch
                end
            end
            %   try
            if isfield(NIRS.Dt,'AUX')
                %ISS SYSTEM or create trig via GUI_AUXedit USE
                %'trigger' as marker name column 1
                %NIRX
                for iaux = 1:numel(NIRS.Dt.AUX)
                    tmp = AUX(iaux).marker{end,2};
                    idstimAUX{iaux} = [];
                    if iscell(tmp)
                        tmp = tmp{1};
                    end
                    if isempty(tmp)
                        tmp = {'nul'};
                    end
                    for itypestim = 1:numel(trigger)
                        if strcmp(tmp(1),'S') %S suffix nirx but not ISS
                            idtrigAUX{iaux} = ['S',sprintf('%3.0f',trigger(itypestim))];
                        else
                            idtrigAUX{iaux} = [sprintf('%3.0f',trigger(itypestim))];
                        end
                        
                        idstimAUX{iaux} = [idstimAUX{iaux};strmatch(deblank(idtrigAUX{iaux}),deblank(AUX(iaux).marker(:,2)))];
                    end
                    %PRESEGMENTATION ALLREALY DONE.
                    if isfield(NIRS.Dt.AUX(iaux).pp(moduleaux),'sync_timesec')
                    tstart= NIRS.Dt.AUX(iaux).pp(moduleaux).sync_timesec{f}
                    tstop = tstart+size(d,2)*1/NIRS.Cf.dev.fs
                    indstim*1/NIRS.Cf.dev.fs
                    valtimestim =        AUX(iaux).ind_dur_ch(idstimAUX{iaux},1).*AUX(iaux).infoBV.SamplingInterval/1000000
                    triginside = find(valtimestim>tstart & valtimestim<tstop)                   
                    tmp = idstimAUX{iaux} ;
                    idstimAUX{iaux} = tmp(triginside);
                    end
                    
                     %itrigger =[itrigger, ones(1,numel(idstim)).*trigger(itypestim)];
                    if numel(indstim)==numel(idstimAUX{iaux}) %gerer les cas stim agree dimention only
                        indstim_AUX{iaux} = AUX(iaux).ind_dur_ch(idstimAUX{iaux},1);
                    else
                        idtocheckAUX = AUX(iaux).ind_dur_ch(idstimAUX{iaux},1)
                        timetrig = idtocheckAUX * AUX(iaux).infoBV.SamplingInterval/1000000
                        timestim = indstim*1/NIRS.Cf.dev.fs
                        figure;hold on;
                        for i=1:numel(timetrig)
                            plot( timetrig(i),1,'+r','displayname',['AUX',num2str(i)])
                        end                 
                        for i=1:numel(timestim)
                            plot(timestim(i),1,'+b','displayname',['NIRS',num2str(i)])
                        end
                        
                        title(['Error in number of trig ',num2str(trigger) ,' between AUX ', NIRS.Dt.AUX(iaux).label,' and NIRS, Please Ajust trig manualy or in VMRK file to make them match before segmentation.'])
                        return
                    end
                    
                    fsAUX{iaux} = 1/(AUX(iaux).infoBV.SamplingInterval/1000000); %uS en seconde
                    pretimeAUX{iaux} = round(fsAUX{iaux}*pretime)-1; 
                    posttimeAUX{iaux} = round(fsAUX{iaux}*posttime);
                    
                end
            end
%              catch
%             msgbox('AUX could not be segment')
%         end
        end
        
        
    
        t = 1/fs:1/fs:1/fs*size(d,2);
        if ~isempty(indstim)
        for istim = 1:numel(indstim)
             %dnorm = zeros(NC,nstim);
            if find((indstim(istim)-pretime) <= 0)
                 disp('Could not be process : Segment pretime before stim out of range use 0 or negative value to take from the trig or la ')
                 return
                % indstim(:) = indstim(:) + size(d,2) ;
                % d =[fliplr(d),d];
            end
           
%             if (max(indstim)+posttime)>size(d,2)
%                  switch job.posttime
%                  case 'end'
%                       posttime = size(d,2)-max(indstim);  
%                      otherwise
%                         msgbox('Warning out of range padding in the last bloc')
%                         d = [d,fliplr(d)]; 
%                         noise = [noise,fliplr(noise)]; 
%                 end
%                    %round(fs*str2num(job.posttime));
%             end
%                switch job.pretime
%                  case 'start'
%                       pretime = indstim(istim)-1;  
%                       if 0
%                             pretimeAUX{iaux} = round(fsAUX{iaux}*pretime)-1;
%                             posttimeAUX{iaux} = round(fsAUX{iaux}*posttime);
%                       end
%                  otherwise
%              	spretime = 
%                end
              
               
          %  nstim = (sposttime + spretime + 1);   
            %Segment the data NIRS
            dnorm = d(:,indstim(istim)-spretime:indstim(istim)+sposttime);
            dnoise = noise(:,indstim(istim)-spretime:indstim(istim)+sposttime);
            %figure;imagesc(dnoise)
            
            
            
            bloc = sprintf('%s%02.0f','b',istim);
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            if 1
                outfile = fullfile(dir2,[prefix fil1 bloc ext1]);
                fileOutRoot_vmrk= fullfile(dir2,[prefix fil1 bloc '.vmrk']);
                fileOutRoot_vhdr = fullfile(dir2,[prefix fil1 bloc '.vhdr']) ;
                if isfield(NIRS.Dt,'EEG')
                    %WRITE IN SUB FILE.... or add sync_timesec
                    [dirEEG,filEEG,extEEG]= fileparts(NIRS.Dt.EEG.pp(end).p{f});
                    
                    outfileEEG = fullfile(dirEEG,[filEEG bloc '.dat']);
                    stimlist=indstim_EEG;
                    %write and segment EEG file in homologous bloc
                    if (max(stimlist)+posttimeEEG)>size(EEG.data,1)
                        msgbox('Warning out of range padding in the last bloc aux')
                        EEG.data= [EEG.data;flipud(EEG.data)];
                        %  figure;plot(AUX(iAUX).data)
                    end
                                       
                    if flagsegmentEEG
                        datasegR = EEG.data(indstim_EEG(istim)-pretimeEEG:indstim_EEG(istim)+posttimeEEG,:);
                        fwrite_EEG(outfileEEG,EEG,indstim_EEG(istim)-pretimeEEG,indstim_EEG(istim)+posttimeEEG);
                        bloc                       
                        NIRS.Dt.EEG.pp(moduleEEG+1).p{ifile,1}=outfileEEG;                            
                    else
                        NIRS.Dt.EEG.pp(moduleEEG+1).p{ifile} = NIRS.Dt.EEG.pp(moduleEEG).p{f};
                        NIRS.Dt.EEG.pp(moduleEEG+1).sync_timesec{ifile} = (indstim_EEG(istim)-pretimeEEG)*1/fsEEG;                        
                    end
                    
                    
                end
                if isfield(NIRS.Dt,'Video')
                    if ~isempty(NIRS.Dt.Video.pp(moduleVideo).offset)
                        offset = NIRS.Dt.Video.pp(moduleVideo).offset{f};
                    else
                        offset = 0;
                    end
                    switch  NIRS.Dt.Video.syncref                       
                        case 'EEG'
                            try
                            NIRS.Dt.Video.pp(moduleVideo+1).p{ifile} = NIRS.Dt.Video.pp(moduleVideo).p{f};
                            NIRS.Dt.Video.pp(moduleVideo+1).sync_timesec{ifile} = (indstim_EEG(istim)-pretimeEEG)*1/fsEEG + offset;
                            catch
                                msgbox('Video segmentation could not be done using EEG trigger')
                                return
                            end
                     
                        case 'NIRS'
                            NIRS.Dt.Video.pp(moduleVideo+1).p{ifile} = NIRS.Dt.Video.pp(moduleVideo).p{f};
                            NIRS.Dt.Video.pp(moduleVideo+1).sync_timesec{ifile} = (indstim(istim)-pretime)*1/fs  + offset;
                         case 'AUX'
                            NIRS.Dt.Video.pp(moduleVideo+1).p{ifile} = NIRS.Dt.Video.pp(moduleVideo).p{f};
                            indstim_AUXval = indstim_AUX{1}; %synchronisation will refere to file 1 
                            pretimeAUXval =  round(fsAUX{iaux}*pretime)-1; %pretimeAUX{1};
                            fsAUXval = fsAUX{1};
                            
                            NIRS.Dt.Video.pp(moduleVideo+1).sync_timesec{ifile} = (indstim_AUXval(istim)- pretimeAUXval)*1/fsAUXval  + offset;
                    end
                    
                end
                if isfield(NIRS.Dt,'Audio')
                     offset = NIRS.Dt.Audio.pp(moduleAudio).offset{f};
                    switch  NIRS.Dt.Audio.syncref                       
                        case 'EEG'
                            try
                            NIRS.Dt.Audio.pp(moduleAudio+1).p{ifile} = NIRS.Dt.Audio.pp(moduleAudio).p{f};
                            NIRS.Dt.Audio.pp(moduleAudio+1).sync_timesec{ifile} = (indstim_EEG(istim)-pretimeEEG)*1/fsEEG + offset;
                            catch
                                msgbox('Video segmentation could not be done using EEG trigger')
                                return
                            end
                     
                        case 'NIRS'
                            NIRS.Dt.Audio.pp(moduleAudio+1).p{ifile} = NIRS.Dt.Audio.pp(moduleAudio).p{f};
                            NIRS.Dt.Audio.pp(moduleAudio+1).sync_timesec{ifile} = (indstim(istim)-pretime)*1/fs  + offset;
                         case 'AUX'
                            NIRS.Dt.Audio.pp(moduleAudio+1).p{ifile} = NIRS.Dt.Audio.pp(moduleAudio).p{f};
                            indstim_AUXval = indstim_AUX{1}; %synchronisation will refere to file 1 
                            pretimeAUXval = pretimeAUX{1};
                            fsAUXval = fsAUX{1};
                            
                            NIRS.Dt.Audio.pp(modulemoduleAudio+1).sync_timesec{ifile} = (indstim_AUXval(istim)- pretimeAUXval)*1/fsAUXval  + offset;
                    end
                    
                end
               %round(fsAUX{iaux}*pretime)-1
                try
                if isfield(NIRS.Dt,'AUX')
                    for iAUX = 1:numel(NIRS.Dt.AUX)
                        
                        [dirAUX,filAUX,extAUX]= fileparts(NIRS.Dt.AUX(iAUX).pp(moduleaux).p{f});
                        outfileAUX{iAUX} = fullfile(dirAUX,[filAUX bloc '.dat']);
                        %write and segment EEG file in homologous bloc
                        stimlist=indstim_AUX{iAUX}
                        if (max(stimlist)+posttimeAUX{iAUX})>size(AUX(iAUX).data,1)
                            msgbox('Warning out of range padding in the last bloc aux')
                            AUX(iAUX).data= [AUX(iAUX).data;flipud(AUX(iAUX).data)];
                            %  figure;plot(AUX(iAUX).data)
                        end
                        if flagsegmentAUX
                            datasegR = AUX(iAUX).data(stimlist(istim)-pretimeAUX{iAUX}:stimlist(istim)+posttimeAUX{iAUX},:);
                            fwrite_EEG(outfileAUX{iAUX},AUX(iAUX),stimlist(istim)-pretimeAUX{iAUX},stimlist(istim)+posttimeAUX{iAUX});
                            NIRS.Dt.AUX(iAUX).pp(moduleaux+1).p{ifile,1}=outfileAUX{iAUX};
                        else
                            NIRS.Dt.AUX(iAUX).pp(moduleaux+1).p{ifile,1}=NIRS.Dt.AUX(iAUX).pp(moduleaux).p{f};
                            NIRS.Dt.AUX(iAUX).pp(moduleaux+1).sync_timesec{ifile,1} = (stimlist(istim)-pretimeAUX{iAUX})*1/fsAUX{iAUX};
                        end
                        bloc;
                    end
                    
                end
                catch
                        msgbox('AUX could not be segment')
                    end
                %dnorm=log(dnorm)
                %WRITE THE RESEGMENTATION
                NIRS.Dt.fir.sizebloc{ifile}=size(dnorm,2)
                fwrite_NIR(outfile,dnorm);
                
                
                
                rDtptmp{ifile} = rDtp{f};   %Nom temporaire pour gere la division du fichier et garder la trace du fichier original
                auxraw{ifile} = NIRS.Dt.fir.aux5{f}
                chok(:,ifile) = NIRS.Cf.H.C.ok(:, f);
                %reajuster tout les trigs between pretime and postime of the
                %selected segment
                aux5all = NIRS.Dt.fir.aux5{f}
                idtrigtokeep = find(aux5all(:,2)>(indstim(istim) - pretime) & aux5all(:,2)<(indstim(istim) +posttime))
                if ~isempty( idtrigtokeep)
                    aux5new(:,1) = aux5all(idtrigtokeep,1);
                    aux5new(:,2) = aux5all(idtrigtokeep,2) - indstim(istim) +pretime;
                    aux5temp{ifile} = aux5new
                    clear aux5new;
                else
                    aux5temp{ifile} =[itrigger(istim), pretime];
                end
                
                % bloc.aux5{ifile}=[itrigger(istim), indstim(istim)];
                
                NIRS.Dt.fir.pp(lst+1).p{ifile,1} = outfile;                
                      
                
                % NIRS.Dt.fir.aux5{ifile,1} = NIRS.Dt.fir.aux5{ifile,1}
            
                fprintf('%s\n',outfile);
                %             %write new .vmrk file
                %             infilewmrk = fullfile(dir1,[fil1 '.vmrk']);
                %             outfilevmrk = fullfile(dir1,[prefix fil1 '.vmrk']);
                %             copyfile(infilewmrk,outfilevmrk);
                
                ChannelLabels = ConvertmlIDsrs2label(NIRS)
                nirs_boxy_write_vhdr(fileOutRoot_vhdr,... %Output file
                    outfile,... %DataFile
                    fileOutRoot_vmrk,... %MarkerFile,...
                    'nirs_run_segment',... %Function that created the header
                    '',... %Channel Resolution
                    '',... %Channel Units
                    ChannelLabels,... %names given as a column of cells
                    1/NIRS.Cf.dev.fs*1000000,... %SamplingInterval in microseconds
                    nstim); %SamplingInterval in microseconds
                
                SD.Markers{1,1}.Type='New Segment';
                SD.Markers{1,1}.Description='';
                SD.Markers{1,1}.Position=1;
                SD.Markers{1,1}.Size=1;
                SD.Markers{1,1}.ChNumber=0; %all channels
                temp_markers{1,1} = SD.Markers{1,1};
                nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
                    outfile,... %DataFile
                    temp_markers); 
                
                aux5new = aux5temp{ifile};                
                [label,ind_dur_ch] = read_vmrk_all(fileOutRoot_vmrk);
                for itrig = 1:size(aux5new,1)
                     label_2 = [{'trigger'}, {['S ',num2str(aux5new(itrig,1))]}];
                     label = [label;label_2];
                     ind_dur_ch_2 = [aux5new(itrig,2),1,0];
                     ind_dur_ch = [ind_dur_ch;ind_dur_ch_2];
                end
                clear aux5new
                ind_dur_ch_2 = mat2d2ind_dur_ch(double(dnoise));
                label_2 = cell(size(ind_dur_ch_2,1),2);
                label_2(:,1)={'bad_step'};
                label_2(:,2)={'manual'};
                ind_dur_ch = [ind_dur_ch;ind_dur_ch_2];
                label = [label;label_2];
                write_vmrk_all(fileOutRoot_vmrk,ind_dur_ch,label);
                ifile = ifile +1;
            end      
         end
      end

    end
    NIRS.Dt.fir.pp(lst+1).pre = 'Segmentation';
    NIRS.Dt.fir.pp(lst+1).job = job;
       for f=1:numel(rDtp)
            if DelPreviousData
                delete(rDtp{f,1});
                [dir1,fil1,ext] = fileparts(rDtp{f,1})
                delete(fullfile(dir1,[fil1,'.vmrk']));
                delete(fullfile(dir1,[fil1,'.vhdr']));
                try
                    infileAC = fullfile(dir1,[fil1 'AC' '.nir']);
                    delete(infileAC)
                    infilePH = fullfile(dir1,[fil1 'PH' '.nir']);
                    delete(infilePH)
                catch
                end
            end
       end
        %add outfile name to NIRS
        if 1%
            NIRS.Cf.H.C.ok = chok;
            NIRS.Dt.fir.aux5ini = auxraw; %trig present dans le fichier raw
            NIRS.Dt.fir.aux5 = aux5temp;
            NIRS.Dt.fir.pp(lst).p = rDtptmp;
        end
        save(fullfile(dir2,'NIRS.mat'),'NIRS');
        job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');
end
    out.NIRSmat = job.NIRSmat;