function out = nirs_run_segment(job)
%Use to segment NIRS event around trigger
%Synchronised according to mutual trigger information EEG, AUX and Video Field
%
%filename prefix 
prefix = 's'; %for "segment"
DelPreviousData  = job.DelPreviousData;

 %write offset time in the AUX and EEG structure and do not create a copy of the data segmented 0 keep whole original segment 1 segment
flagsegmentAUX = 0;
flagsegmentEEG = 0;
for filenb = 1:size(job.NIRSmat,1)
    NIRS = [];
    load(job.NIRSmat{filenb,1});
    [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});
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
        [ind_dur_ch1] = read_vmrk_find(vmrk_path,'bad_step');
        [ind_dur_ch2] = read_vmrk_find(vmrk_path,'Bad Interval');
        ind_dur_ch = [ind_dur_ch1 ind_dur_ch2];
        nsample = size(d,2);
         noise = ind_dur_ch2mat(ind_dur_ch, nsample,NC)';
  
        if isfield(NIRS.Dt,'EEG')
            try
                fileeeg=NIRS.Dt.EEG.pp(moduleEEG).p{f}; %OPEN EEG HERE
                 EEG.file = fileeeg;
                [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]=fopen_EEG(fileeeg);
            catch
            end
            
          if 0
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
                AUX(iaux).file = fileaux;
                [AUX(iaux).data,AUX(iaux).infoBV,AUX(iaux).marker,AUX(iaux).ind_dur_ch]=fopen_EEG(fileaux);
            end
        end
        catch
           % msgbox('AUX could not be segment')
           disp('Error AUX could not be segment')
        end
        
        trigger = job.trigger;

 
    
      
        
       
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
                idstim= aux5((aux5(:,1) == trigger(itypestim)),2);
                indstim  = [indstim  ;idstim ];
                ind = [1; idstim; size(d,2)];
                itrigger =[itrigger, ones(1,numel(idstim)).*trigger(itypestim)];
                disp(['Search trigger: ', num2str(trigger(itypestim))])  
            end
            disp(['Find ', num2str(numel(indstim)), ' NIRS trigger ', sprintf('S%3.0f, ',trigger),sprintf('\n'),...
                'Time: ', sprintf('%.2f ',1/fs*indstim),'seconds to segment',sprintf('\n'),...
                'Sample: ',sprintf('%.0f ', indstim )]);
                
            %disp(['Find ', num2str(numel(indstim_EEG)), ' EEG trigger ',idtrigEEG , sprintf('\n'),'Time: ', sprintf('%.2f, ',EEG.infoBV.SamplingInterval/1000000*indstim_EEG),'seconds to sync', sprintf('\n'),'Sample: ',sprintf('%d, ',indstim_EEG) ])
             try
            if numel(indstim)==1
                %trig could be ajuste to start and end of the bloc no need
                %to define fix start stop cause only one trigger exist. 
            else
                if job.m_SegmentTrig==0 % use all trig
                    switch job.pretime
                      case 'start'
                          disp('Multiple trigs are not recommend only one trig with PreTime <Start>  consider to change the menu Option trig: <Use all> for <Only first>.');
                     end
                     switch job.posttime
                      case 'end'
                          disp('Multiple trigs are not recommend only one trig with PostTime <End>  consider to change the menu Option trig: <Use all> for <Only first>.');
                     end
                elseif job.m_SegmentTrig==1 % use first trig
                    indstim = indstim(1);
                   disp(['Use only first trig ', num2str(trigger(1)) ,' to segment at Time: ', sprintf('%.2f ',1/fs*indstim) ])
                end
                 
            end
             catch
                 disp('Failed NIRS segment trig could not be found to segment please verify your data or your parameter')
                 out.NIRSmat = job.NIRSmat;
                 return
             end
             
              if strfind(job.pretime,'start')
                  
                       if numel(indstim)==1
                            pretime = indstim - 1 ;
                       else
                          disp('Warning using pretime ''start'', only one trigger are recommand for this special option, multi-segment imply equal pretime.')
                          pretime =indstim(1) - 1 ;
                       end
              elseif strfind(job.pretime,'Comments')
                    job.pretime(9:end)
                    if ~isfield(NIRS.Dt.fir,'comments')
                        disp('Error, fields comments do not exist no segmentation')
                        return
                    end
                 tmp=   NIRS.Dt.fir.comments{1}
           
                  findevent = strfind( tmp(:,1),strtrim(job.pretime(9:end)))
                  for idevt=1:numel(findevent)
                      if ~isempty(findevent{idevt})
                            pretime = tmp{idevt,2}
                      end
                  end
              else %specific time 
                      pretime  = round(fs*str2num(job.pretime));
                      
              end
              if strcmp(job.posttime,'end')                 
                      if numel(indstim)==1 
                           posttime = size(d,2)- indstim;
                      else
                          disp(['Verify trigger present in file : '  vmrk_path]);
                          disp('Warning using posttime ''end '', only one trigger are recommand for this special option, multi-segment imply equal posttime.');
                          posttime = size(d,2)- indstim(1);
                      end
              elseif strfind(job.posttime, 'trig')
                  triggerend = str2num(job.posttime(5:end));
                   idstimend= aux5(aux5(:,1) == triggerend,2)+1;
                  posttime = (idstimend - indstim(1));
                   disp(['Using posttime at ' job.posttime, ' at time: ', num2str(idstimend/fs) , 's find a duration segment of ', num2str(posttime),'sample']);
                 
              else  %specific time number with time duration in second
                     posttime = round(fs*str2num(job.posttime));
              end
           
            if isfield(NIRS.Dt,'EEG')
                try
                    tmp = EEG.marker{end,2};
                    idstimEEG = [];
                    tmp = EEG.marker{end,2};
                    if iscell(tmp)
                        tmp = tmp{1};
                    end
                    if isempty(tmp)
                        tmp = {'nul'};
                    end
                    for itypestim = 1:numel(trigger)
                        if strcmp(tmp(1),'S')|strcmp(NIRS.Cf.dev.n,'NIRx') %S suffix nirx but not ISS
                            idtrigEEG = ['S',sprintf('%3.0f',trigger(itypestim))];
                        else
                            idtrigEEG = [sprintf('%3.0f',trigger(itypestim))];
                        end
                        
                        idstimEEG = [idstimEEG;strmatch(idtrigEEG,EEG.marker(:,2))];
                    end
 
                if job.m_SegmentTrig==1
               
                else
                      if isfield(NIRS.Dt.EEG.pp(moduleEEG),'sync_timesec')
                        tstart= NIRS.Dt.EEG.pp(moduleEEG).sync_timesec{f};
                        tstop = tstart+size(d,2)*1/NIRS.Cf.dev.fs     ;               
                        valtimestim =       EEG.ind_dur_ch(idstimEEG,1).*EEG.infoBV.SamplingInterval/1000000;
                        triginside = find(valtimestim>tstart & valtimestim<tstop)      ;             
                         tmp = idstimEEG ;
                         idstimEEG = tmp(triginside);
                      end
                end                         
                 
                         disp(['Find ', num2str(numel(idstimEEG)), ' AUX file ',EEG.file, ' trigger ', sprintf('S%3.0f, ',trigger),sprintf('\n'),...
                            'Time: ', sprintf('%.2f ',EEG.infoBV.SamplingInterval/1000000.* EEG.ind_dur_ch(idstimEEG)),'seconds to segment',sprintf('\n'),...
                            'Sample: ',sprintf('%.0f ', EEG.ind_dur_ch(idstimEEG,1)) ])
                        
                        
                    idtrigEEG = ['S',sprintf('%3.0f',trigger)];
                    if numel(indstim)==numel(idstimEEG) %gerer les cas stim agree dimention only
                        indstim_EEG = EEG.ind_dur_ch(idstimEEG,1);
                        disp(['Find ', num2str(numel(indstim_EEG)), ' EEG trigger ',idtrigEEG , sprintf('\n'),'Time: ', sprintf('%.2f, ',EEG.infoBV.SamplingInterval/1000000*indstim_EEG),'seconds to sync', sprintf('\n'),'Sample: ',sprintf('%d, ',indstim_EEG) ])
                    else
                        %msgbox('unequal EEG trigger identification')
                        if job.m_SegmentTrig==1 % use first trig
                            indstim_EEG = EEG.ind_dur_ch(idstimEEG(1),1);    
                            disp(['Use only first trig ', num2str(numel(indstim_EEG)), ' EEG trigger ',idtrigEEG , sprintf('\n'),'Time: ', sprintf('%.2f, ',EEG.infoBV.SamplingInterval/1000000*indstim_EEG),'seconds to sync', sprintf('\n'),'Sample: ',sprintf('%d, ',indstim_EEG) ])
                        else
                            indstim_EEG = EEG.ind_dur_ch(idstimEEG,1);
                            disp('Error unequal EEG trigger identification review vmrk')
                            disp(['Find ', num2str(numel(indstim_EEG)), ' EEG trigger ',idtrigEEG , sprintf('\n'),'Time: ', sprintf('%.2f, ',EEG.infoBV.SamplingInterval/1000000*indstim_EEG),'seconds to sync', sprintf('\n'),'Sample: ',sprintf('%d, ',indstim_EEG) ])  
                        end
                    end
                    fsEEG = 1/(EEG.infoBV.SamplingInterval/1000000); %uS en seconde
                    pretimeEEG = round(fsEEG*pretime/fs)-1;
                    posttimeEEG = round(fsEEG*posttime/fs);
                catch 
                    disp('Error EEG could not be open and synchronised')
                end
            end
           
            if isfield(NIRS.Dt,'AUX')
                    try
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
                        if strcmp(tmp(1),'S')|strcmp(NIRS.Cf.dev.n,'NIRx') %S suffix nirx but not ISS
                            idtrigAUX{iaux} = ['S',sprintf('%3.0f',trigger(itypestim))];
                        else
                            idtrigAUX{iaux} = [sprintf('%3.0f',trigger(itypestim))];
                        end
                        
                        idstimAUX{iaux} = [idstimAUX{iaux};strmatch(deblank(idtrigAUX{iaux}),deblank(AUX(iaux).marker(:,2)))];
                    end
                    %PRESEGMENTATION ALLREALY DONE.
                      if job.m_SegmentTrig==1
                          
                      else
                        if isfield(NIRS.Dt.AUX(iaux).pp(moduleaux),'sync_timesec')
                        tstart= NIRS.Dt.AUX(iaux).pp(moduleaux).sync_timesec{f};
                        tstop = tstart+size(d,2)*1/NIRS.Cf.dev.fs*AUX(iaux).infoBV.SamplingInterval/1000000;
                        indstim*1/NIRS.Cf.dev.fs;
                        valtimestim =        AUX(iaux).ind_dur_ch(idstimAUX{iaux},1).*AUX(iaux).infoBV.SamplingInterval/1000000;
                        triginside = find(valtimestim>tstart & valtimestim<tstop)        ;           
                        tmp = idstimAUX{iaux} ;
                        idstimAUX{iaux} = tmp(triginside);
                        end
                      end
                           disp(['Find ', num2str(numel(idstimAUX{iaux})), ' AUX file ',AUX(iaux).file, ' trigger ', sprintf('S%3.0f, ',trigger),sprintf('\n'),...
                            'Time: ', sprintf('%.2f ',AUX(iaux).infoBV.SamplingInterval/1000000.* AUX(iaux).ind_dur_ch(idstimAUX{iaux})),'seconds to segment',sprintf('\n'),...
                            'Sample: ',sprintf('%.0f ', AUX(iaux).ind_dur_ch(idstimAUX{iaux},1)) ])
                     %itrigger =[itrigger, ones(1,numel(idstim)).*trigger(itypestim)];
                    if numel(indstim)==numel(idstimAUX{iaux}) %gerer les cas stim agree dimension only
                        indstim_AUX{iaux} = AUX(iaux).ind_dur_ch(idstimAUX{iaux},1);
                         % disp(['Find AUX onset time: ', sprintf('%.2f ',AUX.infoBV.SamplingInterval/1000000*indstim_AUX{iaux}),'seconds to sync' ])
                        disp(['Find ', num2str(numel( indstim_AUX{iaux})), ' AUX trigger ', sprintf('S%3.0f, ',trigger) ,...
                            sprintf('\n'),'Time: ', sprintf('%.2f, ',AUX(iaux).infoBV.SamplingInterval/1000000*indstim_AUX{iaux}),'seconds to sync', sprintf('\n'),'Sample: ',sprintf('%d, ',indstim_AUX{iaux}) ]);

                    else
%                         idtocheckAUX = AUX(iaux).ind_dur_ch(idstimAUX{iaux},1);
%                         timetrig = idtocheckAUX * AUX(iaux).infoBV.SamplingInterval/1000000;
%                         timestim = indstim*1/NIRS.Cf.dev.fs;
%                         if 0
%                         figure;hold on;
%                         for i=1:numel(timetrig)
%                             plot( timetrig(i),1,'+r','displayname',['AUX',num2str(i)])
%                         end                 
%                         for i=1:numel(timestim)
%                             plot(timestim(i),1,'+b','displayname',['NIRS',num2str(i)])
%                         end
%                         title(['Error in number of trig ',num2str(trigger) ,' between AUX ', NIRS.Dt.AUX(iaux).label,' and NIRS, Please Ajust trig manualy or in VMRK file to make them match before segmentation.']);
%                         end
                        if job.m_SegmentTrig==1 % use first trig
                            temp = idstimAUX{iaux};
                            indstim_AUX{iaux}= AUX(iaux).ind_dur_ch(temp(1),1);
                            disp(['Use only first trig ', num2str(numel( indstim_AUX{iaux})), ' AUX trigger ', sprintf('S%3.0f, ',trigger) ,...
                            sprintf('\n'),'Time: ', sprintf('%.2f, ',AUX(iaux).infoBV.SamplingInterval/1000000*indstim_AUX{iaux}),'seconds to sync', sprintf('\n'),'Sample: ',sprintf('%d, ',indstim_AUX{iaux}) ]);
                        else 
                            indstim_AUX{iaux} = AUX(iaux).ind_dur_ch(idstimAUX{iaux},1);
                            disp('Error unequal AUX trigger identification review vmrk')
                            disp(['Find ', num2str(numel( indstim_AUX{iaux})), ' AUX trigger ', sprintf('S%3.0f, ',trigger) ,...
                            sprintf('\n'),'Time: ', sprintf('%.2f, ',AUX(iaux).infoBV.SamplingInterval/1000000*indstim_AUX{iaux}),'seconds to sync', sprintf('\n'),'Sample: ',sprintf('%d, ',indstim_AUX{iaux}) ]);
                        end
                    end
                    
                    fsAUX{iaux} = 1/(AUX(iaux).infoBV.SamplingInterval/1000000); %uS en seconde
                    pretimeAUX{iaux} = round(fsAUX{iaux}*pretime/fs)-1;
                    posttimeAUX{iaux} = round(fsAUX{iaux}*posttime/fs);
                    
                end
                catch
                    disp('Error AUX could not be open and synchronised')
                end
            end   
        end
        disp(['PreTime: ',num2str(pretime*1/fs), ' s ,PostTime: ' ,num2str(posttime*1/fs), ' s Apply to create ', num2str(numel(indstim)),' segments at raw time: ', num2str(indstim'*1/fs) ]);
        nstim = posttime + pretime + 1;
         
        t = 1/fs:1/fs:1/fs*size(d,2);
        if ~isempty(indstim)
        for istim = 1:numel(indstim)
             dnorm = zeros(NC,nstim);
            if find((indstim(istim)-pretime) <= 0)
                 disp('Warning could not be process : Segment pretime before stim out of range use 0 or negative value to take from the trig or la ')
                 return
                % indstim(:) = indstim(:) + size(d,2) ;
                % d =[fliplr(d),d];
            end
           
            if (max(indstim)+posttime)>size(d,2)              
               % msgbox('Warning out of range padding in the last bloc')
                disp('Warning out of range mirror padding for fNIRS data in the last bloc to apply PostTime lenght')
                d = [d,fliplr(d)]; 
                noise = [noise,fliplr(noise)];                
            end
            %Segment the data NIRS
            dnorm = d(:,indstim(istim)-pretime:indstim(istim)+posttime);
            dnoise = noise(:,indstim(istim)-pretime:indstim(istim)+posttime);
            %figure;imagesc(dnoise)
            
            
            
            bloc = sprintf('%s%02.0f','b',istim);
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            if 1
                outfile = fullfile(dir2,[prefix fil1 bloc ext1]);
                fileOutRoot_vmrk= fullfile(dir2,[prefix fil1 bloc '.vmrk']);
                fileOutRoot_vhdr = fullfile(dir2,[prefix fil1 bloc '.vhdr']) ;
                if isfield(NIRS.Dt,'EEG')
                    try
                    %WRITE IN SUB FILE.... or add sync_timesec
                    [dirEEG,filEEG,extEEG]= fileparts(NIRS.Dt.EEG.pp(moduleEEG).p{f});                    
                    outfileEEG = fullfile(dirEEG,[filEEG bloc '.dat']);
                    stimlist=indstim_EEG;                      
                    if flagsegmentEEG %not used create copy of EEG file segmented
                        %write and segment EEG file in homologous bloc
                        if (max(stimlist)+posttimeEEG)>size(EEG.data,1)
                            %msgbox('Warning out of range padding in the last bloc aux')
                            disp('Warning out of range mirror padding for EEG data in the last bloc to apply PostTime lenght')
                            EEG.data= [EEG.data;flipud(EEG.data)];
                            %  figure;plot(AUX(iAUX).data)
                        end
                        datasegR = EEG.data(indstim_EEG(istim)-pretimeEEG:indstim_EEG(istim)+posttimeEEG,:);
                        fwrite_EEG(outfileEEG,EEG,indstim_EEG(istim)-pretimeEEG,indstim_EEG(istim)+posttimeEEG);
                        bloc                       
                        NIRS.Dt.EEG.pp(moduleEEG+1).p{ifile,1}=outfileEEG;                             
                    else  %apply reference to file and indicate offset time
                        NIRS.Dt.EEG.pp(moduleEEG+1).p{ifile} = NIRS.Dt.EEG.pp(moduleEEG).p{f};
                        NIRS.Dt.EEG.pp(moduleEEG+1).sync_timesec{ifile} = (indstim_EEG(istim)-pretimeEEG)*1/fsEEG;  
                        disp(['EEG sync time ', num2str(NIRS.Dt.EEG.pp(moduleEEG+1).sync_timesec{ifile} ) ' sec'])
                    end
                    catch
                        disp('Error Segment EEG could not be open and synchronised')
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
                                disp('Error Video segmentation could not be done using EEG trigger')                                
                            end
                     
                        case 'NIRS'
                            try 
                            NIRS.Dt.Video.pp(moduleVideo+1).p{ifile} = NIRS.Dt.Video.pp(moduleVideo).p{f};
                            NIRS.Dt.Video.pp(moduleVideo+1).sync_timesec{ifile} = (indstim(istim)-pretime)*1/fs  + offset;
                            catch
                                disp('Error Video segmentation could not be done using NIRS trigger') 
                            end
                         case 'AUX'
                             try
                                NIRS.Dt.Video.pp(moduleVideo+1).p{ifile} = NIRS.Dt.Video.pp(moduleVideo).p{f};
                                indstim_AUXval = indstim_AUX{1}; %synchronisation will refere to file 1 
                                pretimeAUXval = pretimeAUX{1};
                                fsAUXval = fsAUX{1};                            
                                NIRS.Dt.Video.pp(moduleVideo+1).sync_timesec{ifile} = (indstim_AUXval(istim)- pretimeAUXval)*1/fsAUXval  + offset;
                             catch
                                 disp('Error Video segmentation could not be done using AUX trigger') 
                             end
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
                                %msgbox('Audio segmentation could not be done using EEG trigger')
                                disp('Error Audio segmentation could not be done using EEG trigger')
                            end
                     
                        case 'NIRS'
                            try
                            NIRS.Dt.Audio.pp(moduleAudio+1).p{ifile} = NIRS.Dt.Audio.pp(moduleAudio).p{f};
                            NIRS.Dt.Audio.pp(moduleAudio+1).sync_timesec{ifile} = (indstim(istim)-pretime)*1/fs  + offset;
                            catch
                                disp('Error Audio segmentation could not be done using fNIRS trigger')
                            end
                        case 'AUX'
                            try
                            NIRS.Dt.Audio.pp(moduleAudio+1).p{ifile} = NIRS.Dt.Audio.pp(moduleAudio).p{f};
                            indstim_AUXval = indstim_AUX{1}; %synchronisation will refere to file 1 
                            pretimeAUXval = pretimeAUX{1};
                            fsAUXval = fsAUX{1};                            
                            NIRS.Dt.Audio.pp(modulemoduleAudio+1).sync_timesec{ifile} = (indstim_AUXval(istim)- pretimeAUXval)*1/fsAUXval  + offset;
                            catch
                                disp('Error Audio segmentation could not be done using AUX trigger')
                            end
                    end
                    
                end
                
              
                if isfield(NIRS.Dt,'AUX')
                    try
                    for iAUX = 1:numel(NIRS.Dt.AUX)
                        
                        [dirAUX,filAUX,extAUX]= fileparts(NIRS.Dt.AUX(iAUX).pp(moduleaux).p{f}); %AUX(f).file             
                        outfileAUX{iAUX} = fullfile(dirAUX,[filAUX bloc '.dat']);
                        %write and segment EEG file in homologous bloc
                        stimlist=indstim_AUX{iAUX};
                      
                        if flagsegmentAUX
                            if (max(stimlist)+posttimeAUX{iAUX})>size(AUX(iAUX).data,1)
                            %msgbox('Warning out of range padding in the last bloc aux')
                            disp('Warning aux data out of range regard to NIRS size, mirror padding in the last bloc aux')
                            AUX(iAUX).data= [AUX(iAUX).data;flipud(AUX(iAUX).data)];
                            %  figure;plot(AUX(iAUX).data)
                            end
                            datasegR = AUX(iAUX).data(stimlist(istim)-pretimeAUX{iAUX}:stimlist(istim)+posttimeAUX{iAUX},:);
                            fwrite_EEG(outfileAUX{iAUX},AUX(iAUX),stimlist(istim)-pretimeAUX{iAUX},stimlist(istim)+posttimeAUX{iAUX});
                            NIRS.Dt.AUX(iAUX).pp(moduleaux+1).p{ifile,1}=outfileAUX{iAUX};
                        else
                            NIRS.Dt.AUX(iAUX).pp(moduleaux+1).p{ifile,1}=NIRS.Dt.AUX(iAUX).pp(moduleaux).p{f}; 
                            NIRS.Dt.AUX(iAUX).pp(moduleaux+1).sync_timesec{ifile,1} = (stimlist(istim)-pretimeAUX{iAUX})*1/fsAUX{iAUX};
                             disp(['AUX sync time ', num2str(NIRS.Dt.AUX(iAUX).pp(moduleaux+1).sync_timesec{ifile,1}) ' sec'])
                        end
                        bloc;
                    end
                    
               
                catch
                       % msgbox('AUX could not be segment')
                       disp('Error AUX could not be segmented');
                end
                end
                %dnorm=log(dnorm)
                %WRITE THE RESEGMENTATION
                try
                    NIRS.Dt.fir.sizebloc{ifile}=size(dnorm,2);
                catch
                    NIRS.Dt.fir.sizebloc(ifile)=size(dnorm,2);
                end
                fwrite_NIR(outfile,dnorm);
                
                
                
                rDtptmp{ifile} = rDtp{f};   %Nom temporaire pour gere la division du fichier et garder la trace du fichier original
                auxraw{ifile} = NIRS.Dt.fir.aux5{f};
                chok(:,ifile) = NIRS.Cf.H.C.ok(:, f);
                %reajuster tout les trigs between pretime and postime of the
                %selected segment
                aux5all = NIRS.Dt.fir.aux5{f};
                idtrigtokeep = find(aux5all(:,2)>(indstim(istim) - pretime) & aux5all(:,2)<(indstim(istim) +posttime));
                if ~isempty( idtrigtokeep)
                    aux5new(:,1) = aux5all(idtrigtokeep,1);
                    aux5new(:,2) = aux5all(idtrigtokeep,2) - indstim(istim) +pretime;
                    aux5temp{ifile} = aux5new;
                    clear aux5new;
                else
                    aux5temp{ifile} =[itrigger(istim), pretime];
                end
                
                % bloc.aux5{ifile}=[itrigger(istim), indstim(istim)];
                
                NIRS.Dt.fir.pp(lst+1).p{ifile,1} = outfile;                
                      
                
                % NIRS.Dt.fir.aux5{ifile,1} = NIRS.Dt.fir.aux5{ifile,1}
            
                fprintf('%s\n',['Create: ', outfile]);
                %             %write new .vmrk file
                %             infilewmrk = fullfile(dir1,[fil1 '.vmrk']);
                %             outfilevmrk = fullfile(dir1,[prefix fil1 '.vmrk']);
                %             copyfile(infilewmrk,outfilevmrk);
                try 
                    info = read_vhdr_brainvision((fullfile(dir1,[fil1,'.vhdr'])));
                    ChannelLabels = info.label;
                catch
                    ChannelLabels = ConvertmlIDsrs2label(NIRS);
                end
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
                [dir1,fil1,ext] = fileparts(rDtp{f,1});
                delete(fullfile(dir1,[fil1,'.vmrk']));
                delete(fullfile(dir1,[fil1,'.vhdr']));
                try
                    infileAC = fullfile(dir1,[fil1 'AC' '.nir']);
                     fid=fopen(infileAC)~=-1;
                    if fid~=-1
                        fclose(fid);
                        delete(infileAC);
                    end
                      infilePH = fullfile(dir1,[fil1 'PH' '.nir']);
                      fid=fopen(infilePH);
                    if fid~=-1
                        fclose(fid);
                        delete(infilePH);
                    end
                catch
                end
                disp(['Delete previous .nir data file: ',rDtp{f,1}])
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

    