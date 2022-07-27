function out = nirs_run_normalize(job)
% Normalize the raw data using the time and duration of significant subintervals.
% The normalized data is recorded in a new .nir binary file.
%

%filename prefix
prefix = 'n'; %for "normalize"
DelPreviousData  = job.DelPreviousData;

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
    NIRS = [];
    load(job.NIRSmat{filenb,1});
    [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});
    %use last step of preprocessing
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    fprintf('%s\n','File processed convert raw intensity in delta optical density dOD=log(I/Io)');
    if isfield(job.normtype,'b_choicenormstim')
        fprintf('%s\n', ['Define piecewise normalization around trigger: ',num2str(job.normtype.b_choicenormstim.trigger),', interval I include from PreTime= ', job.normtype.b_choicenormstim.pretime ,' sec to PostTime= ', job.normtype.b_choicenormstim.posttime,' sec'])
        if job.normtype.b_choicenormstim.m_choiceNan
             fprintf('%s\n', ['Exclude artifacts in the Io definition']);
        else
            fprintf('%s\n', ['Keep artifacts in the Io definition']);
        end
        if job.normtype.b_choicenormstim.m_NormType==0
            fprintf('%s\n','Io define from trigger-Pretime to trigger onset'); 
        elseif job.normtype.b_choicenormstim.m_NormType==1
             fprintf('%s\n','Io define from trigger-Pretime to trigger+Posttime onset');             
        end
    end
    if  isfield(job.normtype,'b_choiceglobal')
        fprintf('%s\n', ['Define Io as the whole file average'])
    end
    if isfield(job.normtype,'b_choiceinternan')
        fprintf('%s\n', 'Define Io as the clean period between 2 artifact identify in yellow')   
    end
    ifile = 1; %Utiliser si on remets chaque stim normaliser dans des fichier différent
    for f=1:numel(rDtp) %Loop over all files of a NIRS.mat 
        d = fopen_NIR(rDtp{f},NC);
        dnorm = zeros(size(d));
        % Check negative value (could be induce by the offset ajuste) 
        % Add offset to avoid negative log ! 
       [val,ind] =min(d');
       idchneg = find(val<0);
     %  figure;plot( d(idchneg,:)')
       
       if ~isempty(idchneg)        
           valadjust = val(idchneg)';      
           d(idchneg,:) =  d(idchneg,:) - valadjust*ones(1,size(d,2)) + 0.2; %(not zero add 0.2 usual smaller value NIRx)
       end

    % figure;imagesc(dnan)
    if isfield(job.normtype,'b_choiceglobal'); %Global normalization       
            if job.normtype.b_choiceglobal.m_choiceNan
                [dir1,fil1,ext]=fileparts(rDtp{f});
                vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
                [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
                nch = size(d,1);
                nsample = size(d,2);
                noise = ind_dur_ch2mat(ind_dur_ch, nsample,nch)';
                idnoise = find(double(noise));
                dnan = d;
                if ~isempty(idnoise)
                    dnan(idnoise)=nan;
                end
                %figure;imagesc(dnan);
                for Idx = 1:NC
                    meansub = nanmean(dnan(Idx,:)); %Normalize whole blocks 
                    dnorm(Idx,:) = log10(d(Idx,:)./meansub); %dOD
                end                
            else
                for Idx = 1:NC
                    meansub = nanmean(d(Idx,:)); %Normalize whole blocks 
                    dnorm(Idx,:) = log10(d(Idx,:)./meansub); %dOD
                end
            end           
        elseif isfield(job.normtype,'b_choiceinternan') %Inter NAN normalization
             [dir1,fil1,ext]=fileparts(rDtp{f});
            vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
            [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
            nch = size(d,1);
            nsample = size(d,2);
            noise = ind_dur_ch2mat(ind_dur_ch, nsample,nch)';
            idnoise = find(double(noise));
            dnan = d;
            if ~isempty(idnoise) 
                dnan(idnoise)=nan;
            end                       
            for Idx = 1:NC
                segnanid = isnan(dnan(Idx,1:end-1))-isnan(dnan(Idx,2:end));
                idstartint = find( segnanid==+1);
                idstopint = find( segnanid==-1);          
                if numel(idstartint) == numel( idstopint) %start stop 
                    segnan = [[1,find( segnanid==+1) ];[find( segnanid==-1),size(d,2)]];
                elseif numel(idstartint) < numel(idstopint)
                    numel(segnanid);;
                    segnan = [[1,idstartint ];[idstopint(1:end-1),size(d,2)]];
                elseif numel(idstartint) > numel(idstopint)
                    numel(segnanid);
                    segnan = [[1,idstartint(1:end-1) ];[idstopint,size(d,2)]];
                end
                for iseg = 1:size(segnan,2)
                    meansub = nanmean(dnan(Idx,segnan(1,iseg):segnan(2,iseg))); %Normalize using mean between nan periode
                    dnorm(Idx,segnan(1,iseg):segnan(2,iseg)) = log10(d(Idx,segnan(1,iseg):segnan(2,iseg))./meansub);
                end
            end
        elseif isfield(job.normtype,'b_choicenormstim') %Stim normalization
            trigger = job.normtype.b_choicenormstim.trigger;
            pretime = round(fs*str2num(job.normtype.b_choicenormstim.pretime))-1;
            posttime = round(fs*str2num(job.normtype.b_choicenormstim.posttime));
            NIRS.Dt.fir.aux5bloc = NIRS.Dt.fir.aux5{f};
            aux5 = NIRS.Dt.fir.aux5{f};
            if trigger == 0 %If 0, consider every trigger 
                indstim = aux5(:,2);
                ind = [indstim; size(d,2)];
                itrigger = 1;
            else
                indstim = [];
                itrigger = [];
                for itypestim = 1:numel(trigger)
                    idstim= aux5((aux5(:,1) == trigger(itypestim)),2);
                    indstim  = [indstim  ;idstim ];
                    ind = [1; idstim; size(d,2)];
                    itrigger =[itrigger, ones(1,numel(idstim)).*trigger(itypestim)];
                end
                
            end
            if job.normtype.b_choicenormstim.m_choiceNan
                [dir1,fil1,ext]=fileparts(rDtp{f});
                vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
                [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
                nch = size(d,1);
                nsample = size(d,2);
                noise = ind_dur_ch2mat(ind_dur_ch, nsample,nch)';
                idnoise = find(double(noise));
                dnan = d;
                if ~isempty(idnoise)
                    dnan(idnoise)=nan;
                end  
            else
                dnan = d;
            end 
            
            for istim = 1:numel(indstim)                
                
                if find((indstim(istim)-pretime) <= 0)
                    disp(['Pretime to short '])
                    istart = 1;
                else
                    istart = indstim(istim)-pretime;
                end
               if find((indstim(istim)+posttime) > size(d,2))
                    disp(['Prepostime to long '])
                    istop = size(d,2);
               else
                   istop = indstim(istim)+posttime;
                end
                               
%                 if (max(indstim)+posttime)>size(d,2)
%                     disp('Warning out of range padding in the last bloc')
%                     d = [d,fliplr(d)];
%                 end
                %figure;imagesc(dnan)
                if job.normtype.b_choicenormstim.m_NormType==0      %I/Io Io = pretime
                   for Idx = 1:NC
                        meansub= nanmean(dnan(Idx,istart:indstim(istim)),2);
                        dnorm(Idx, istart:istop) = log10(d(Idx,istart:istop)./meansub);
                   end
                elseif job.normtype.b_choicenormstim.m_NormType==1  %I/Io Io = pretime to posttime
                    for Idx = 1:NC
                        meansub= nanmean(dnan(Idx,istart:istop),2);
                        dnorm(Idx, istart:istop) = log10(d(Idx,istart:istop)./meansub);
                    end      
                end               
            end
         end
        
        [dir1,fil1,ext1] = fileparts(rDtp{f});
        infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
        infilevhdr = fullfile(dir1,[fil1 '.vhdr']);
        if ~exist(dir2,'dir'), mkdir(dir2); end;
        outfile = fullfile(dir2,[prefix fil1 ext1]);
        outfilevmrk = fullfile(dir2,[prefix fil1  '.vmrk']);
        outfilevhdr = fullfile(dir2,[prefix fil1  '.vhdr']);
        fwrite_NIR(outfile,dnorm);
        copyfile(infilevmrk,outfilevmrk);

        
      try 
                 try
                    info = read_vhdr_brainvision((fullfile(dir1,[fil1,'.vhdr'])));
                    ChannelLabels = info.label;
                catch
                    ChannelLabels = ConvertmlIDsrs2label(NIRS);
                end
                    SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
                    nirs_boxy_write_vhdr(outfilevhdr,... %Output file
                        outfile,... %DataFile
                        outfilevmrk,... %MarkerFile,...
                        'nirs_run_normalize',... %Function that created the header
                        '',... %Channel Resolution
                        '',... %Channel Units
                        ChannelLabels,... %names given as a column of cells
                        SamplingInterval,...
                        size(d,2)); %SamplingInterval in microseconds
      catch
      end
      
        NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
        fprintf('%s\n',outfile);
        
        
        if DelPreviousData
            delete(rDtp{f,1});
            delete(infilevmrk);
            delete(infilevhdr);
            try
                infileAC = fullfile(dir1,[fil1 'AC' '.nir']);
                delete(infileAC)
                infilePH = fullfile(dir1,[fil1 'PH' '.nir']);
                delete(infilePH)
            catch
            end
            disp(['Delete previous .nir data file: ',rDtp{f,1}]);
        end
        NIRS.Dt.fir.pp(lst+1).pre = 'Normalization';
        NIRS.Dt.fir.pp(lst+1).job = job;
    end
    
    save(fullfile(dir2,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');
    
end


out.NIRSmat = job.NIRSmat;


