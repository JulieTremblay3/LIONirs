%___________________________________________________________________
% Copyright (C) 2023 LION Lab, Centre de recherche CHU Sainte-Justine 
% www.lionlab.umontreal.ca
%___________________________________________________________________
function out = nirs_run_E_prewhitening(job)

prefix = 'w'; %for "whitening"
DelPreviousData  = job.DelPreviousData;

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects

    NIRS = [];
    load(job.NIRSmat{filenb,1});

    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    
    [dir2,~,~] = fileparts(job.NIRSmat{filenb,1});
    
    for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat 
        d = fopen_NIR(rDtp{f,1},NC);
        
        tic
        % CAN CHANGE THESE VALUES -- These are defaults based on our data.
        % Could revise so these are user-defined inputs
        maxMdlOrder = 8;
        pmax = round(maxMdlOrder*fs);
        pctThresh = 0.05;
        
        BIC = nan(pmax,NC);
        Ps = nan(1,NC);
        dw = nan(size(d));
        for i = 1:NC
            y = d(i,:)';
            fprintf('Getting best BIC channel %d of %d\n',i,NC)
            for p = 1:pmax
                mdl = regARIMA(p,0,0);
                
                [~,~,LogL] = estimate(mdl,y,'Display','off');
                [~,BIC(p,i)] = aicbic(LogL,i+1,length(y)-i);
               
                if p>2
                    pctDiff = 100*(diff(BIC(1:p,i))./BIC(1:p-1,i));
                    if pctDiff(end)<pctThresh
                        Ps(i) = p;
                        fprintf('Model order = %d\n',p)
                        break
                    end
                end
            end

            fprintf('Getting innovations channel %d of %d\n',i,NC)
            y = d(i,:)';
            mdl = regARIMA(Ps(i),0,0);
            estMdl = estimate(mdl,y,'Display','off');
            dw(i,:) = infer(estMdl,y);
        end

        % Get bad intervals
        mrks = []; 
        ind = [];
        dur = [];
        mrk_type = 'bad_step';
        mrk_type_arr = cellstr(mrk_type);
        padtime = max(Ps);

        [dir1,fil1,~] = fileparts(rDtp{f});
        vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
        [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);

        % Remove bad intervals from innovations
        if ~isempty(ind_dur_ch)
            for Idx = 1:NC %Loop over all channels
                mrks = find(ind_dur_ch(:,3)==Idx | ind_dur_ch(:,3)==0);
                ind = ind_dur_ch(mrks,1);
                indf = ind + ind_dur_ch(mrks,2);
                for i = 1:numel(ind)
                    if ind(i)-padtime < 1
                        ind(i) = padtime+1;
                    end
                    if indf(i)+padtime > size(d,2)
                        indf(i) = size(d,2)-padtime;
                    end
                    dw(Idx,ind(i):indf(i)+padtime) = NaN;
                end
            end
            nanCheck = zeros(size(dw));
            for i = 1:NC
                nanCheck(i,:) = isnan(dw(i,:));    
            end
            nanCheck = sum(nanCheck)>0;
            dw = dw(:,~nanCheck);
            fprintf('Removed %d bad intervals\n',numel(ind))

            % Detrend for offset adjustment  
            intensnorm = dw';
            X = 1:1:size(intensnorm,1);
            Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
            Mb2 =  intensnorm(1,:)'; %offset    
            A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
            dw = (intensnorm - A)';
            disp('Detrended')
         else
           disp('No bad intervals found');
        end
        [dir1,fil1,ext1] = fileparts(rDtp{f});
        infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
        try
        infilevhdr = fullfile(dir1,[fil1 '.vhdr']);    
        catch
        end
        outfile = fullfile(dir2,[prefix fil1 ext1]);
        outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
        outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);

        % Save .nir file
        fwrite_NIR(outfile,dw); 
        fprintf('%s\n',outfile);

        % Write new .vmrk file 
        try
        copyfile(infilevmrk,outfilevmrk);
        catch
        end
        try
        copyfile(infilevhdr,outfilevhdr);
        try
            info = read_vhdr_brainvision((fullfile(dir1,[fil1,'.vhdr'])));
            ChannelLabels = info.label;
        catch
            ChannelLabels = ConvertmlIDsrs2label(NIRS);
        end
        SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
        nirs_boxy_write_vhdr(outfilevhdr,... %Output file
                    fileOut,... %DataFile
                    outfilevmrk,... %MarkerFile,...
                    'nirs_run_prewhitening',... %Function that created the header
                    '',... %Channel Resolution
                    '',... %Channel Units
                    ChannelLabels,... %names given as a column of cells
                    SamplingInterval,...
                    size(dw,2)); %SamplingInterval in microseconds
        catch
        end

        if DelPreviousData
            try
            delete(rDtp{f,1});                    
            delete(infilevmrk);
            delete(infilevhdr);
            disp(['Delete previous .nir data file: ',rDtp{f,1}]);
            catch                        
            end
        end
        
        %add outfile name to NIRS
        if f == 1
            NIRS.Dt.fir.pp(lst+1).pre = 'Prewhitened';
            NIRS.Dt.fir.pp(lst+1).job = job;
        end 
        NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;

        save(fullfile(dir2,'NIRS.mat'),'NIRS');
        job.NIRSmat{1} = fullfile(dir2,'NIRS.mat');
        
        toc
    end    
end
 out.NIRSmat = job.NIRSmat;

