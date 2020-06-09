function out = nirs_run_E_detrend(job)
%Remode linear trend
%filename prefix 
prefix = 'd'; %for "detrend"
DelPreviousData  = job.DelPreviousData;


for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
%     try
        NIRS = [];
        load(job.NIRSmat{filenb,1});
        [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});

        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        fprintf('%s\n','File processed');
        for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat 
%             try
                d = fopen_NIR(rDtp{f,1},NC);              
                [dir1,fil1,~] = fileparts(rDtp{f});
                intensnorm = d';
                X = 1:1:size(intensnorm,1);
                Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
                Mb2 =  intensnorm(1,:)'; %offset    
                A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
                d = (intensnorm - A)';       
                
                [dir1,fil1,ext1] = fileparts(rDtp{f});
                infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
                try
                infilevhdr = fullfile(dir1,[fil1 '.vhdr']);    
                catch
                end
                
                if ~exist(dir2,'dir'), mkdir(dir2); end
                outfile = fullfile(dir2,[prefix fil1 ext1]);
                outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
                outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);

               
                fwrite_NIR(outfile,d); 
                fprintf('%s\n',outfile);
 
                %write new .vmrk file 
                try
                copyfile(infilevmrk,outfilevmrk);
                catch;end
                try
                    ChannelLabels = ConvertmlIDsrs2label(NIRS);
                    SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
                    nirs_boxy_write_vhdr(outfilevhdr,... %Output file
                        outfile,... %DataFile
                        outfilevmrk,... %MarkerFile,...
                        'nirs_E_detrend',... %Function that created the header
                        '',... %Channel Resolution
                        '',... %Channel Units
                        ChannelLabels,... %names given as a column of cells
                        SamplingInterval,...
                        size(d,2)); %SamplingInterval in microseconds
                catch
                end
                 if DelPreviousData
                    delete(rDtp{f,1});
                    delete(infilevmrk)
                    delete(infilevhdr)
                end
                %add outfile name to NIRS
                if f == 1
                    NIRS.Dt.fir.pp(lst+1).pre = 'Detrend';
                    NIRS.Dt.fir.pp(lst+1).job = job;
                end 
                NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
        end
            save(fullfile(dir2,'NIRS.mat'),'NIRS');
            job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');

end
out.NIRSmat = job.NIRSmat;
                        
                        
                        