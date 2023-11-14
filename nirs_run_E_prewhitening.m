%___________________________________________________________________
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine 
% www.lionlab.umontreal.ca
%___________________________________________________________________
function out = nirs_run_E_prewhitening(job)

%job.NIRSmat = 'D:\Data\ELAN\BB044\Clean\NIRS.mat';
%job.DelPreviousData = 0;
prefix = 'w'; %for "whitening"
DelPreviousData  = job.DelPreviousData;

pmax = 10; %number of iteration 

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
    NIRS = [];
    %load(job.NIRSmat);
    
    load(job.NIRSmat{filenb,1});

    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;

    [dir2,~,~] = fileparts(job.NIRSmat{filenb,1});
    %[dir2,tmp,tmp] = fileparts(job.NIRSmat);

    for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat 

        d = fopen_NIR(rDtp{f,1},NC);
        %d = fopen_NIR('Y:\NIRx-actiCHAMP_Analysed\ELAN\Resting2021\Rae\Data\BB044\Clean\maBB044.nir',NC);
 

        %LogL = nan(pmax,NC);
        BIC = nan(pmax,NC);
        for i = 1:NC
            disp(i)
            y = d(i,:)';
            for p = 1:pmax
                tic
                mdl = regARIMA(p,0,0);
                [~,~,LogL] = estimate(mdl,y,'Display','off');
                [~,BIC(p,i)] = aicbic(LogL,i+1,length(y)-i);
                if p>5
                    if max(BIC(p-5:p,i)) == BIC(p,i)
                        break
                    end
                end
                toc
            end
        end

        dw = nan(size(d));
        for i = 1:NC
            disp(i)
            y = d(i,:)';
            p = find(min(BIC(:,i)));
            mdl = regARIMA(p,0,0);
            estMdl = estimate(mdl,y,'Display','off');
            dw(i,:) = infer(estMdl,y);
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

        fwrite_NIR(outfile,dw); 
        fprintf('%s\n',outfile);

         %write new .vmrk file 
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
    end
end
 out.NIRSmat = job.NIRSmat;

