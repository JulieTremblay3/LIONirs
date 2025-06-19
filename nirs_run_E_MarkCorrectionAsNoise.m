function out = nirs_run_E_MarkCorrectionAsNoise(job)
% USED THE SUBSTRACT COMPONENT TO WERE COMPONENT WERE corrected tp 
% 
prefix = 'y'; %add yellow for substracted artefact"

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    load(job.NIRSmat{filenb,1});
    [pathstr, name, ext] = fileparts(job.NIRSmat{filenb,1});
    infileCorrectionApply = fullfile(pathstr,'CorrectionApply.mat');
    lst = numel(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; 
    try
        load(infileCorrectionApply);
 

    NC = NIRS.Cf.H.C.N;
    for f=1:size(rDtp,1)
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
            [ind_dur_ch] = read_vmrk_find( infilevmrk ,'Bad Interval');
             d = fopen_NIR(rDtp{f,1},NC);     
            if isempty(ind_dur_ch)
                noise = zeros(size(d));
            else
                noise = mat2d2ind_dur_ch(ind_dur_ch);
            end
            for idx = 1:numel(PARCORR)
                if PARCORR(idx).file == f
                 
                    if strcmp(PARCORR(idx).type, 'PARAFAC') | strcmp(PARCORR(idx).type, 'PCA') 
                        listgood =  [PARCORR(idx).listgood] ;
                             noise(listgood, PARCORR(idx).indt) = 1;                             
                    end         
             
                end
                   
                
            end
            
            tmp = noise(1:NC/2,:) + noise(NC/2+1:end,:);
            noise = [tmp;tmp]; %ensure both wavelenght have the same noise. 
           % figure;imagesc(noise);

    end
    
    [label,ind_dur_ch] = read_vmrk_all(infilevmrk);
    ind_dur_ch_2 = mat2d2ind_dur_ch(noise);
    label_2 = cell(size(ind_dur_ch_2,1),2);
    label_2(:,1)={'bad_step'};
    label_2(:,2)={'corrected'};
    ind_dur_ch = [ind_dur_ch;ind_dur_ch_2];
    label = [label;label_2];



     [dir1,fil1,ext1] = fileparts(rDtp{f});
    outfilevmrk  =  fullfile(dir1,[prefix fil1 '.vmrk']); 
    write_vmrk_all( outfilevmrk,ind_dur_ch,label);

            infilevhdr = fullfile(dir1,[fil1 '.vhdr']);
            outfilevhdr = fullfile(dir1,[prefix fil1 '.vhdr']);
            copyfile(infilevhdr,outfilevhdr);
           
        outfile = fullfile(dir1,[prefix fil1 ext1]);
        if job.DelPreviousData
            delete(rDtp{f,1});
            delete(infilevmrk);
            delete(infilevhdr);
            disp(['Delete previous .nir data file: ',rDtp{f,1}])
        end

         if f == 1
            NIRS.Dt.fir.pp(lst+1).pre = 'MarkCorrectedInYellow';
            NIRS.Dt.fir.pp(lst+1).job = job;
         end
         NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
        fwrite_NIR(outfile,d);
    save(fullfile(dir1,'NIRS.mat'),'NIRS');
 
       catch
        disp(['No correction found, noise is not modify ', infileCorrectionApply])
    end
       NIRSmat = job.NIRSmat{filenb,1};
end


out.NIRSmat = {NIRSmat};


