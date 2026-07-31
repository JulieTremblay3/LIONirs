function out = nirs_run_E_autoregressivemodel(job)
% %use the fonction fillgaps to replace nan by extrapoled data.
% 
% %filename prefix
prefix = 'ar'; %for "autoregressivemodel
DelPreviousData  = job.DelPreviousData;
% 
% 
% 
% 
 for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
     tic
     NIRS = [];
     load(job.NIRSmat{filenb,1});
        [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});


     %use last step of preprocessing
     lst = length(NIRS.Dt.fir.pp);
     rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
     NC = NIRS.Cf.H.C.N;
     fs = NIRS.Cf.dev.fs;

     for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat
         d = fopen_NIR(rDtp{f,1},NC);
         samp_length = size(d,2);
         time = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:samp_length*1/NIRS.Cf.dev.fs;
                [dir1,fil1,ext1] = fileparts(rDtp{f});
               
        [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(lst).p{f});
        [dir1,fil1,ext1] = fileparts(rDtp{f});
        infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
        infilevhdr = fullfile(dir1,[fil1 '.vhdr']);


           
        rejectionThreshold = 0.3
        mrk_type_arr='bad_step'
        [ind_dur_ch] = read_vmrk_find(infilevmrk,mrk_type_arr);
    
        [cleanOD, noiseMask, rejectedChannels] = processFIRS_AR_Interpolation(d', ind_dur_ch, rejectionThreshold,job.e_nbsampleprediction, job.e_autoregressivemodelorder);
                 
            
                %if NewDirCopyNIRS
                 %   dir2 = [dir1 filesep NewNIRSdir];
                    if ~exist(dir2,'dir'), mkdir(dir2); end 
                    outfile = fullfile(dir2,[prefix fil1 ext1]);   
                    outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
                    outfilevhdr = fullfile(dir2,[prefix fil1  '.vhdr']);
                    fwrite_NIR(outfile,cleanOD');                 
                    try
                       copyfile(infilevmrk,outfilevmrk);
                       copyfile(infilevhdr,outfilevhdr); 
                    catch
                    end

        
        
        
                if DelPreviousData
                    delete(rDtp{f,1});
                    delete(infilevmrk);
                    delete(infilevhdr);
                    disp(['Delete previous .nir data file: ',rDtp{f,1}]);
                end

                %add outfile name to NIRS
                if f == 1
                    NIRS.Dt.fir.pp(lst+1).pre = 'AutoregressiveModel';
                    NIRS.Dt.fir.pp(lst+1).job = job;
                end 
                NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile; 
                fprintf('%s\n',outfile);

     end

         save(fullfile(dir2,'NIRS.mat'),'NIRS');   
        job.NIRSmat{filenb,1}=fullfile(dir2,'NIRS.mat');


 end

out.NIRSmat = job.NIRSmat;
 
