function out = nirs_run_writeSNIRF(job)
%NIRS_RUN_READSNIRF This function takes in a job and creates an analysis
%step from it. It will read the snirf file and write a NIRS.mat construct
%that can be used in the next steps of an analysis.
%   prjname: .prj file for project
%   fileIn: raw data file in .snirf format
%   fileOut: .nir, .mat, .vhdr, .vmkr

% NOTE: While based on nirs_run_readNIRxscout, DATA.d in nirs_run_readSNIRF
% corresponds to the transposed of DATA.d in nirs_run_readNIRxscout. That
% means that while original code use the value of DATA.d as is, code reused
% from readNIRxscout has the DATA.d variable transposed (and the size
% function has the 1 argument switched to 2).

load(job.NIRSmat{1})     % Matrice NIRS
%Big loop over all subjects

    for Idx=1:size(job.NIRSmat,1)      
        nbch = numel(NIRS.Cf.H.C.n);
        dall = [];
        Sessionid = numel(NIRS.Dt.fir.pp); %prendre la dernière session

    for i = 1:numel(NIRS.Dt.fir.pp(1,Sessionid).p)
       [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(1,Sessionid).p{i});
       d = fopen_NIR([pathstr,filesep,name,ext],nbch);       
        if ~isdir(job.f_writeNIRSdir{1})
           mkdir(job.f_writeNIRSdir{1})
       end
       if ~isdir([job.f_writeNIRSdir{1},filesep,name])
           mkdir([job.f_writeNIRSdir{1},filesep,name])
       end       
        outfile= [job.f_writeNIRSdir{1},filesep,name,filesep,name,'.snirf'];
        d = d'; %data time point x channels

            SD = [];
            SD.Lambda = NIRS.Cf.dev.wl;
            SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
            SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
            SD.nSrcs = size(SD.SrcPos,1);
            SD.nDets = size(SD.DetPos,1);
            if numel(NIRS.Cf.H.C.id(2,:)) == size(d,2)
                ml = [NIRS.Cf.H.C.id(2,:)',NIRS.Cf.H.C.id(3,:)', ones(numel(NIRS.Cf.H.C.id(2,:)),1),NIRS.Cf.H.C.wl']; %srs,det,1,wav
                t = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(d,1)*1/NIRS.Cf.dev.fs;
            end       
            
             aux5 =  NIRS.Dt.fir.aux5{1}; %trigger info to set in s structure
             %Example; [2, 9;
             %2, 1733;
             %1,1]
             %Trig 2 sample of time 2 et 1733, 
             %Trig 1 sample 
             nbmaxtrig =  max(aux5(:,1));
             s = zeros(numel(t), nbmaxtrig);
             for i = 1:size(aux5,1)
                 s(aux5(i,2),aux5(i,1))=1;
             end
                 
              SD.MeasList = ml;
            
             %need to have and aux file
            aux = zeros(size(d,1),1);
    end

         nirs= struct('d',d,'SD',SD,'t',t','s',s,'aux',aux);
        % Save SNIRF: Convert .nirs format data to SnirfClass object, save it to .snirf file (HDF5)
        fprintf('Saving %s ...\n', outfile);
        snirf_saved = SnirfClass(nirs);
        tic; snirf_saved.Save(outfile); toc
    end

out.NIRSmat = job.NIRSmat;

