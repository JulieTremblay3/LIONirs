function out = nirs_run_writenirshomer(job)
% Write toolbox data from the last node as data as .nirs format 
% Field present in the structure
% d : Array of raw intensity 
% t : Array of data time point
% s : Array of stimulus onset N time points x N conditions
% SD structure
% SD.Lambda Wavelength used for data acquisition
% SD.nSrcs
% SD.nDets
% SD.SrcPos Array of source coordinates nSrcs x 3
% SD.DetPos Array of detector coordinates nDets x 3
% SD.MeasList %list of measurement channel: source idx/ detector id / unused/ wavelength
 %NIRS_exportoption,NIRSname
load(job.NIRSmat{1})     % Matrice NIRS
%Big loop over all subjects

    for Idx=1:size(job.NIRSmat,1)      
        nbch = numel(NIRS.Cf.H.C.n);
        dall = [];
        Sessionid = numel(NIRS.Dt.fir.pp); %prendre la dernière session

    for i = 1:numel(NIRS.Dt.fir.pp(1,Sessionid).p)
       [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(1,Sessionid).p{i});
       %nbf = sprintf('%02.0f',i)%nom du fichier data .nir
       d = fopen_NIR([pathstr,filesep,name,ext],nbch);       
       outfile= [pathstr,filesep,name,'.nirs'];
       
       
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
             t = t';
             aux = zeros(size(s,1),1);
             SD.MeasList = ml;
             save(outfile,'d','SD','t','s','aux','-mat');
             fprintf('Saving %s ...\n', outfile);
    end
    end



out.NIRSmat = job.NIRSmat;