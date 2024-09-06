function out = nirs_run_E_HyperScanCombineNIRS(job)
%Combine the data of 2 subjects (NIRS.mat) in one file, 
% detector, and source label will be kept identical for the first-subject and 
% replaced as additional sources and detectors for the second subject; for example, 
% second-subject detectors 1-16 will become detectors 17 to 32,(if the first subject have 16 detectors...
% do not support combine eeg or auxiliary file 
%filename prefix 
prefix = 'h'; %for "hyperscancombine"
disp('WORK IN PROGRESS FUNCTION');

NIRS = [];
try    
 load(job.NIRSmat{2,1});
catch
    disp('You must add 2 synchronized NIRS file to combine them for hyperscanning ')
end
NIRS2 = NIRS;
 NIRS = [];
load(job.NIRSmat{1,1});
NIRS1 = NIRS;
%default fields will be define according to file 1, the rest will be modify
%to combine both file 

    %use last step of operation
    lst = length(NIRS1.Dt.fir.pp);
    rDtp1 = NIRS1.Dt.fir.pp(lst).p; % path for files to be processed
    NC1 = NIRS1.Cf.H.C.N;
    fs = NIRS1.Cf.dev.fs;

    rDtp2 = NIRS2.Dt.fir.pp(lst).p; % path for files to be processed
    NC2 = NIRS2.Cf.H.C.N;
    if fs~= NIRS2.Cf.dev.fs;
        disp('Error time frequency for both file are different')
    end

    %Channel to be combine in one project 
    nbSrc = NIRS.Cf.H.S.N; %sources
    nbDet = NIRS.Cf.H.D.N; %detector

    NIRS.Cf.H.C.n = [NIRS1.Cf.H.C.n; NIRS2.Cf.H.C.n];
    temp1 = NIRS1.Cf.H.C.id; %channel first subject  %trik first half subject 1 HbO HbR
    temp2 = NIRS2.Cf.H.C.id; %channel second subject
    temp2(1,:)= temp2(1,:)+ max(temp1(1,:));
    temp2(2,:) = temp2(2,:)+ max(temp1(2,:));
    temp2(3,:) = temp2(3,:)+ max(temp1(3,:));
    NIRS.Cf.H.C.id = [temp1,temp2];
    NIRS.Cf.H.C.wl =  [NIRS1.Cf.H.C.wl,NIRS2.Cf.H.C.wl];
    NIRS.Cf.H.C.gp =  [NIRS1.Cf.H.C.gp;NIRS2.Cf.H.C.gp];
    NIRS.Cf.H.C.ok =  [NIRS1.Cf.H.C.ok;NIRS2.Cf.H.C.ok];
    NIRS.Cf.H.C.N  = [ NIRS1.Cf.H.C.N +  NIRS2.Cf.H.C.N];

    NIRS.Cf.H.S.N = [NIRS1.Cf.H.S.N + NIRS2.Cf.H.S.N];
    NIRS.Cf.H.S.r.o.mm.p = [NIRS1.Cf.H.S.r.o.mm.p,  NIRS2.Cf.H.S.r.o.mm.p];

    NIRS.Cf.H.D.N = [NIRS1.Cf.H.D.N + NIRS2.Cf.H.D.N];
    NIRS.Cf.H.D.r.o.mm.p = [NIRS1.Cf.H.D.r.o.mm.p,  NIRS2.Cf.H.D.r.o.mm.p];
  
    %sort wavelengt first:  HbO first subject, HbO second subject, HbR
    %first subject, HbR second subject
    [B, idxwl] = sort( NIRS.Cf.H.C.wl) ;

    NIRS.Cf.H.C.n = NIRS.Cf.H.C.n(idxwl);
    NIRS.Cf.H.C.id = NIRS.Cf.H.C.id(:,idxwl);
    NIRS.Cf.H.C.wl =  NIRS.Cf.H.C.wl(idxwl);
    NIRS.Cf.H.C.gp =     NIRS.Cf.H.C.gp(idxwl);
    NIRS.Cf.H.C.ok =   NIRS.Cf.H.C.ok(idxwl);
    NIRS.Cf.H.C.N = NIRS.Cf.H.C.N;
     if isempty(job.f_HyperScan_outdir)
         disp('Enter an output folder to save the sessions with both subject otherwise it will be save in the first subject folder');
          [filepath,fil,ext]=fileparts(rDtp1{1});
         job.f_HyperScan_outdir{1} =   filepath;
     end
      
    dirout = job.f_HyperScan_outdir{1};
    if ~isdir(dirout)
        mkdir(dirout);
        disp(['Folder: ' ,dirout,' Created']);
    end
    for f=1:numel(rDtp1) %Loop over all files of a NIRS.mat
         [filepath,fil,ext]=fileparts(rDtp1{f});
         outfile = fullfile(dirout,[prefix fil ext]);
          outfile_vmrk = fullfile(dirout,[prefix fil '.vmrk']);
            outfile_vhdr = fullfile(dirout,[prefix fil '.vhdr']);
        try    
        d1 = fopen_NIR(rDtp1{f},NC1);     
        catch
            disp(['Failed to open file: ' rDtp1{f} ,' if you move directory location use utility folder adjustement '])      
        end
        d2 = fopen_NIR(rDtp2{f},NC2); 
        disp(['Open file A: ' rDtp1{f}])
        disp(['Open file B: ' rDtp2{f}])

        nsample = min([size(d2,2), size(d1,2)]);
      
        d = [d1(:,1:nsample);d2(:,1:nsample)];
        d = d(idxwl,:);
       
       nbChanneloffset = size(d1,1)/2; 
       nbChanneloffsetHBR = size(d1,1);
            %ajust the channel number for the second subject in the ind_dur_ch information according to the fact that usualy channel are presented 
        
        fwrite_NIR(outfile,d);


              
        NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;  

                %transfert file noise 
                [dir1,fil1,~] = fileparts(rDtp1{f});
                vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
                [label1,ind_dur_ch1] = read_vmrk_all(vmrk_path);
                 %HBO channel subject 1 unchanged

                 %HBR  
                 ind_dur_ch1(find(ind_dur_ch1(:,3)>nbChanneloffset),3) =  ind_dur_ch1(find(ind_dur_ch1(:,3)>nbChanneloffset),3) -  nbChanneloffset +  nbChanneloffsetHBR;
             
                [dir2,fil2,~] = fileparts(rDtp2{f});
                vmrk_path = fullfile(dir2,[fil2 '.vmrk']);
                [label2,ind_dur_ch2] = read_vmrk_all(vmrk_path);
                %HBR channel subject 2
                ind_dur_ch2(find(ind_dur_ch2(:,3)>nbChanneloffset),3) =  ind_dur_ch2(find(ind_dur_ch2(:,3)>nbChanneloffset),3) +  nbChanneloffsetHBR;
                %HBO channel subject 2
                ind_dur_ch2(find(ind_dur_ch2(:,3)<nbChanneloffset & ind_dur_ch2(:,3)>0),3) =  ind_dur_ch2(find(ind_dur_ch2(:,3)<nbChanneloffset & ind_dur_ch2(:,3)>0),3) +  nbChanneloffset ; %add offset channel execept for the zero all channel implied

                ind_dur_ch_all = [ind_dur_ch1;ind_dur_ch2];
                label_all = [label1;label2];
                write_vmrk_all( outfile_vmrk,ind_dur_ch_all,label_all);
               


                % try 
                %     info1 = read_vhdr_brainvision((fullfile(dir1,[fil1,'.vhdr'])));
                %      info2 = read_vhdr_brainvision((fullfile(dir2,[fil2,'.vhdr'])));
                %     ChannelLabels = info.label;               
                % 
                % catch
                % end
                    ChannelLabels = ConvertmlIDsrs2label(NIRS);
               
                    nirs_boxy_write_vhdr( outfile_vhdr,... %Output file
                    outfile,... %DataFile
                    outfile_vmrk,... %MarkerFile,...
                    'nirs_run_HyperScanCombine',... %Function that created the header
                    '',... %Channel Resolution
                    '',... %Channel Units
                    ChannelLabels,... %names given as a column of cells
                    1/NIRS.Cf.dev.fs*1000000,... %SamplingInterval in microseconds
                    nsample); %info.datapoint
                                   

        %Video dual file not supported first file value by default 
        %Aux dual file not supported  first file value by default 
        %EEG dual file not supported  first file value by default 
              
              
    end
    NIRS.Dt.fir.pp(lst+1).pre = 'HyperScanCombine';
    NIRS.Dt.fir.pp(lst+1).job = job;
    save(fullfile(dirout,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} =fullfile(dirout,'NIRS.mat');

    out.NIRSmat = job.NIRSmat;

       