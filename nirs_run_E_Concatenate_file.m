function out = nirs_run_E_Concatenate_file(job)
% Add raw file one after the other
% This function is built to merge several bloc in the same NIRS.mat
% Look for Concatenate_NIRS to merge different NIRS.mat restriction of only test for one bloc... used concatenate file before.
prefix = 'All';


for filenb=1:1 %only one NIRS.mat merge file inside
    NIRS = [];
    load(job.NIRSmat{filenb,1});
    [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});
    
    %use last step of preprocessing
    lst = length(NIRS.Dt.fir.pp);
    rDtptmp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    idselected = job.e_Concatenate_blocid;
    if idselected==0
        rDtp =rDtptmp; %keep all same order
        idselected = 1 : numel(rDtp);
    else
        rDtp =rDtptmp(idselected,1);
    end
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    dall = [];
    ind_dur_chtmp =[];
    labeltmp=[] ;
    if isfield(NIRS.Dt,'EEG')
        moduleEEG =  numel(NIRS.Dt.EEG(end).pp);
    end
    if isfield(NIRS.Dt,'AUX')
        moduleaux = numel(NIRS.Dt.AUX(end).pp);
    end
    if isfield(NIRS.Dt,'Video')
        NIRS.Dt = rmfield(NIRS.Dt,'Video');
    end
    %retro compatibilité find size block
    if 1 %~isfield(NIRS.Dt.fir,'sizebloc')
        for f=1:size(rDtp,1)
         [pathstr, name, ext] = fileparts(rDtp{f,1});
         d = fopen_NIR(rDtp{f,1},NC);
         NIRS.Dt.fir.sizebloc{f}=size(d,2); 
         end
    end
%     if numel(NIRS.Dt.fir.sizebloc)~=numel(NIRS.Dt.fir.aux5)
%         for f=1:size(rDtp,1)
%             [pathstr, name, ext] = fileparts(rDtp{f,1});
%             d = fopen_NIR(rDtp{f,1},NC);
%             NIRS.Dt.fir.sizebloc{f} = size(d,2); 
%         end
%     end
    
    if isfield(NIRS.Dt.fir,'aux5')
        aux5 = [];
        offsetsizebloc = 0;
        for f=1:numel(rDtp) %size(rDtp,1)
            val =    NIRS.Dt.fir.aux5{f};
            if f>1
                NIRS.Dt.fir.sizebloc{f-1}-1;
                offsetsizebloc = offsetsizebloc +  NIRS.Dt.fir.sizebloc{f-1};
                val(:,2)= val(:,2)+ offsetsizebloc ;
            end
            aux5 =[aux5;val ];
            
        end
        NIRS.Dt.fir.aux5 = {aux5};
        clear aux5
    end
    if isfield(NIRS.Dt.fir,'aux4')
        aux = [];
        for f=1:size(rDtp,1);aux =[aux,NIRS.Dt.fir.aux4{f} ]; end
        NIRS.Dt.fir.aux4 = {aux};        clear aux
    end
    if isfield(NIRS.Dt.fir,'aux3')
        aux = [];
        for f=1:size(rDtp,1);aux =[aux,NIRS.Dt.fir.aux3{f} ]; end
        NIRS.Dt.fir.aux3 = {aux};        clear aux
    end
    if isfield(NIRS.Dt.fir,'aux2')
        aux = [];
        for f=1:size(rDtp,1);aux =[aux,NIRS.Dt.fir.aux2{f} ]; end
        NIRS.Dt.fir.aux2 = {aux};        clear aux
    end
    if isfield(NIRS.Dt.fir,'aux1')
        aux = [];
        for f=1:size(rDtp,1);aux =[aux,NIRS.Dt.fir.aux1{f} ]; end
        NIRS.Dt.fir.aux1 = {aux};        clear aux
    end
    if isfield(NIRS.Dt.fir,'aux5i')
        aux = [];
        for f=1:size(rDtp,1);aux =[aux,NIRS.Dt.fir.aux5i{f} ]; end
        NIRS.Dt.fir.aux5i = {aux};        clear aux
    end
    
    if isfield(NIRS.Dt.fir,'sizebloc')
        totsizebloc = 0;
        for f=1:size(rDtp,1)
            totsizebloc=totsizebloc+NIRS.Dt.fir.sizebloc{f};
        end
        newsizebloc = totsizebloc;
    end
    if size(rDtp,2)>size(rDtp,1)
        rDtptmp = rDtp';
        clear rDtp;
        rDtp = rDtptmp;clear rDtptmp;
    end
    offsetsizebloc = 0;
    for f=1:size(rDtp,1) %For every path of NIRS.mat file
        [pathstr, name, ext] = fileparts(rDtp{f,1});
        d = fopen_NIR(rDtp{f,1},NC); %Load whole bloc
        infilevmrk = fullfile(pathstr,[name,'.vmrk']);
        infilevhdr = fullfile(pathstr,[name,'.vhdr']);
        if job.m_Concatenate_option == 1 % detrend each bloc data segment
            intensnorm = d';
            X = 1:1:size(intensnorm,1);
            Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
            Mb2 =  intensnorm(1,:)'; %offset
            A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
            d = (intensnorm - A)';
        end
        noise = logical(zeros(size(d)));
        %          [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
        %     if ~isempty(ind_dur_ch)
        %         maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
        %         badind = find(maxpoint>size(noise,1));
        %         if ~isempty(badind)
        %             disp(['Warning file ' vmrk_path ' marker : ' num2str(badind') ' are out of range'])
        %             ind_dur_ch(badind,2)=size(noise,2)- ind_dur_ch(badind,1);
        %         end
        %         for Idx = 1:size(noise,2)
        %             mrks = find(ind_dur_ch(:,3)==Idx);
        %             ind = ind_dur_ch(mrks,1);
        %             indf = ind + ind_dur_ch(mrks,2) - 1;
        %             if ~isempty(ind)
        %                 try
        %                     for i = 1:numel(ind)
        %                         noise(ind(i):indf(i),Idx) = 1;
        %                     end
        %                 catch
        %                     msgbox('Noise reading problem')
        %                 end
        %             end
        %         end
        %     end
        
        
        dall = [dall, d];
        [label_all,ind_dur_ch_all] = read_vmrk_all(infilevmrk);
        if f>1
            offsetsizebloc = offsetsizebloc+NIRS.Dt.fir.sizebloc{f-1};
            ind_dur_ch_all(:,1)=ind_dur_ch_all(:,1)+ offsetsizebloc;
        end
        ind_dur_chtmp = [ind_dur_chtmp;ind_dur_ch_all];
        labeltmp = [labeltmp;label_all];
    end
    if isfield(NIRS.Dt,'AUX') %Concatenate auxilary
        for iaux = 1:numel(NIRS.Dt.AUX)
            AUX{iaux}.data = [];
            AUX{iaux}.marker = [];
            AUX{iaux}.ind_dur_ch = [];
            offsetsizebloc = 0;
            
            
            for f=1:size(rDtp,1)
                fileAUX = NIRS.Dt.AUX(iaux).pp(moduleaux).p{idselected(f)};
                tstartf = NIRS.Dt.AUX(iaux).pp(moduleaux).sync_timesec{idselected(f)};
                tstopf = tstartf+ size( d,2)*1/NIRS.Cf.dev.fs;
                [data,infoBV,marker,ind_dur_ch]= fopen_EEG(fileAUX,tstartf , tstopf);
                if f>1
                    offsetsizebloc = offsetsizebloc+NIRS.Dt.fir.sizebloc{idselected(f-1)};
                    ind_dur_ch(:,1)=ind_dur_ch(:,1)+ offsetsizebloc;
                end
                if job.m_Concatenate_option == 1 % detrend each bloc data segment
                    intensnorm = data;
                    X = 1:1:size(intensnorm,1);
                    Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
                    Mb2 =  intensnorm(1,:)'; %offset
                    A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
                    data = (intensnorm - A);
                end
                AUX{iaux}.data=[AUX{iaux}.data;data];
                AUX{iaux}.marker =  [AUX{iaux}.marker;marker];
                AUX{iaux}.ind_dur_ch = [AUX{iaux}.ind_dur_ch;ind_dur_ch];
            end
            [dir1,fil1,ext1] = fileparts(fileAUX);
            outfileAUX = fullfile(dir1,[prefix fil1 ext1]);
            
            AUX{iaux}.infoBV=infoBV;
            AUX{iaux}.infoBV.DataPoints = size(AUX{iaux}.data,1);
            AUX{iaux}.outfile =  outfileAUX;
            NIRS.Dt.AUX(iaux).pp(moduleaux+1).p{1} = outfileAUX;
            NIRS.Dt.AUX(iaux).pp(moduleaux+1).sync_timesec{1} = 0;
            fwrite_EEG(outfileAUX, AUX{iaux},1,AUX{iaux}.infoBV.DataPoints );
        end
    end
    if isfield(NIRS.Dt,'EEG') %Concatenate EEG
        for iEEG = 1:numel(NIRS.Dt.EEG)
            EEG{iEEG}.data = [];
            EEG{iEEG}.marker = [];
            EEG{iEEG}.ind_dur_ch = [];
            offsetsizebloc = 0;
            for f=1:size(rDtp,1)
                fileEEG = NIRS.Dt.EEG(iEEG).pp(end).p{idselected(f)}
                 tstartf = NIRS.Dt.EEG(iEEG).pp(moduleaux).sync_timesec{idselected(f)};
                tstopf = tstartf+ size( d,2)*1/NIRS.Cf.dev.fs;
                
                [data,infoBV,marker,ind_dur_ch]= fopen_EEG(fileEEG, tstartf , tstopf);
                if f>1
                    offsetsizebloc = offsetsizebloc+NIRS.Dt.fir.sizebloc{f-1}-1;
                    ind_dur_ch(:,1)=ind_dur_ch(:,1)+ offsetsizebloc;
                end
                EEG{iEEG}.data=[EEG{iEEG}.data;data];
                EEG{iEEG}.marker =  [EEG{iEEG}.marker;marker];
                EEG{iEEG}.ind_dur_ch = [EEG{iEEG}.ind_dur_ch;ind_dur_ch];
            end
            [dir1,fil1,ext1] = fileparts(fileEEG);
            outfileEEG = fullfile(dir1,[prefix fil1 ext1]);
            
            EEG{iEEG}.infoBV=infoBV;
            EEG{iEEG}.infoBV.DataPoints = size(EEG{iEEG}.data,1);
            EEG{iEEG}.outfile =  outfileEEG;
            NIRS.Dt.EEG(iEEG).pp(moduleEEG+1).p{1} = outfileEEG;
            NIRS.Dt.EEG(iEEG).pp(moduleEEG+1).sync_timesec{1} = 0;
            fwrite_EEG(outfileEEG, EEG{iEEG},1,EEG{iEEG}.infoBV.DataPoints );
        end
    end
    
    
    [dir1,fil1,ext1] = fileparts(rDtp{1,1});
    NIRS.Dt.fir.sizebloc{1} =newsizebloc;
    %.nir format
    % if NewDirCopyNIRS
    %     dir2 = [dir1 filesep NewNIRSdir];
    if ~exist(dir2,'dir'), mkdir(dir2); end
    outfile = fullfile(dir2,[prefix fil1 ext1]);
    outfile2 = fullfile(dir2,[prefix fil1 '_std' ext1]);
    outfile3 = fullfile(dir2,[prefix fil1 '_events' ext1]);
    outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
    outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);
    fileOut_nir = fullfile(dir2,[prefix fil1 '.nir']);
    %         else
    %             outfile = fullfile(dir1,[prefix fil1 ext1]);
    %             outfile2 = fullfile(dir1,[prefix fil1 '_std' ext1]);
    %             outfile3 = fullfile(dir1,[prefix fil1 '_events' ext1]);
    %             outfilevmrk = fullfile(dir1,[prefix fil1 '.vmrk']);
    %             outfilevhdr = fullfile(dir1,[prefix fil1 '.vhdr']);
    %             fileOut_nir = fullfile(dir1,[prefix fil1 '.nir']);
    %         end
    
    
    fwrite_NIR(outfile,dall);
    write_vmrk_all(outfilevmrk,ind_dur_chtmp,labeltmp);
    
    %write outvhdr file
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
        'nirs_Concatenate_File',... %Function that created the header
        '',... %Channel Resolution
        '',... %Channel Units
        ChannelLabels,... %names given as a column of cells
        SamplingInterval,...
        NIRS.Dt.fir.sizebloc{1}); %SamplingInterval in microseconds
    
    clear dall
    %si plusieur fichier le premier correspondant au .nirs
    
    [dir1,fil1,ext1] = fileparts(rDtp{f,1});
    NIRS.Dt.fir.pp(lst+1).p = {outfile};
    NIRS.Dt.fir.pp(lst+1).pre = 'Concatenate File';
    NIRS.Dt.fir.pp(lst+1).job = job;

    NIRS.Cf.H.C.ok(:,1) = (sum(NIRS.Cf.H.C.ok,2)>=size(NIRS.Cf.H.C.ok,2)/2);  %Ensure tant channel in less than 50 % of the bloc are kept as remove channel
    save(fullfile(dir2,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');
    
end


out.NIRSmat = job.NIRSmat;
