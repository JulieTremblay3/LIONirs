function out = nirs_run_E_Concatenate_nirsmat(job)
% Add nirs.mat raw data file on after the other
% This function is built to merge several NIRS.mat, use Concatenate file
% before to support one big bloc by NIRS.mat
%%% NOT ADAPTED YET

 
prefix = 'All'; 
[~,~,ext] =fileparts(job.f_nirsmatinfo{1});
if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
    [data, text, rawData] = xlsread(job.f_nirsmatinfo{1});
elseif strcmp(ext,'.txt')
    [data, text, rawData] = readtxtfile_asxlsread(job.f_nirsmatinfo{1});
end
    
NIRSDtp = rawData(2:end,1);
ListDtp = rawData(2:end,2);
[pathoutlist, namelist, ext] = fileparts(job.f_nirsmatinfo{1});

temp = ListDtp{1}
if strcmp(temp(end-3:end),'zone')% USE ZONE TEMPLATE FIRST TO LOOK IN OTHER FILE    
    dall = [];
    ind_dur_chtmp =[];
    labeltmp=[] ;
    aux5 = []; %keep as important structure built in trigger
    sizebloc = 0;
    for filenb=1:size(NIRSDtp,1) %size(job.NIRSmat,1) %For every specified NIRS.mat file
        NIRS = [];            
        load(fullfile(NIRSDtp{filenb},'NIRS.mat'));    
        lst = length(NIRS.Dt.fir.pp);
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs; 
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        
    if strfind(NIRS.Dt.fir.pp(lst).pre,'Epoch averaging') %do only first file
        nmax =1;
    else
        nmax=numel(rDtp); %do all file
    end
       if filenb==1
            NIRSref = NIRS;        
            ZoneTemplate = load(fullfile(pathoutlist,ListDtp{1}),'-mat');    
            NCTemplate = NIRS.Cf.H.C.N;
       end
     ZoneHoy = load(fullfile(pathoutlist,ListDtp{filenb}),'-mat');

     for f=1:nmax %for each bloc of the NIRS.mat file (last processing step)
        d = fopen_NIR(rDtp{f,1},NC); %Load whole bloc   
         dzone = zeros(NCTemplate,size(d,2));
        [pathstr,name,ext]=fileparts(rDtp{f,1});
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
        
        for izone = 1:numel(ZoneTemplate.zone.label)
            chtemplateHbO = ZoneTemplate.zone.plotLst{izone};
            chtemplateHbR = ZoneTemplate.zone.plotLst{izone}+NCTemplate/2;
            for izonenow = 1:numel(ZoneHoy.zone.label)
                if strcmp(ZoneTemplate.zone.label{izone},ZoneHoy.zone.label{izonenow})
                    chHbO = ZoneHoy.zone.plotLst{izonenow};
                    chHbR = ZoneHoy.zone.plotLst{izonenow}+NC/2 ;
                    break
                end
            end
            dzone(chtemplateHbO,:) = ones(numel(chtemplateHbO),1)*nanmean(d(chHbO,:),1);
            dzone(chtemplateHbR,:) = ones(numel(chtemplateHbR),1)*nanmean(d(chHbR,:),1);

        end
        dall = [dall, dzone];
        %figure;plot(dzone')
       find(double(sum(abs(dall),2)>0));
          [label_all,ind_dur_ch_all] = read_vmrk_all(infilevmrk);   
        %transfert in matrix to adjust channel
        
        ind_dur_ch_all(:,1)=ind_dur_ch_all(:,1)+ sizebloc;       
        ind_dur_chtmp = [ind_dur_chtmp;ind_dur_ch_all];
        labeltmp = [labeltmp;label_all];
        
       if isfield(NIRS.Dt.fir,'aux5')           
            val =    NIRS.Dt.fir.aux5{f};              
            val(:,2)= val(:,2)+  sizebloc;               
            aux5 =[aux5;val ];
       end
          sizebloc  = sizebloc + size(d,2);
      end
    end
        newsizebloc = size(dall,2);
        [dir1,fil1,ext1] = fileparts(rDtp{1,1});
        fil1 = namelist;                
        outfile = fullfile(pathoutlist,[prefix fil1 ext1]);
        outfilevmrk = fullfile(pathoutlist,[prefix fil1 '.vmrk']);
        outfilevhdr = fullfile(pathoutlist,[prefix fil1 '.vhdr']);
        fileOut_nir = fullfile(pathoutlist,[prefix fil1 '.nir']);           
        fwrite_NIR(outfile,dall);
        write_vmrk_all(outfilevmrk,ind_dur_chtmp,labeltmp);
        
       
        NIRS = NIRSref;
        if isfield(NIRS.Dt.fir,'pp')
            NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'pp');
            NIRS.Dt.fir.pp(1).p = {outfile};
            NIRS.Dt.fir.pp(1).pre = 'Concatenate File';
            NIRS.Dt.fir.pp(1).job = job;   
         end
         if isfield(NIRS.Dt.fir,'aux5')
              NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux5');
              NIRS.Dt.fir.aux5={aux5};
         end
     
         NIRS.Dt.fir.sizebloc{1} = newsizebloc;
             NIRS.Cf.H.C.ok = double(sum(abs(dall),2)>0);
             
          %clean up old aux file format use export file format for auxiliary 
          if isfield(NIRS.Dt.fir,'aux1');NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux1');end      
          if isfield(NIRS.Dt.fir,'aux2'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux2');end     
          if isfield(NIRS.Dt.fir,'aux3'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux3');end     
          if isfield(NIRS.Dt.fir,'aux4'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux4');end     
          if isfield(NIRS.Dt.fir,'aux5ini'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux5ini');end   
          if isfield(NIRS.Dt.fir,'aux5i'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux5i');end     
          if isfield(NIRS.Dt.fir,'aux5bloc'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux5bloc');end     

          nirs_boxy_write_vhdr(outfilevhdr,... %Output file
            fileOut_nir,... %DataFile
            outfilevmrk,... %MarkerFile,...
            'nirs_convert_boxy',... %Function that created the header
            '',... %Channel Resolution
            '',... %Channel Units
            NIRS.Cf.H.C.n,... %names given as a column of cells
            1/NIRS.Cf.dev.fs*1e6,... %SamplingInterval in microseconds
            NIRS.Dt.fir.sizebloc{1}); %SamplingInterval in microseconds 
    
    
    
    
else %CASE LIST CHANNEL 
    
    
    
    

dall = [];
ind_dur_chtmp =[];
labeltmp=[] ;
aux5 = []; %keep as important structure built in trigger
sizebloc = 0;
for filenb=1:size(NIRSDtp,1) %size(job.NIRSmat,1) %For every specified NIRS.mat file
    NIRS = [];            
    load(fullfile(NIRSDtp{filenb},'NIRS.mat'));    
    lst = length(NIRS.Dt.fir.pp);
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs; 
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    if strfind(NIRS.Dt.fir.pp(lst).pre,'Epoch averaging') %do only first file
        nmax =1;
    else
        nmax=numel(rDtp); %do all file
    end
        %open channel list and check associate channel
        ML_new= [NIRS.Cf.H.C.id(2:3,:)',...
        ones(size(NIRS.Cf.H.C.id,2),1),...
        [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
        fid = fopen(fullfile(pathoutlist,ListDtp{filenb}));
        try
        chlist = textscan(fid, '%s%s');
        catch
        disp(['ERROR Please ensure that channel list ', fullfile(pathoutlist,ListDtp{filenb}) ' exist'])       
        end
        fclose(fid);
        DetL= chlist{1};
        SrsL= chlist{2};
        list = zeros(numel(chlist{1}),1);
        for i=1:numel(chlist{1})
            NIRS.Cf.dev.n;
            NIRSDtp{filenb};
            switch NIRS.Cf.dev.n
                case 'ISS Imagent'
                    SDdetL = StrBoxy2SDDet_ISS(DetL{i});
                    SDsrsL = StrBoxy2SDPairs(SrsL{i});
                case 'NIRx' 
                    SDdetL = StrBoxy2SDDet(DetL{i});
                     tmp = SrsL{i};
                     SDsrsL =str2num(tmp(2:end));
                case 'NIRS FILE HOMER'
                     SDdetL = StrBoxy2SDDet(DetL{i});
                     tmp = SrsL{i};
                     SDsrsL =str2num(tmp(2:end));
                otherwise
                    SDdetL = StrBoxy2SDDet_ISS(DetL{i});
                    SDsrsL = StrBoxy2SDPairs(SrsL{i});
            end 
            L1 = find(ML_new(:,1)==SDsrsL & ML_new(:,2)==SDdetL & ML_new(:,4)==1);
            L2 = find(ML_new(:,1)==SDsrsL & ML_new(:,2)==SDdetL & ML_new(:,4)==2);
            if isempty(L1)
                sprintf(['check ', DetL{i},' ', SrsL{i}]);
                listname{i,1} = [DetL{i} ' ' SrsL{i}];
            else
                listHBOch(i,1)= L1;
                listHBRch(i,1)= L2;
                listname{i,1} = [DetL{i} ' ' SrsL{i}];
                listnameHbO{i,1} = [DetL{i} ' ' SrsL{i},'HbO'];
                listnameHbR{i,1} = [DetL{i} ' ' SrsL{i},'HbR'];
            end
        end
        chlst = [listHBOch ;listHBRch]; %list finale for this subject using channel list
        if filenb==1
            NIRSref = NIRS;        
            NIRSref.Cf.H.C.n =  [listnameHbO;listnameHbR] ; %NIRS.Cf.H.C.n(chlst);
            NIRSref.Cf.H.C.N = numel(chlst);
            NIRSref.Cf.H.C.id =  NIRS.Cf.H.C.id(:,chlst);
            NIRSref.Cf.H.C.wl =  NIRS.Cf.H.C.wl(:,chlst);
            NIRSref.Cf.H.C.gp =  NIRS.Cf.H.C.gp(chlst,:);
            NIRSref.Cf.H.C.ok =  NIRS.Cf.H.C.ok(chlst,:);
        end
        try 
        NIRSref.Cf.H.prj = NIRS.Cf.H.prj;
        catch;end
     for f=1:nmax %for each bloc of the NIRS.mat file (last processing step)
        d = fopen_NIR(rDtp{f,1},NC); %Load whole bloc   
        [pathstr,name,ext]=fileparts(rDtp{f,1});
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
        
         
       %Normalized do not consider artifact channels. 
       idexclude= find(NIRS.Cf.H.C.ok(:,f)==0);
       if job.m_Concatenate_Exclude == 0
           d(idexclude,:) = nan;
       elseif job.m_Concatenate_Exclude == 1 %keep all
           
       end
       
       
       
       
       if job.m_Concatenate_Normalized==0 %noting
    
            dall = [dall, d(chlst,:)];
       elseif job.m_Concatenate_Normalized==1 %minmax  
            tmp = d(chlst,:);
            minval = nanmin(nanmin(tmp(:)));
            maxval = nanmax(nanmax(tmp(:)));
            Newmax = 1;
            Newmin = 0;
            dminmax =   (tmp - minval)./(maxval - minval) * (Newmax- Newmin) +  Newmin;
            figure;plot(dminmax');
      
            
       elseif job.m_Concatenate_Normalized==2 %zscore
           dall = [dall, d(chlst,:)];
           tmp = d(chlst,:);
           meanval = nanmean(tmp(:));
           stdval = nanstd(tmp(:));
           dzscore= (d(chlst,:)-meanval)./ stdval;
           dall = [dall, dzscore];
       end
            
         
         
        [label_all,ind_dur_ch_all] = read_vmrk_all(infilevmrk);   
        %transfert in matrix to adjust channel
        
        ind_dur_ch_all(:,1)=ind_dur_ch_all(:,1)+ sizebloc;       
        ind_dur_chtmp = [ind_dur_chtmp;ind_dur_ch_all];
        labeltmp = [labeltmp;label_all];
        
       if isfield(NIRS.Dt.fir,'aux5')           
            val =    NIRS.Dt.fir.aux5{f};              
            val(:,2)= val(:,2)+  sizebloc;               
            aux5 =[aux5;val ];
       end
          sizebloc  = sizebloc + size(d,2);
         
     end
    
end
        newsizebloc = size(dall,2);
        [dir1,fil1,ext1] = fileparts(rDtp{1,1});
        fil1 = namelist;                
        outfile = fullfile(pathoutlist,[prefix fil1 ext1]);
        outfilevmrk = fullfile(pathoutlist,[prefix fil1 '.vmrk']);
        outfilevhdr = fullfile(pathoutlist,[prefix fil1 '.vhdr']);
        fileOut_nir = fullfile(pathoutlist,[prefix fil1 '.nir']);           
        fwrite_NIR(outfile,dall);
        write_vmrk_all(outfilevmrk,ind_dur_chtmp,labeltmp);
        
       
        NIRS = NIRSref;
        if isfield(NIRS.Dt.fir,'pp')
            NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'pp');
            NIRS.Dt.fir.pp(1).p = {outfile};
            NIRS.Dt.fir.pp(1).pre = 'Concatenate File';
            NIRS.Dt.fir.pp(1).job = job;   
         end
         if isfield(NIRS.Dt.fir,'aux5')
              NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux5');
              NIRS.Dt.fir.aux5={aux5};
         end
     
         NIRS.Dt.fir.sizebloc{1} = newsizebloc;
             NIRS.Cf.H.C.ok = ones(numel(chlst),1);
          %clean up old aux file format use export file format for auxiliary 
          if isfield(NIRS.Dt.fir,'aux1');NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux1');end      
          if isfield(NIRS.Dt.fir,'aux2'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux2');end     
          if isfield(NIRS.Dt.fir,'aux3'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux3');end     
          if isfield(NIRS.Dt.fir,'aux4'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux4');end     
          if isfield(NIRS.Dt.fir,'aux5ini'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux5ini');end   
          if isfield(NIRS.Dt.fir,'aux5i'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux5i');end     
          if isfield(NIRS.Dt.fir,'aux5bloc'); NIRS.Dt.fir = rmfield(NIRS.Dt.fir, 'aux5bloc');end     

          nirs_boxy_write_vhdr(outfilevhdr,... %Output file
            fileOut_nir,... %DataFile
            outfilevmrk,... %MarkerFile,...
            'nirs_convert_boxy',... %Function that created the header
            '',... %Channel Resolution
            '',... %Channel Units
            NIRS.Cf.H.C.n,... %names given as a column of cells
            1/NIRS.Cf.dev.fs*1e6,... %SamplingInterval in microseconds
            NIRS.Dt.fir.sizebloc{1}); %SamplingInterval in microseconds
end
if isfield(NIRS.Dt,'AUX')
   NIRS.Dt = rmfield(NIRS.Dt,'AUX');
end
if isfield(NIRS.Dt,'EEG')
   NIRS.Dt = rmfield(NIRS.Dt,'EEG');
end
if isfield(NIRS.Dt,'Video')
   NIRS.Dt = rmfield(NIRS.Dt,'Video');
end
if isfield(NIRS.Dt,'Audio')
   NIRS.Dt = rmfield(NIRS.Dt,'Audio') ;
end


          
        job.NIRSmat{1} =fullfile(pathoutlist,'NIRS.mat');          
        save(fullfile(pathoutlist,'NIRS.mat'),'NIRS')
        out.NIRSmat = job.NIRSmat;


