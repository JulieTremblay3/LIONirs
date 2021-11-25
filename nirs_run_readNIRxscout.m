function out = nirs_run_readNIRxscout(job)
%prjname: .prj file for project
%fileIn: raw data file in .nirs in homer format
%fileOut: .nir, .mat, .vhdr, .vmkr

%PARAMETRE D'AQUISITION
NIRS.Cf.dev.n = 'NIRx'; %'ISS Imagent' %'NIRSx';
%Use ISS Imagent label for now naming A a1b2
%All export channel list have been implement with ISS naming and need to be
%adpated and tested for NIRx label. 

    pathout = job.output_path;   
    if isdir(pathout)
        choice = questdlg(['Directory ', pathout,' already exist do  you want to continue and overwrite this folder ?'], ...
        'WARNING',...
            'Yes Overwrite', 'Cancel','Cancel');
    switch choice
         case 'Cancel'
            return
        case 'Yes Overwrite'
            try 
                 rmdir(pathout, 's');
                 mkdir(pathout);
            catch
                msgbox(['Please remove ' , pathout,' manually'])
                return
            end
    end
    else
        mkdir(pathout)
    end

%PARAMETRE DU MONTAGE (HELMET)
NIRS.Cf.H.n = 'SainteJustine prj';
%A completer pour mettre toute info du casque, inutilisé dans la présente
%toolbox
LoadedStruct = load(job.prjfile{1},'-mat');
LoadedStruct = Load_PrjStruct(LoadedStruct,false);
DispHelm = get_Helmet( LoadedStruct );
sMtg = get_Mtg( DispHelm );
vHoles = get_vHoles( DispHelm );
NIRS.Cf.H.prj = job.prjfile{1};

%ASSOCIATION SOURCE ML 1 = SOURCE HELMET 018001


%FIDUCIAL
NIRS.Cf.H.F.r.o.mm = sMtg.matFiducials(1:3,:)';
%Source JUSTE POUR LE IMAGINC 16 CANAUX AVEC LE CASQUE STE-JUSTINE
listsrs = [018001,020003,022005,024007,026009,028011,030013,032015,017002,019004,021006,023008,025010,027012,029014,031016];
%Source POUR LE IMAGINC 32 CANAUX AVEC LE PRJ STE-JUSTINE
listsrs = [018001,020003,022005,024007,026009,028011,030013,032015,017002,019004,021006,023008,025010,027012,029014,031016,...
    050033, 052035, 054037, 056039, 058041, 060043, 062045, 064047, 049034, 051036, 053038, 055040, 057042, 059044,061046, 063048,...
    082065, 084067, 086069, 088071, 090073, 092075, 094075, 096079, 081066, 083068, 085070, 087072, 089074, 091076, 093078, 095080,...
    114097, 116099, 118101, 120103, 122105, 124107, 126109, 128111, 113098, 115100, 117102, 119104, 121106, 123108, 125110, 127112];


NIRS.Cf.H.S.N =[16]; %Nb of source
for i=1:numel(listsrs)
    p = find(listsrs(i)==sMtg.v_HolesMtg);
    if numel(p)==1
        NIRS.Cf.H.S.r.o.mm.p(:,i) =[vHoles(p).Coord.x,vHoles(p).Coord.y,vHoles(p).Coord.z]*100; %position mm
    else
        NIRS.Cf.H.S.r.o.mm.p(:,i) = [0,0,0];
    end
end
%Detector
listdet = [1:32]*1000000; %maximal num of detector
NIRS.Cf.H.D.N =[32]; %Nb of detecteur
for i=1:numel(listdet)
    p = find(listdet(i)==sMtg.v_HolesMtg);
    if numel(p)==1
        NIRS.Cf.H.D.r.o.mm.p(:,i) =[vHoles(p).Coord.x,vHoles(p).Coord.y,vHoles(p).Coord.z]*100; %position mm
    else
        NIRS.Cf.H.D.r.o.mm.p(:,i) = [0,0,0];
    end
end




%PARAMETRE
NIRS.Dt.s.age =job.age1 ;
NIRS.Dt.ana.T1 = ''; %FICHIER DE IRM

%FILE
%NIRS.Dt.fir.pp.p %read and process tdms file 
fprintf('%s\n','File open: ');
    inputrawscout = job.inputNIRxscout;
for Idx_File=1:numel(inputrawscout)

    [dir1,fil1,ext1] = fileparts(inputrawscout{Idx_File});
    fprintf('%s\n', inputrawscout{Idx_File})
    filehdr = fullfile(dir1,[fil1,ext1]);
    info = NIRxreadHDR(filehdr);
    filewav1 = fullfile(dir1,[fil1,'.wl1']);
    filewav2 = fullfile(dir1,[fil1,'.wl2']);
    fileevt = fullfile(dir1,[fil1,'.evt']);
    evt = load(fileevt);    
    WV1 = load(filewav1);
    WV2 = load(filewav2);
%      figure;hold on
%     for i=1:size(WV1,2)
%         plot(WV1(:,i),'displayname',num2str(i))
%     end
%     if isfield(job.c_shortdistance,'b_shortdistanceavailable') %detetor (nb-1)+8
%         info.nbDet = 23;
%     end 
    i = 1;
    for iwv = 1:2
        for  iSrs = 1:info.nbSrs  
            for iDet = 1:info.nbDet
                ml(i,1) = iSrs;
                ml(i,2) = iDet;
                ml(i,3) = 1;
                ml(i,4) = iwv; 
                i = i+ 1;
            end
        end
    end
    %NIRSx as Sr1 Det1, Sr1 Det2...
    DATA.ml = ml;
    
    if size(WV1,1)==size(WV2,1)
        DATA.d =  [WV1,WV2];
    else
        msgbox('WARNING nb of point for wavelength one and wavelenght 2 are different');
        maxnb = min([size(WV1,1),size(WV2,1)]);
        DATA.d = [WV1(1:maxnb,:),WV2(1:maxnb,:)];
    end
%     figure;  
%      hold on 
%      bar(mean(WV2),'r','displayname','850')
% id = find(ml(:,2)==job.c_shortdistance.b_shortdistanceavailable.e_shortdistancedet & ml(:,4)==1)
%  id = find(ml(:,2)==16 & ml(:,4)==1)
%     figure;imagesc(WV2(:,id)') 
%      figure;plot(WV2(:,id))
%      id = find(ml(:,1)==1 &ml(:,4)==1)
%      id2 = find(ml(:,1)==9&ml(:,4)==1)
%      figure;hold on
%      plot(WV2(:,id),'r')
%      plot(WV2(:,id2),'b')
     
% %     bar(mean(WV1),'b','displayname','760')
    clear WV1 WV2 ml
    NIRS.Cf.dev.wl =info.wavelenght;
    NIRS.Cf.dev.fs = info.SamplingRate;

% Handle response
    %ML
    %DISTANCE GEOMETRIQUE
    for i=1:size(DATA.ml,1)
        temp{i,1} = ['S',num2str(DATA.ml(i,1)),'_D',num2str(DATA.ml(i,2))];  %Cell channel identification
        NIRS.Cf.H.C.id(1,i) = i;
        NIRS.Cf.H.C.id(2,i) = DATA.ml(i,1); %Source
        NIRS.Cf.H.C.id(3,i) = DATA.ml(i,2);  %et detecteur
        NIRS.Cf.H.C.wl(i) = DATA.ml(i,4);
        x1 =  NIRS.Cf.H.S.r.o.mm.p(1,DATA.ml(i,1));
        x2 =  NIRS.Cf.H.D.r.o.mm.p(1,DATA.ml(i,2));
        y1 =  NIRS.Cf.H.S.r.o.mm.p(2,DATA.ml(i,1));
        y2 =  NIRS.Cf.H.D.r.o.mm.p(2,DATA.ml(i,2));
        z1 =  NIRS.Cf.H.S.r.o.mm.p(3,DATA.ml(i,1));
        z2 =  NIRS.Cf.H.D.r.o.mm.p(3,DATA.ml(i,2));
        NIRS.Cf.H.C.gp(i,1) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);%À CALCULER DISTANCE GÉOMETRIQUE
        pos(i,1) = (x2-x1)/2+x1;
        pos(i,2) = (y2-y1)/2+y1;
        pos(i,3) = (z2-z1)/2+z1; 
        pos(i,4) = NIRS.Cf.H.C.gp(i,1);
    end
    
  %  pos(id,4) = sqrt( (DetPos(1)-SrcPos(1))^2 + (DetPos(2)-SrcPos(2))^2+ (DetPos(3)-SrcPos(3))^2);           %'di
    
    Dist = NIRS.Cf.H.C.gp;
    %Enlever les canaux qui ne sont pas dans les distances
    %idbad = find(NIRS.Cf.H.C.gp > job.distmax  | NIRS.Cf.H.C.gp < job.distmin );
    %Enlever les canaux qui ne sont pas dans le mask
   
    % idbad = find([info.mask,info.mask ]==0);  
    idbad = find(NIRS.Cf.H.C.gp<job.distmin|NIRS.Cf.H.C.gp>job.distmax);
    
    %ADD SHORT DISTANCE CHANNEL
  
    if isfield(job.c_shortdistance,'b_shortdistanceavailable')
        DETSHORT = job.c_shortdistance.b_shortdistanceavailable.e_shortdistancedet;
        idbad = [idbad; find(DATA.ml(:,2)==DETSHORT)];
        SSHORT = job.c_shortdistance.b_shortdistanceavailable.e_shortdistancesrs; %[2,3,5,7,11,10,13,15];
        tokeep = [];
        for isrs=1:numel(SSHORT)
         idkeep =  find(DATA.ml(idbad,1)==SSHORT(isrs)& DATA.ml(idbad,2)== DETSHORT); %detector 16
        % Keep all detector 8 source for short distance regression whenever the
        % distance on the helmet
         tokeep = [tokeep; idkeep]; 
        end  
        idbad(tokeep) = [];
    end
    
    
     idall = ones(size(DATA.ml,1),1); 
     idall(idbad) = [];
     %idall(idbad)=0;
     
     idgood = find(idall);
     for i = 1:numel(idgood)
        NIRS.Cf.H.C.n{i,1}  = temp{idgood(i) ,1};
     end
    NIRS.Cf.H.C.N = numel(NIRS.Cf.H.C.n);
   
   %  idbad = find(std(DATA.d)./mean(DATA.d)>job.m_removestd_tdms )
    NIRS.Cf.H.C.gp(idbad) = [];
    NIRS.Cf.H.C.wl(idbad) = [];
    NIRS.Cf.H.C.id(:,idbad)=[];
    DATA.ml(idbad,:)=[];
    try
    DATA.d(:,idbad)=[];
    catch
        msgbox('Do you record short distance probe ?')
    end
    %TRIG A VERIFIER ACTUELEMENT PRENDRE LE TRIG MANUELLE DE L'ORDI
    NIRS.Dt.fir.rons = [];
    
    NIRS.Cf.H.C.ok(:,Idx_File) = ones(size(DATA.ml,1),1);
    
   % NIRS.Cf.H.C.ok(find(stdbad),Idx_File)=0;
   % NIRS.Cf.H.C.ok(find(lowbad),Idx_File)=0;
   
    pos(idbad,:) = [];
    %Trigger value
       
    if isfield(job.c_shortdistance,'b_shortdistanceavailable')
        SHORTmlid = find(DATA.ml(:,2)== DETSHORT & DATA.ml(:,4)== 1);  
        if isempty(SHORTmlid)
            disp('Please verify short distance they could not be include.')
            job.c_shortdistance = rmfield(job.c_shortdistance,'b_shortdistanceavailable');
        end
    end
    if isfield(job.c_shortdistance,'b_shortdistanceavailable')
        distancemat = zeros(size(pos,1)/2,numel(SHORTmlid ));       
        %DISTANCE BETWEEN CHANNEL FOR CLUSTERING
        for i= 1:size(pos,1)/2
            for j = 1:numel(SHORTmlid )         
             x1 = pos(i,1); 
             y1 = pos(i,2); 
             z1 = pos(i,3); 
             ishort = SHORTmlid (j);
             x2 = pos(ishort,1); 
             y2 = pos(ishort,2); 
             z2 = pos(ishort,3);          
             distancemat(i,j) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);
            end
        end
    %take the closest from the regressor to create zone each channel is
    %take only once                  
       zone.SD.Lambda = NIRS.Cf.dev.wl;
       zone.SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
       zone.SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
       zone.SD.nSrcs  = numel(zone.SD.SrcPos)/3;
       zone.SD.nDets  = numel(zone.SD.DetPos)/3;
       zone.SD.SrcAmp = ones(zone.SD.nSrcs ,2);
       zone.SD.DetAmp = ones(zone.SD.nDets,2);
       zone.color = [lines(8);lines(8)];     
       zone.ml = DATA.ml;
       zone.peak = [];
       zone.pos = pos;
       %identify the 8 first izone for each short channel position 
       [val,id] = sort(distancemat);
       %define de zone (associate to the closest distance of the short
       %detector
       ishort = 1;
       for izone = 1:8
           id(ishort, 1);
           zone.label{izone} = ['Short', num2str(DATA.ml(id(1,ishort),1))] ;
           zone.plot{izone} = [];
           zone.plotLst{izone} =  [];
           ishort = ishort+1;
       end
       %associate all other channel to the closest zone 
         for izone = 1:8
           for ich = 1:size(distancemat,1)
              [val,izone]= min(distancemat(ich,:));
              if val>0
                zone.plot{izone} = [zone.plot{izone};[DATA.ml(ich,1),DATA.ml(ich,2)]];
                zone.plotLst{izone} = [zone.plotLst{izone}, ich];
              end          
           end
         end
       ishort = 1;
       for izone = 9:1:16
           id(ishort, 1);
           zone.label{izone} = ['Regressor Short', num2str(DATA.ml(id(1,ishort),1))] ;
           zone.plot{izone} = [DATA.ml(id(1,ishort),1),DATA.ml(id(1,ishort),2)];
           zone.plotLst{izone} =  [id(1,ishort)];
           ishort = ishort+1;
       end             
        save( [pathout,filesep, 'SHORTDISTANCE.zone'],'zone');
        disp(['Save zone: ',  pathout,filesep, 'SHORTDISTANCE.zone'])
        clear zone
         distancemat = zeros(size(pos,1)/2,numel(SHORTmlid ));       
        %DISTANCE BETWEEN CHANNEL FOR CLUSTERING
        for i= 1:size(pos,1)/2
            for j = 1:numel(SHORTmlid )         
             x1 = pos(i,1); 
             y1 = pos(i,2); 
             z1 = pos(i,3); 
             ishort = SHORTmlid (j);
             x2 = pos(ishort,1); 
             y2 = pos(ishort,2); 
             z2 = pos(ishort,3);          
             distancemat(i,j) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);
            end
        end
    %take the closest from the regressor to create zone each channel is
    %take only once                  
       zone.SD.Lambda = NIRS.Cf.dev.wl;
       zone.SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
       zone.SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
       zone.SD.nSrcs  = numel(zone.SD.SrcPos)/3;
       zone.SD.nDets  = numel(zone.SD.DetPos)/3;
       zone.SD.SrcAmp = ones(zone.SD.nSrcs ,2);
       zone.SD.DetAmp = ones(zone.SD.nDets,2);
       zone.color = [lines(8);lines(8)];     
       zone.ml = DATA.ml;
       zone.peak = [];
       zone.pos = pos;
       %identify the 8 first izone for each short channel position 
       [val,id] = sort(distancemat);
       %define de zone (associate to the closest distance of the short
       %detector
         zone.label{1} = ['Regressor allSD'] ;
         zone.plot{1} = [DATA.ml(SHORTmlid ,1),DATA.ml(SHORTmlid ,2)];
         zone.plotLst{1} =  [SHORTmlid ];
           
         zone.label{2} = ['allSD'] ;
         zone.plot{2} = [DATA.ml(DATA.ml(: ,4)==1,1), DATA.ml(DATA.ml(: ,4)==1,2)];
         zone.plotLst{2} =  [find(DATA.ml(: ,4)==1)];
         save( [pathout, filesep,'AllShortDistance.zone'],'zone')
        disp(['Save zone: ',  pathout, filesep,'AllShortDistance.zone']);
              
          clear zone
    end
    
     zone.SD.Lambda = NIRS.Cf.dev.wl;
     zone.SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
     zone.SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
     zone.SD.nSrcs  = numel(zone.SD.SrcPos)/3;
     zone.SD.nDets  = numel(zone.SD.DetPos)/3;
     zone.SD.SrcAmp = ones(zone.SD.nSrcs ,2);
     zone.SD.DetAmp = ones(zone.SD.nDets,2);
     zone.color = [lines(8);lines(8)];     
     zone.ml = DATA.ml;
     zone.peak = [];
     zone.pos = pos;
    zone.label{1} = ['Regressor All'];
    zone.plot{1} = [DATA.ml(DATA.ml(: ,4)==1,1), DATA.ml(DATA.ml(: ,4)==1,2)];
    zone.plotLst{1} =  [find(DATA.ml(: ,4)==1)];
    zone.label{2} = ['All'];
    zone.plot{2} = [DATA.ml(DATA.ml(: ,4)==1,1), DATA.ml(DATA.ml(: ,4)==1,2)];
    zone.plotLst{2} =  [find(DATA.ml(: ,4)==1)];
    save( [pathout, filesep,'Global.zone'],'zone');
    disp(['Save zone: ',  pathout, filesep,'Global.zone']);
if numel(evt)>0 % Verify that the event matrix is not empty.
    NIRS.Dt.fir.aux5{Idx_File} = [ones(size(evt,1),1),evt(:,1)]; %valeur 1 internal trigger
    NIRS.Dt.fir.aux5{Idx_File} = [ones(size(evt,1),1),evt(:,1)]; %valeur 1 internal trigger
    trigdec = evt(:,2)*2^0+  evt(:,3)*2^1+evt(:,4)*2^2+evt(:,5)*2^3; % Convert the trigger identifier byte to decimal.
    NIRS.Dt.fir.aux5{Idx_File} = [trigdec ,evt(:,1)]; % Concatenate the decimal trigger identifier with the time collumn of the event file.
else
     NIRS.Dt.fir.aux5{Idx_File} = [1 1]; % Create a row of two ones if no event is provided.
end
if isfield(job.c_nameconvention_NIRxscout   ,'b_defaultname_NIRxscout')
    fileout = fil1;    
elseif isfield(job.c_nameconvention_NIRxscout  ,'b_manualname_NIRxscout')
     manualname = job.c_nameconvention_NIRxscout.b_manualname_NIRxscout.e_manualname_NIRxscout;
     outname = strsplit( manualname,',');
    if numel(inputrawscout)==numel(outname)
        fileout = outname{Idx_File};   
        disp(['Name ',fileout, ' is used to file ', fullfile(dir1,fil1), '.'])
    else 
        fileout = [outname{1}, sprintf('%03.0f',Idx_File)];
        disp(['Name ',fileout, ' is used to file ', fullfile(dir1,fil1), ' match the name and number of file open if you want to be more specific.'])
    end
elseif isfield(job.c_nameconvention_NIRxscout  ,'b_foldername_NIRxscout')
       [str,rem] = strtok( fliplr(fullfile(dir1,fil1)),'\/');
       fileout = fliplr(strtok(rem,'\/'));
        disp(['Name ',fileout, ' is used to file ', fullfile(dir1,fil1)]);       
end
    



    if ~isempty(evt)
        fid = fopen(fullfile(pathout,[fileout,'_evt1trig.m']),'w');
        fprintf(fid,'%s\n',['filename = [{''',fileout,'''}];' ]);
        fprintf(fid,'%s\n', 'Trigvalue = 1');
        if numel(evt)==1
            fprintf(fid,'%s\n', ['timingfile{1} = [',num2str(evt(i)/NIRS.Cf.dev.fs),'];']);
        else
            for i=1:numel(evt(:,1))
                if i==1 
                    fprintf(fid,'%s\n', ['timingfile{1} = [',num2str(evt(i,1)/NIRS.Cf.dev.fs),';']);
                elseif i==numel(evt(:,1))
                    fprintf(fid,'%s\n', [num2str(evt(i,1)/NIRS.Cf.dev.fs),'];']);
                else
                    fprintf(fid,'%s\n', [num2str(evt(i,1)/NIRS.Cf.dev.fs),';']); 
                end
            end 
        end
        fclose(fid);
    end
     preprocess.STDamp =  job.STD_amp_choice.STD_amp;
     preprocess.STDampenable = job.STD_amp_choice.STD_enable;  
     preprocess.STDmenu = job.STD_amp_choice.STD_menu;
     preprocess.DCamp = job.DC_amp_choice.DC_amp;
     preprocess.DCampenable = job.DC_amp_choice.DC_enable;
     fs=  NIRS.Cf.dev.fs;
     
     
    [dout,ok] = fast_preprocessing(DATA.d',  preprocess, fs);
     NIRS.Cf.H.C.ok(:,Idx_File) = ok;
     NIRS.Dt.fir.sizebloc{Idx_File} = size(DATA.d,1); 
    %WRITE OUTPUT FILES
    NIRS.Dt.s.p = [pathout,'/'];%RAW NIRS PATH    
    fileOut_nir = fullfile(pathout,[fileout,'.nir']);
    fileOutRoot_vhdr = fullfile(pathout,[fileout,'.vhdr']);
    fileOutRoot_vmrk = fullfile(pathout,[fileout,'.vmrk']);
    
    fwrite_NIR(fileOut_nir,DATA.d');
    disp(['File processed: ',fileOut_nir])
    
    nirs_boxy_write_vhdr(fileOutRoot_vhdr,... %Output file
        fileOut_nir,... %DataFile
        fileOutRoot_vmrk,... %MarkerFile,...
        'nirs_run_readNIRxscout',... %Function that created the header
        '',... %Channel Resolution
        '',... %Channel Units
        NIRS.Cf.H.C.n,... %names given as a column of cells
        1/NIRS.Cf.dev.fs*1e6,... %SamplingInterval in microseconds
        NIRS.Dt.fir.sizebloc{Idx_File}); %SamplingInterval in microseconds
    
    %dummy segment, add real marker read
    SD.Markers{1,Idx_File}.Type='New Segment';
    SD.Markers{1,Idx_File}.Description='';
    SD.Markers{1,Idx_File}.Position=1;
    SD.Markers{1,Idx_File}.Size=1;
    SD.Markers{1,Idx_File}.ChNumber=0; %all channels
    temp_markers{1,1} = SD.Markers{1,Idx_File};
    idmarker = 1;
    for istimulus=1:size(evt,1)
        SD.Markers{idmarker,1}.Type='trigger';
        label = sprintf( 'S%3.0f', bin2dec(num2str(fliplr(evt(istimulus,2:9)))));
        SD.Markers{idmarker,1}.Description=[label];
        SD.Markers{idmarker,1}.Position=evt(istimulus,1);
        SD.Markers{idmarker,1}.Size=1;
        SD.Markers{idmarker,1}.ChNumber=0; %all channels
        idmarker =  idmarker + 1;
    end
 %   temp_markers = SD.Markers{Idx_File,: };
    nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
        fileOut_nir,... %DataFile
        SD.Markers);
    
    
    %FIRST STEP OF PROCESSING PP1
    NIRS.Dt.fir.pp(1).p{Idx_File,1} = fileOut_nir;
    NIRS.Dt.fir.pp(1).pre = 'READ_RAW_NIRxScout';
    NIRS.Dt.fir.pp(1).job = job;
end         
       
        save(fullfile(pathout,'NIRS.mat'),'NIRS');
        job.NIRSmat{1} = fullfile(pathout,'NIRS.mat');
        out.NIRSmat = job.NIRSmat;
end

function info = NIRxreadHDR(file)
        fid=fopen(file);
        while ~feof(fid)
            line = fgetl(fid);
            if strfind(line,'SamplingRate')
                [tok,rem]=strtok(line,'=');
                info.SamplingRate = str2num(rem(2:end)); 
            end
            if strfind(line,'Wavelengths')
                [tok,rem]=strtok(line,'=');
                info.wavelenght = str2num(rem(3:end-1)); 
            end
            if strfind(line,'Device')
                [tok,rem]=strtok(line,'=');
                info.Device = rem(2:end); 
            end
            if strfind(line,'Sources')
                [tok,rem]=strtok(line,'=');
                info.nbSrs= str2num(rem(2:end)); 
            end
              if ~isempty(strfind(line,'Detectors'))* isempty(strfind(line,'ShortDetectors'))
                [tok,rem]=strtok(line,'=');
                info.nbDet= str2num(rem(2:end)); 
                %  info.nbDet = 23
              end
            if strfind(line,'S-D-Key=')
            1;
            end
             if strfind(line,'S-D-Mask="#')
                maskdata = [];
                 while isempty(strfind(line, '#"'))
                    line = fgetl(fid);
                    maskdata = [maskdata, str2num(line)];
                 end
                 info.mask = maskdata;
             end
        end
       fclose(fid);
end

function  dec = bi2dehome(bin)

dec = bin(1)*2^1+bin(2)*2^2+bin(3)*2^3+bin(4)*2^4;

end
