function out = nirs_run_readNIRxSPORT(job)
%prjname: .prj file for project
%fileIn: raw data file in .nirs in homer format
%fileOut: .nir, .mat, .vhdr, .vmkr

%PARAMETRE D'AQUISITION
NIRS.Cf.dev.n = 'NIRx'; %'ISS Imagent' %'NIRSx';
%Use ISS Imagent label for now naming A a1b2
%All export channel list have been implement with ISS naming and need to be
%adpated and tested for NIRx label. 

    pathout = job.output_path;   
       
    k = [strfind( pathout,'/') , strfind( pathout,'\')];
    if ~isempty(k)
        pathout(k)= filesep;
    end
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
                disp(['Please remove: ' , pathout]);
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
    inputrawscout = job.inputNIRxSport;
for Idx_File=1:numel(inputrawscout)

    [dir1,fil1,ext1] = fileparts(inputrawscout{Idx_File});
    fprintf('%s\n', inputrawscout{Idx_File})
    filehdr = fullfile(dir1,[fil1,ext1]);
    info = NIRxreadHDR(filehdr);

    %based on info about nbSrc and nbDet create a ml file 
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


    fil1 = fil1(1:end-7);
    if isfield(job.c_PruningNIRSport,'Pruning_nirsport') %default case .wvl data already containt the exact number of channel corresponding to the channel mask in the .hdr         
   
          filewav1 = fullfile(dir1,[fil1,'.wl1']);%fil1
          filewav2 = fullfile(dir1,[fil1,'.wl2']);%fil1,
          WV1 = load(filewav1);
          WV2 = load(filewav2);
          disp(['Channel mask containt ', num2str(sum(info.Channelmask(:))),' channels']  )
          disp(['Load file ', filewav1, ' with ', num2str(size(WV1,2)),' channel'])
    
          %

          %use the channel mask in for bad channel in the mlfile 
          id = 0;
          mlSPORT = [];
          idbad =  [];
          for i=1:size(info.Channelmask,1)
              for j=1:size(info.Channelmask,2)
                  id = id+1;
                  if info.Channelmask(i,j)==1
                      mlSPORT =[mlSPORT; [i,j, 1, 1]];
                  else
                      idbad = [  idbad ; id];
                  end
              end
          end
          idbad =  [idbad; idbad + (size(mlSPORT,1) + numel(idbad))];
          DATA.ml(idbad,:) = [];% adjust ml to actual already exported in wv according to nirsport config see ncfg file 
          idbad = []; %no additionnal pruning necessary
    elseif isfield(job.c_PruningNIRSport,'distmax_nirsport') 
        try
        filewav1 = fullfile([dir1,filesep],['.wl1']);%fil1
        filewav2 = fullfile([dir1,filesep],['.wl2']);%fil1,
        WV1 = load(filewav1); 
        WV2 = load(filewav2);
        disp(['Load file ', filewav1, ' with ', num2str(size(WV1,2)),' channel'])   
        catch
            disp('To prune distance you first need to use the nirs converter on the .roh or .zip file to export .wl1 and .wl2 with all recorded channel (as example 16det x 16src) = 256 ch ')
        end 
    elseif  isfield(job.c_PruningNIRSport,'channelmask_nirsport') %custom mask 
        try
        filewav1 = fullfile([dir1,filesep, fil1],['.wl1']);%fil1
        filewav2 = fullfile([dir1,filesep, fil1],['.wl2']);%fil1,
        WV1 = load(filewav1);
        WV2 = load(filewav2);
        disp(['Load file ', filewav1, ' with ', num2str(size(WV1,2)),' channel'])   
        catch
            disp('To prune using a custom mask you first need to use the nirs converter on the .roh or .zip file to export .wl1 and .wl2 with all recorded channel (as example 16det x 16src) = 256 ch ')
        end 
        Channelmask = reshape(str2num(job.c_PruningNIRSport.channelmask_nirsport),info.nbDet, info.nbSrs)
        %use the channel mask in for bad channel in the mlfile 
          id = 0;
          mlSPORT = [];
          idbad =  [];
          for i=1:size(Channelmask ,1)
              for j=1:size(Channelmask ,2)
                  id = id+1;
                  if Channelmask (i,j)==1
                      mlSPORT =[mlSPORT; [i,j, 1, 1]];
                  else
                      idbad = [  idbad ; id];
                  end
              end
          end
          idbad =  [idbad; idbad + (size(mlSPORT,1) + numel(idbad))];
        %  DATA.ml(idbad,:) = [];% adjust ml to actual already exported in wv according to nirsport config see ncfg file 
       






    end
    

          if size(WV1,1)==size(WV2,1)
              DATA.d =  [WV1,WV2];
          else
              msgbox('WARNING nb of point for wavelength one and wavelenght 2 are different');
              maxnb = min([size(WV1,1),size(WV2,1)]);
              DATA.d = [WV1(1:maxnb,:),WV2(1:maxnb,:)];
          end
          clear WV1 WV2 ml


    try
       fileevt = fullfile(dir1,[fil1,'_lsl.tri']);
       evtLSL = load(fileevt);     
       evt = evtLSL(:,2:3);
       disp(['LSL trigger file ', fileevt, ' is read with ' num2str(size(evt,1)), ' event'])
    catch
        evt = [];
       disp(['LSL trigger file ', fileevt, 'could not be read'])
    end
       
  
    NIRS.Cf.dev.wl =info.wavelenght;
    disp(['Wavelenght used: ' , num2str(    NIRS.Cf.dev.wl), 'nm' ])
    NIRS.Cf.dev.fs = info.SamplingRate;
    disp(['Sampling Rate used: ' , num2str(NIRS.Cf.dev.fs), 'Hz' ])

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
    
    
    Dist = NIRS.Cf.H.C.gp;
 
    %Enlever les canaux qui ne sont pas dans les distances tri des distance
    %maximal 
    if isfield(job.c_PruningNIRSport,'distmax_nirsport')         
        idbad = find(NIRS.Cf.H.C.gp > str2num(job.c_PruningNIRSport.distmax_nirsport) );   
        DATA.d(:,idbad)= [];
        channelmask = ones(16,16);
        channelmask(idbad(1:end/2))=0;
        disp(['Channel mask = ' , num2str(channelmask(:)')])
    elseif isfield(job.c_PruningNIRSport,'channelmask_nirsport')  
         DATA.d(:,idbad)= [];
    end
    
   
   
    
     idall = ones(size(DATA.ml,1),1); 
     idall(idbad) = [];     
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
    % try
    % %DATA.d(:,idbad)=[];
    % catch
    %     msgbox('Do you record short distance probe ?')
    % end
    %TRIG A VERIFIER ACTUELEMENT PRENDRE LE TRIG MANUELLE DE L'ORDI
    NIRS.Dt.fir.rons = [];
    
    NIRS.Cf.H.C.ok(:,Idx_File) = ones(size(DATA.ml,1),1);
    
   % NIRS.Cf.H.C.ok(find(stdbad),Idx_File)=0;
   % NIRS.Cf.H.C.ok(find(lowbad),Idx_File)=0;
   
    pos(idbad,:) = [];
    %Trigger value
       

    

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
    disp(['Save zone with all channels grouped: ',  pathout, filesep,'Global.zone']);
if numel(evt)>0 % Verify that the event matrix is not empty.
%    NIRS.Dt.fir.aux5{Idx_File} = [ones(size(evt,1),1),evt(:,1)]; %valeur 1 internal trigger
%    NIRS.Dt.fir.aux5{Idx_File} = [ones(size(evt,1),1),evt(:,1)]; %valeur 1 internal trigger
%    trigdec = evt(:,2)*2^0+  evt(:,3)*2^1+evt(:,4)*2^2+evt(:,5)*2^3; % Convert the trigger identifier byte to decimal.
    NIRS.Dt.fir.aux5{Idx_File} = [evt(:,2) ,evt(:,1)]; % TRig id and sample of the trig
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
   
 





     NIRS.Cf.H.C.ok(:,Idx_File) = 1;
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
        label = sprintf( 'S%3.0f', evt(istimulus,2));
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
    

%% ACCELEROMETER DATA
        if 1
            try
             Accdata = load(fullfile([dir1,filesep, fil1],[fil1, '.acc']));
            [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]= fopen_EEG(fileOut_nir ); %load default parameter for format
            idAUX = 1;


        %figure;plot(1/100:1/100:1/100*2000, Accdata(1:100*20,7))
        NIRS.Dt.AUX(idAUX).pp.p{size(Idx_File,1),1}=[];
        NIRS.Dt.AUX(idAUX).pp.sync_timesec{size(Idx_File,1),1}=[];
        NIRS.Dt.AUX(idAUX).label = [];
       
        %col 2 sample time equivalent in NIRS, acc freq 100Hz
        %col 7 xcoordinate
        
       

        %transfert to nirs sample rate using average value
       tic
        if 1
            fileoutAUX =    fullfile(pathout, [fileout, ' ACC',sprintf('%03.0f',Idx_File)]);
            NirsSample = Accdata(:,2);
            AccNirs = zeros(max(NirsSample),9);
            for i=1:max(NirsSample)
            AccNirs(i,:) = mean(Accdata(NirsSample==i,7:15));
            end
        end
        toc

         AUX.data = [ AccNirs] ;   
        AUX.ind_dur_ch=EEG.ind_dur_ch;
        AUX.ind_dur_ch(:,1) = EEG.ind_dur_ch(:,1);
        AUX.marker=EEG.marker;
        AUX.infoBV = EEG.infoBV; 
        AUX.infoBV.SamplingInterval = EEG.infoBV.SamplingInterval ; 
        label = 'Accelerometer';
        AUX.infoBV.DataType = 'TIMEDOMAIN';
        AUX.infoBV.DataFormat = 'BINARY';
        AUX.infoBV.DataOrientation='VECTORIZED';
        AUX.infoBV.NumberOfChannels = '9';
        AUX.infoBV.DataType = 'TIMEDOMAIN';
        AUX.infoBV.name_ele = {'x1pos','y1pos','z1pos','x2pos','y2pos','z2pos','x3pos','y3pos','z3pos'}; 

        AUX.infoBV.coor_r = [1 1 1 1 1 1 1 1 1 ];
        AUX.infoBV.coor_theta  = [90 90 90 90 90 90 90 90 90];
        AUX.infoBV.coor_phi = [45 45 45 45 45 45 45 45 45];
        AUX.infoBV.DataPoints = size(AUX.data,1 ); 


       fwrite_EEG(fileoutAUX,AUX,1,AUX.infoBV.DataPoints );
        NIRS.Dt.AUX(idAUX).pp(1,numel(NIRS.Dt.AUX(1).pp)).p{Idx_File,1}=fileoutAUX;
        NIRS.Dt.AUX(idAUX).pp(1,numel(NIRS.Dt.AUX(1).pp)).sync_timesec{Idx_File,1}=0;
        NIRS.Dt.AUX(idAUX).label = label;

        disp('Verify with nirx which accelerometer is used, ajust to fNIRS time frequency ')
            catch
                disp('No auxilary load')
            end
        end



    
    %FIRST STEP OF PROCESSING PP1
    NIRS.Dt.fir.pp(1).p{Idx_File,1} = fileOut_nir;
    NIRS.Dt.fir.pp(1).pre = 'READ_RAW_NIRSport';
    NIRS.Dt.fir.pp(1).job = job;
end         
       
        save(fullfile(pathout,'NIRS.mat'),'NIRS');
        job.NIRSmat{1} = fullfile(pathout,'NIRS.mat');
        out.NIRSmat = job.NIRSmat;
end

function info = NIRxreadHDR(file)
        fid=fopen(file);
        if fid==-1
            disp(['Could not open: ',file])
        end
        while ~feof(fid)
            line = fgetl(fid);
            if strfind(line,'Sampling rate')
                [tok,rem]=strtok(line,'=');
                info.SamplingRate = str2num(rem(2:end)); 
            end
            if strfind(line,'Wavelengths') 
                [tok,rem]=strtok(line,'=');
               
            end
            if strfind(line,'Device')
                [tok,rem]=strtok(line,'=');
                info.Device = rem(2:end); 
            end
            if strfind(line,'Sources')
                [tok,rem]=strtok(line,'=');
                info.nbSrs= str2num(rem(2:end)); 
            end
            if strfind(line,'Channel Mask')
                Channelmask = [];
                itmax = 1;
                 line =  fgetl(fid);
                while isempty(strfind(line, '#')) & itmax<64                 
                    Channelmask = [Channelmask; str2num(line)];
                    itmax = itmax + 1;
                      line =  fgetl(fid);
                end
                info.Channelmask = Channelmask; 
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
       info.wavelenght = [760,850];   %disapear in NIRSPORT .hdr format DEVICE  use value [760,850] from manufacturer 
       fclose(fid);
end

function  dec = bi2dehome(bin)

dec = bin(1)*2^1+bin(2)*2^2+bin(3)*2^3+bin(4)*2^4;

end
