function out = nirs_run_readhomerfile(job)
%Usage: generate NIRS.mat for homer .nirs file,
%Input files
%prjname: .prj file for project
%fileIn: raw data file in .nirs in homer format
%fileOut: .nir, .mat, .vhdr, .vmkr

%PARAMETRE D'AQUISITION 
NIRS.Cf.dev.n = 'nirs';



%PARAMETRE DU MONTAGE (HELMET)
NIRS.Cf.H.n = 'SainteJustine prj';


%Verify job.NIRS.mat location 
 %dir2 will be the new location 
dirout = job.output_path ; 
    if isdir(dirout)
        choice = questdlg(['Directory ', dirout,' already exist do  you want to continue and overwrite this folder ?'], ...
        'WARNING',...
            'Yes Overwrite', 'Cancel','Cancel');
    switch choice
         case 'Cancel'
            return
        case 'Yes Overwrite'
            try 
                 rmdir(dirout, 's');
                 mkdir(dirout);
            catch
                msgbox(['Please remove ' , dirout,' manually'])
                return
            end
    end
    else
        mkdir(dirout);
    end
    
if isfield(job.c_createImportProjectSnirf,'b_createProject')
    LoadedStruct = Createproject_fromsnirf(job.inputrawhomer{1});
    SaveStruct = Create_PrjStruct( LoadedStruct );
    [filepath,name,ext] =fileparts(job.inputrawhomer{1});
    save(fullfile(dirout,[name,'.prj']), 'SaveStruct');
    prjfile = fullfile(dirout,[name,'.prj']);
    disp(['Create: ', prjfile]);
elseif isfield(job.c_createImportProjectSnirf,'b_importProject')
    prjfile = job.c_createImportProjectSnirf.b_importProject.prjfile{1};
    LoadedStruct = load( prjfile,'-mat');
    LoadedStruct = Load_PrjStruct(LoadedStruct,false);
end
DispHelm = get_Helmet( LoadedStruct );
sMtg = get_Mtg( DispHelm );
vHoles = get_vHoles( DispHelm );

%ASSOCIATION SOURCE ML 1 = SOURCE HELMET 018001
Nb_Sources_Disp = sMtg.Gen_Params.Nb_Sources_Disp/2;
Nb_Detect_Disp = sMtg.Gen_Params.Nb_Detect_Disp;
%FIDUCIAL
NIRS.Cf.H.F.r.o.mm = sMtg.matFiducials(1:3,:)';
%Source JUSTE POUR LE IMAGINC 16 CANAUX AVEC LE CASQUE STE-JUSTINE
listsrs = [018001,020003,022005,024007,026009,028011,030013,032015,017002,019004,021006,023008,025010,027012,029014,031016,...
    050033,052035,054037,056039,058041,060043,062045,064047, 049034,051036,053038,055040,057042,059044,061046,063048,...
82065,84067,86069,88071,90073,92075,94077,96079,81066,83068,85070,87072,89074,91076,93078,95080,...
114097,116099,118101,120103,122105,124107,126109,128111,113098,115100,117102,119104,121106,123108,125110,127112];

NIRS.Cf.H.S.N =[Nb_Sources_Disp]; %Nb of source
for i=1:Nb_Sources_Disp
    p = find(listsrs(i)==sMtg.v_HolesMtg);
    if numel(p)==1
        NIRS.Cf.H.S.r.o.mm.p(:,i) =[vHoles(p).Coord.x,vHoles(p).Coord.y,vHoles(p).Coord.z]*100; %position cm
    else
         NIRS.Cf.H.S.r.o.mm.p(:,i) = [0,0,0];
    end    
end
%Detector
listdet = [1:Nb_Detect_Disp]*1000000;
NIRS.Cf.H.D.N =[Nb_Detect_Disp]; %Nb of detecteur
for i=1:numel(listdet)
    p = find(listdet(i)==sMtg.v_HolesMtg);
    if numel(p)==1
        NIRS.Cf.H.D.r.o.mm.p(:,i) =[vHoles(p).Coord.x,vHoles(p).Coord.y,vHoles(p).Coord.z]*100; %position cm
    else
         NIRS.Cf.H.D.r.o.mm.p(:,i) = [0,0,0];
    end
end
 



%PARAMETRE  
NIRS.Dt.s.age = job.age1;
NIRS.Dt.ana.T1 = ''; %FICHIER DE IRM
NIRS.Cf.H.prj = prjfile;


    
fprintf('%s\n','File processed');
for Idx_File=1:numel(job.inputrawhomer)
     [dir1,fil1,ext1] = fileparts(job.inputrawhomer{Idx_File});
     disp(['Read file: ' job.inputrawhomer{Idx_File}])
     
     NIRS.Dt.s.p = [dir1,'/'];%RAW NIRS PATH
     fileOut_nir = fullfile(dirout ,[fil1,'.nir']);
     fileOutRoot_vhdr = fullfile(dirout ,[fil1,'.vhdr']);
     fileOutRoot_vmrk = fullfile(dirout ,[fil1,'.vmrk']);
 
     DATA = load(fullfile(dir1,[fil1,ext1]),'-mat');
   
    if isfield(DATA.SD,'MeasList')    
        DATA.ml = DATA.SD.MeasList;
    elseif isfield(DATA,'ml')   
        DATA.ml = DATA.ml;
    end
    [val,id] = sort(DATA.ml(:,4));     
    DATA.ml = DATA.ml(id,:); %Wavelenght 1 et wavelength 2
    DATA.d = DATA.d(:,id);
    NIRS.Cf.dev.wl = DATA.SD.Lambda;
    NIRS.Cf.dev.fs = 1/(DATA.t(2)-DATA.t(1));

    
    %ML 
    %DISTANCE GEOMETRIQUE
    for i=1:size(DATA.ml,1)        
        NIRS.Cf.H.C.n{i,1} = ['S',num2str(DATA.ml(i,1)),'_D',num2str(DATA.ml(i,2)) ];  %Cell channel identification 
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
        NIRS.Cf.H.C.gp(i,1) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);%� CALCULER DISTANCE G�OMETRIQUE
        NIRS.Cf.H.C.ok(i,Idx_File) = 1;
    end
    NIRS.Cf.H.C.N = numel(NIRS.Cf.H.C.n);
    NIRS.Cf.H.C.ok = ones(NIRS.Cf.H.C.N,numel(job.inputrawhomer)); 
   %look to trig in the S structure 
    %add a trig for each marker
    % the structure s indicate the event and is format is s Ntime x N
    % condition each different condition will be formated as a Trigger info [trig value, index of time sample]   
    %Example; [2, 9;
    %2, 1733;
    %1,1]
    %Trig 2 sample of time 2 et 1733, 
    %Trig 1 sample 
    aux5 = [];
   if isfield(DATA, 's')
       for i=1:size(DATA.s,2)
          timesample= find(DATA.s(:,i));
          for js = 1:numel(timesample)
            aux5 = [aux5;i,timesample(js)];
          end
       end
   else
       1;
   end
   if isempty(aux5)
           aux5  =[1 1]; %at least have one trig at the begining of the file
   end
   NIRS.Dt.fir.aux5{1} = aux5;
    NIRS.Dt.fir.rons = [];     
%     if isfield(DATA,'aux')
%         if ~isempty(DATA.aux)        
%             ind =  find(DATA.aux(:,1)>1.5);   
%         end
%     end
%    NIRS.Dt.fir.aux5{Idx_File} = [ones(size(ind),1)*254,ind];
    NIRS.Dt.fir.sizebloc = size(DATA.d,1);       
    fwrite_NIR(fileOut_nir,DATA.d');
  
    NIRS.Cf.H.C.wl;
    
    nirs_boxy_write_vhdr(fileOutRoot_vhdr,... %Output file
                fileOut_nir,... %DataFile
                fileOutRoot_vmrk,... %MarkerFile,...
                'nirs_convert_boxy',... %Function that created the header
                '',... %Channel Resolution
                '',... %Channel Units
                NIRS.Cf.H.C.n,... %names given as a column of cells 
                1/NIRS.Cf.dev.fs * 1000000,... %SamplingInterval in microseconds
                NIRS.Dt.fir.sizebloc); %SamplingInterval in microseconds
      if isfield(DATA,'aux')

        fileoutAUX = fullfile(dirout ,['AUX', fil1 ,'.dat']);
        idAUX = 1;
        NIRS.Dt.AUX(idAUX).pp.p{size(Idx_File,1),1}=[]; %initialised 
        NIRS.Dt.AUX(idAUX).pp.sync_timesec{size(Idx_File,1),1}=[];
        NIRS.Dt.AUX(idAUX).label = [];

         nbAUX = numel(DATA.aux)./numel(DATA.t);
         disp(['Read AUX find ', num2str(nbAUX),' channel']);
        
         reshape(DATA.aux,numel(DATA.t), nbAUX);
         AUX.data =  reshape(DATA.aux,numel(DATA.t), nbAUX);
         AUX.ind_dur_ch=[1 1 0];
         AUX.marker={'Time 0', 'S 1'};
         name_ele = [];
         for i=1:nbAUX
            name_ele = [name_ele, {['AUX', num2str(i)]}];
        
         end
        AUX.infoBV.DataFormat = 'BINARY';
        AUX.infoBV.DataOrientation='VECTORIZED';
        AUX.infoBV.DataType = 'TIMEDOMAIN';
        AUX.infoBV.NumberOfChannels = num2str(nbAUX);
        AUX.infoBV.DataPoints = numel(DATA.t);
        AUX.infoBV.SamplingInterval = round(1/NIRS.Cf.dev.fs*1000000);
        AUX.infoBV.BinaryFormat = 'IEEE_FLOAT_32';
        AUX.infoBV.name_ele = name_ele; 
         AUX.infoBV.coor_r = ones(1,nbAUX);
         AUX.infoBV.coor_theta  = ones(1,nbAUX)*90;
         AUX.infoBV.coor_phi =  ones(1,nbAUX)*45;
         
        
   
        
        fwrite_EEG(fileoutAUX,AUX,1,AUX.infoBV.DataPoints );
        NIRS.Dt.AUX(idAUX).pp(1,numel(NIRS.Dt.AUX(1).pp)).p{Idx_File,1}=fileoutAUX;
        NIRS.Dt.AUX(idAUX).pp(1,numel(NIRS.Dt.AUX(1).pp)).sync_timesec{Idx_File,1}=0;
        NIRS.Dt.AUX(idAUX).label = 'AUX';

    end






     %Create a New Segment marker
    SD.Markers{Idx_File,1}.Type='New Segment';
    SD.Markers{Idx_File,1}.Description='';
    SD.Markers{Idx_File,1}.Position=1;
    SD.Markers{Idx_File,1}.Size=1;
    SD.Markers{Idx_File,1}.ChNumber=0; %all channels       
   
    
    idmarker = 2;    
    for istimulus=1:size(aux5,1)
        SD.Markers{idmarker,1}.Type='trigger';
        SD.Markers{idmarker,1}.Description=['S ' num2str(aux5(istimulus,1))];
        SD.Markers{idmarker,1}.Position= aux5(istimulus,2);
        SD.Markers{idmarker,1}.Size=1;
        SD.Markers{idmarker,1}.ChNumber=0; %all channels
        idmarker =  idmarker + 1;
    end

    nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
                fileOut_nir,... %DataFile
                SD.Markers);
    %FIRST STEP OF PROCESSING PP1
    NIRS.Dt.fir.pp(1).p{Idx_File,1} = fileOut_nir; 
    disp(['Save: ', fileOut_nir])
    NIRS.Dt.fir.pp(1).pre = 'READ_RAW_NIRS';
    NIRS.Dt.fir.pp(1).job = job;   
end

   %DISTANCE GEOMETRIQUE
   %channel in measure list format
   ml = ones(numel(NIRS.Cf.H.C.id(2,:)),4);
   ml(:,1) =  NIRS.Cf.H.C.id(2,:); %= DATA.ml(i,1); %Source
   ml(:,2) =  NIRS.Cf.H.C.id(3,:); %= DATA.ml(i,2);  %et detecteur
   ml(:,3) = 1;
   ml(:,4) =NIRS.Cf.H.C.wl(:); %= DATA.ml(i,4);
   
    for i=1:size(ml,1)
        x1 =  NIRS.Cf.H.S.r.o.mm.p(1,ml(i,1));
        x2 =  NIRS.Cf.H.D.r.o.mm.p(1,ml(i,2));
        y1 =  NIRS.Cf.H.S.r.o.mm.p(2,ml(i,1));
        y2 =  NIRS.Cf.H.D.r.o.mm.p(2,ml(i,2));
        z1 =  NIRS.Cf.H.S.r.o.mm.p(3,ml(i,1));
        z2 =  NIRS.Cf.H.D.r.o.mm.p(3,ml(i,2));
        NIRS.Cf.H.C.gp(i,1) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);%� CALCULER DISTANCE G�OMETRIQUE
        pos(i,1) = (x2-x1)/2+x1;
        pos(i,2) = (y2-y1)/2+y1;
        pos(i,3) = (z2-z1)/2+z1; 
        pos(i,4) = NIRS.Cf.H.C.gp(i,1);
    end
  
    %idbad = find(NIRS.Cf.H.C.gp<job.distmin|NIRS.Cf.H.C.gp>job.distmax);
     zone.SD.Lambda = NIRS.Cf.dev.wl;
     zone.SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
     zone.SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
     zone.SD.nSrcs  = numel(zone.SD.SrcPos)/3;
     zone.SD.nDets  = numel(zone.SD.DetPos)/3;
     zone.SD.SrcAmp = ones(zone.SD.nSrcs ,2);
     zone.SD.DetAmp = ones(zone.SD.nDets,2);
     zone.color = [lines(8);lines(8)];     
     zone.ml = ml;
     zone.peak = [];
     zone.pos = pos;
    zone.label{1} = ['Regressor All'];
    zone.plot{1} = [ml(ml(: ,4)==1,1), ml(ml(: ,4)==1,2)];
    zone.plotLst{1} =  [find(ml(: ,4)==1)];
    zone.label{2} = ['All'];
    zone.plot{2} = [ml(ml(: ,4)==1,1), ml(ml(: ,4)==1,2)];
    zone.plotLst{2} =  [find(ml(: ,4)==1)];
    save( [dirout, filesep,'Global.zone'],'zone');
    disp(['Save: ', fullfile(dirout,'Global.zone')])
    save(fullfile(dirout,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} =fullfile(dirout,'NIRS.mat');



out.NIRSmat = job.NIRSmat;
