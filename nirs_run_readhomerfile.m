function out = nirs_run_readhomerfile(job)
%Usage: generate NIRS.mat for homer .nirs file,
%Input files
%prjname: .prj file for project
%fileIn: raw data file in .nirs in homer format
%fileOut: .nir, .mat, .vhdr, .vmkr

%PARAMETRE D'AQUISITION 
NIRS.Cf.dev.n = 'NIRS FILE HOMER';



%PARAMETRE DU MONTAGE (HELMET)
NIRS.Cf.H.n = 'SainteJustine prj';
%A completer pour mettre toute info du casque, inutilisé dans la présente toolbox

LoadedStruct = load(job.prjfile{1},'-mat');
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
    
LoadedStruct = Load_PrjStruct(LoadedStruct,false);
DispHelm = get_Helmet( LoadedStruct );
sMtg = get_Mtg( DispHelm );
vHoles = get_vHoles( DispHelm );

%ASSOCIATION SOURCE ML 1 = SOURCE HELMET 018001
Nb_Sources_Disp = sMtg.Gen_Params.Nb_Sources_Disp/2;
Nb_Detect_Disp = sMtg.Gen_Params.Nb_Detect_Disp;
%FIDUCIAL
NIRS.Cf.H.F.r.o.mm = sMtg.matFiducials(1:3,:)';
%Source JUSTE POUR LE IMAGINC 16 CANAUX AVEC LE CASQUE STE-JUSTINE
listsrs = [018001,020003,022005,024007,026009,028011,030013,032015,017002,019004,021006,023008,025010,027012,029014,031016,
    050033,052035,054037,056039,058041,060043,062045,064047, 049034,051036,053038,055040,057042,059044,061046,063048
82065,84067,86069,88071,90073,92075,94077,96079,81066,83068,85070,87072,89074,91076,93078,95080,
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

%FILE 
if job.m_inputrawhomer == 0 %write each file in there original form
    
fprintf('%s\n','File processed');
for Idx_File=1:numel(job.inputrawhomer)
     [dir1,fil1,ext1] = fileparts(job.inputrawhomer{Idx_File});
  
     
     NIRS.Dt.s.p = [dir1,'/'];%RAW NIRS PATH
     fileOut_nir = fullfile(dirout ,[fil1,'.nir']);
     fileOutRoot_vhdr = fullfile(dirout ,[fil1,'.vhdr']);
     fileOutRoot_vmrk = fullfile(dirout ,[fil1,'.vmrk']);
 
     DATA = load(fullfile(dir1,[fil1,ext1]),'-mat');
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
        NIRS.Cf.H.C.gp(i,1) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);%À CALCULER DISTANCE GÉOMETRIQUE
        NIRS.Cf.H.C.ok(i,Idx_File) = 1;
    end
    NIRS.Cf.H.C.N = numel(NIRS.Cf.H.C.n);
    NIRS.Cf.H.C.ok = ones(NIRS.Cf.H.C.N,numel(job.inputrawhomer)); 
     %STRUCTURE NON CREER 
    NIRS.Dt.fir.rons = [];     
    if isfield(DATA,'aux')
        if ~isempty(DATA.aux)        
            ind =  find(DATA.aux(:,1)>1.5);   
        end
    end
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
                1/NIRS.Cf.dev.fs,... %SamplingInterval in microseconds
                NIRS.Dt.fir.sizebloc); %SamplingInterval in microseconds
            
     %Create a New Segment marker
    SD.Markers{Idx_File,1}.Type='New Segment';
    SD.Markers{Idx_File,1}.Description='';
    SD.Markers{Idx_File,1}.Position=1;
    SD.Markers{Idx_File,1}.Size=1;
    SD.Markers{Idx_File,1}.ChNumber=0; %all channels       
    temp_markers{1,1} = SD.Markers{Idx_File,1};
    nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
                fileOut_nir,... %DataFile
                temp_markers);
    %FIRST STEP OF PROCESSING PP1
    NIRS.Dt.fir.pp(1).p{Idx_File,1} = fileOut_nir; 
    NIRS.Dt.fir.pp(1).pre = 'READ_RAW_NIRS';
    NIRS.Dt.fir.pp(1).job = job;   
end
elseif job.m_inputrawhomer ==1 || job.m_inputrawhomer ==2 %merge and detrend file
   DATA.d =[];
   DATA.aux = [];
   DATA.s = [];
   %loop all file
    [dir1,fil1,ext1] = fileparts(job.inputrawhomer{1});  
    NIRS.Dt.s.p = [dirout];%RAW NIRS PATH
    fileOut_nir = fullfile(dirout ,[fil1,'.nir']);
    fileOutRoot_vhdr = fullfile(dirout ,[fil1,'.vhdr']);
    fileOutRoot_vmrk = fullfile(dirout ,[fil1,'.vmrk']);

    
    for Idx_File=1:numel(job.inputrawhomer)
        [dir1,fil1,ext1] = fileparts(job.inputrawhomer{Idx_File});       
       
         tmp = load(fullfile(dir1,[fil1,ext1]),'-mat');
        [val,id] = sort(tmp.ml(:,4));
        DATA.ml = tmp.ml(id,:); %Wavelenght 1 et wavelength 2
        intensnorm = tmp.d(:,id);
        if  job.m_inputrawhomer ==2 %detrend here. 
            X = 1:1:size(intensnorm,1);
            Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
            Mb2 =  intensnorm(1,:)'; %offset
            A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
         intensnorm = intensnorm - A;
        figure;subplot(3,1,1);plot(tmp.d(:,id));subplot(3,1,2);plot(A);subplot(3,1,3);plot(intensnorm);

    end 
    
        tmp.s(1) = 255; %add trig 255 on each segment
        DATA.d = [DATA.d; intensnorm];
        DATA.aux = [DATA.aux, tmp.aux];
        DATA.s =  [DATA.s; tmp.s] ;
        NIRS.Cf.dev.wl = tmp.SD.Lambda;
        NIRS.Cf.dev.fs = 1/(tmp.t(2)-tmp.t(1));

    end
    ind =  find(DATA.s(:,1)==255);
     NIRS.Dt.fir.aux5{1} = [ones(size(ind,1),1)*255,ind];
    %write new file 
    
    fwrite_NIR(fileOut_nir,DATA.d');
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
        NIRS.Cf.H.C.gp(i,1) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);%À CALCULER DISTANCE GÉOMETRIQUE
        NIRS.Cf.H.C.ok(i,Idx_File) = 1;
    end
    NIRS.Dt.fir.rons = [];
    if isfield(DATA,'aux')
        if ~isempty(DATA.aux)        
            ind =  find(DATA.aux(:,1)>1.5);   
        end
    end
    %NIRS.Dt.fir.aux5{Idx_File} = [ones(size(ind,1),1)*254,ind];
    NIRS.Dt.fir.sizebloc = size(DATA.d,1);
    NIRS.Cf.H.C.N = numel(NIRS.Cf.H.C.n);
    NIRS.Cf.H.C.ok = ones( NIRS.Cf.H.C.N ,1); 
    nirs_boxy_write_vhdr(fileOutRoot_vhdr,... %Output file
                fileOut_nir,... %DataFile
                fileOutRoot_vmrk,... %MarkerFile,...
                'nirs_convert_boxy',... %Function that created the header
                '',... %Channel Resolution
                '',... %Channel Units
                NIRS.Cf.H.C.n,... %names given as a column of cells 
                1/NIRS.Cf.dev.fs,... %SamplingInterval in microseconds
                NIRS.Dt.fir.sizebloc); %SamplingInterval in microseconds
            
     %Create a New Segment marker
    SD.Markers{Idx_File,1}.Type='New Segment';
    SD.Markers{Idx_File,1}.Description='';
    SD.Markers{Idx_File,1}.Position=1;
    SD.Markers{Idx_File,1}.Size=1;
    SD.Markers{Idx_File,1}.ChNumber=0; %all channels       
    temp_markers{1,1} = SD.Markers{Idx_File,1};
    nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
                fileOut_nir,... %DataFile
                temp_markers);
    %FIRST STEP OF PROCESSING PP1
    NIRS.Dt.fir.pp(1).p{1} = fileOut_nir; 
    NIRS.Dt.fir.pp(1).pre = 'READ_RAW_NIRS';
    NIRS.Dt.fir.pp(1).job = job;    

end
    save(fullfile(dirout,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} =fullfile(dirout,'NIRS.mat');



out.NIRSmat = job.NIRSmat;
