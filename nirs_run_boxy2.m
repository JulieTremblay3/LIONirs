function out = nirs_run_boxy2(job)
%generate NIRS.mat info and convert Boxy files to .nir binary files, 
%Input files
%prjname: .prj file for project
%fileIn: all raw Boxy data

outNIRSmat = {};
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
  

fprintf('%s\n','File processed');

  
    SD.pathoutput = pathout;
    SD.prjname =  job.prjfile{1}; 
%      NIRS = job.cf1;
%     if NIRS.resample > 1
%         NIRS.freq = NIRS.freq/NIRS.resample;
%     end
%     %Extract info from chosen files
%     %Number of BOXY data files for this subject
%     NIRS.N_dfiles = size(job.inputBOXY,1);
%     %NIRS.N_dfiles = size(job.fnames,1);
%     size(job.inputBOXY,1)
%     %path and expname
%     temp_fname = job.subj(1,Idx_subj).fnames(1,:);
%     NIRS.BOXYfiles = temp_fname;
%     %temp_fname = job.fnames(1,:);
%     [subj_BOXY_path expname] = fileparts(temp_fname{1,1});
%     %find subject root directory
%     temp_idx = strfind(subj_BOXY_path,'\');
%     NIRS.subj_path = [subj_BOXY_path(1:temp_idx(end)-1) '\'];
   
    %load project montage
    try
        LS = load(job.prjfile{1},'-mat');          
        %assign shorter names to project structures
        mtg = LS.SaveStruct.m_Helmet.Mtg_Data;
        holes = LS.SaveStruct.m_Helmet.v_Holes;
        %coordinates of fiducials, converted to centimeters
        mat_Fid = 100*LS.SaveStruct.m_DigSubjFids.matFiducials; 
    catch
        msgbox('Problem loading HSJ montage file. Check .prj file is in correct place. Aborting.');
        return
    end
    clear LS %get rid of this big project structure

    %find real world x,y,z coordinates of sources, detectors & electrodes; 
    %sources followed by detectors are in mat_Mtg
    if ~isfield(mtg,'v_pEle')
        mtg.v_pEle = []
    end
    if ~isfield(mtg,'v_HolesEle')
        mtg.v_HolesEle = [];
    end

    [mat_Mtg mat_Ele] = nirs_boxy_get_montage_coordinates(mtg.v_pSrc,...
                            mtg.v_pDet,mtg.v_pEle,holes,mtg.v_HolesMtg);

    %Add to NIRS file describing montage, to be used by boxy_convert
    %Complex code very specific to BOXY ISS Imagent
    %convert distmin and distmax to meters, as they will be converted back 
    %to centimeters later on.
    %ISS Imagent 4 in parallel 
    %WITH 128 source 
    %And 16 detector 
     nb_Det = 16;
     nb_MaxSources = 128;
     SD.Lambda = [830,690];
    SD = nirs_boxy_associate_sources_detectors(SD,job.distmax/100,job.distmin/100,mat_Mtg,...
           nb_Det,mat_Ele,mtg.v_pEle,mtg.v_HolesEle,nb_MaxSources);

    %no need anymore for mtg, mat_Mtg, mat_Ele, holes
    clear mtg mat_Mtg mat_Ele holes

    %if electrodes are not specified in the montage, load a standard reference,
    %scaled to the fiducials... not straightforward

    %Continue building the channels into structure NIRS 
    %Approximate positions of channels by 10-20 or 10-10 electrode system
   
%         SD = nirs_boxy_approx10_20or10_10(SD,1);
%         SD = nirs_boxy_extend_channel_names_lengths(SD);
        
    
    
    %Extend channel names and lengths to channels of both wavelengths
   
    
    %add useful info to NIRS.mat: frequency, sampling interval
    %NIRS.SamplingFrequency = NIRS.freq; %in Hertz
    %NIRS.SamplingInterval =floor(1000000/NIRS.freq); %in microseconds

    %Big Loop over each of the BOXY data files    
    for Idx_file = 1:numel(job.inputBOXY)
        %Output file name
        temp_fName = job.inputBOXY(Idx_file,:);
        fName = temp_fName{1,1};
        [path1 expname1 ext1] = fileparts(fName); %#ok<ASGLU>
        ext1 = ext1(2:end); %skip the initial dot "."
        shortfileOutRoot = [expname1 '_' ext1];
        fileOutRoot = fullfile(SD.pathoutput, shortfileOutRoot);
        fileOut=[fileOutRoot '.nir'];
        fileOutRoot_vhdr = [fileOutRoot '.vhdr'];
        fileOutRoot_vmrk = [fileOutRoot '.vmrk'];
        fileOut_nir = [shortfileOutRoot '.nir'];
        fileOut_vmrk = [shortfileOutRoot '.vmrk'];
        
%         %fast preprocessing option
        preprocess.STDamp =  job.STD_amp_choice.STD_amp;
        preprocess.STDampenable = job.STD_amp_choice.STD_enable;  
        preprocess.STDmenu = job.STD_amp_choice.STD_menu;
        preprocess.DCamp = job.DC_amp_choice.DC_amp;
        preprocess.DCampenable = job.DC_amp_choice.DC_enable;
         SD.resample = 1;
        [SD,ok] = nirs_boxy_convert2(SD,fName,fileOut,Idx_file,preprocess);
        SD.ok(:,Idx_file) = ok;        
        SD.SamplingInterval=floor(1000000/SD.fs);
 
%         NIRS.VMRKfile{Idx_file} = fileOutRoot_vmrk;
        if Idx_file==1 %Juste au premier fichier
            for id=1:size(SD.ml,1) %On ajoute la nomenclature Ste-Justine au nom
                srs = SDPairs2strboxy_ISS(SD.ml(id,1));
                det = SDDet2strboxy_ISS(SD.ml(id,2));
                SD.ChannelLabels{id,1} = [srs, '_', det];
            end
        end
        %Header file
        nirs_boxy_write_vhdr(fileOutRoot_vhdr,... %Output file
                fileOut_nir,... %DataFile
                fileOut_vmrk,... %MarkerFile,...
                'nirs_convert_boxy',... %Function that created the header
                '',... %Channel Resolution
                '',... %Channel Units
                SD.ChannelLabels,... %names given as a column of cells 
                SD.SamplingInterval,...
                SD.sizebloc); %SamplingInterval in microseconds
        
        %Marker file
        temp_markers{1,1} =SD.Markers{Idx_file,1};
        nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
                fileOut_nir,... %DataFile
                temp_markers);

            %Write trigs in the .vmrk file
            try
                aux_ind = find(SD.aux5);                  
                desc = full(SD.aux5(aux_ind));              
                if isempty(NIRS.aux5)%Force un trig 255 au début de l'enregistrement pour tj avoir un trig même si enregistrement n'en contient pas. 
                    desc = 255;
                    aux_ind = 1;
                end
                if 1
                         new_desc = desc(1);
                %change : 
                %aux_ind = new_desc
                aux_ind = [0,aux_ind];
                desc = [0,desc];
                aux_ind = [aux_ind];
                desc = [desc];
                descchange = desc(2:end)-desc(1:end-1);
                auxdt = aux_ind(2:end) - aux_ind(1:end-1);
                aux_ind2 = aux_ind(1);
                for i = 1:numel(desc)-1
                    if auxdt(i)  ~= 1 | descchange(i)~= 0
                        new_desc = [new_desc; desc(i+1)];
                        aux_ind2 = [aux_ind2; aux_ind(i+1)];
                    end
                end
                if Idx_file == 1
                    NIRS.auxtrig = [];
                end
                NIRS.auxtrig{Idx_file} = [new_desc, aux_ind2];
                sizebloc{Idx_file} = NIRS.sizebloc;
                new_desc = cellstr(num2str(new_desc));

                ind_dur_ch = zeros(numel(aux_ind2),3);
                ind_dur_ch(:,1) = aux_ind2;

                write_vmrk(fileOutRoot_vmrk,'trigger',new_desc,ind_dur_ch);
                end
                
                if 1                     
                    if isfield(SD,'aux1')& isfield(SD,'aux2')& isfield(SD,'aux3')& isfield(SD,'aux4')& isfield(SD,'aux5')
                        NIRS.Dt.fir.aux1{Idx_file}=SD.aux1;
                        NIRS.Dt.fir.aux2{Idx_file}=SD.aux2;
                        NIRS.Dt.fir.aux3{Idx_file}= SD.aux3;
                        NIRS.Dt.fir.aux4{Idx_file}=SD.aux4; 
                        NIRS.Dt.fir.aux5i{Idx_file}=SD.aux5i;                        
                    else
                      sprintf('The auxiliairy field could not be read set AUX to false')
                    end
                end
                
             catch
                disp(['No triggers found for file ', num2str(Idx_file)]);
            end
            
           fprintf('%s\n',fileOut_nir);
  
    
    
    
     SD.FidPos = mat_Fid;
    SD.n_Fid = size(mat_Fid,1);
    SD.n_Src = size(SD.SrcPos,1);
    SD.n_Det = size(SD.DetPos,1);
    SD.n_Opt = SD.n_Src+ SD.n_Det;
    
    %Fill out standardized NIRS structure
    NIRS = [];
    NIRS.Cf.dev.n = 'ISS Imagent';
    NIRS.Cf.dev.wl = SD.Lambda;
    NIRS.Cf.dev.fs = SD.fs;   
    NIRS.Dt.s.age = job.age1;
    
    %Data
    NIRS.Dt.fir.pp(1).pre = 'readBOXY';
    NIRS.Dt.fir.pp(1).job = job;
    NIRS.Dt.fir.pp(1).sizebloc{Idx_file} = SD.sizebloc;
    NIRS.Dt.fir.pp(1).p{Idx_file}=fileOut;
   
    NIRS.Cf.dev.wl = SD.Lambda;
    NIRS.Cf.dev.fs = SD.fs; 
    %Information on Helmet
    NIRS.Cf.H.n = 'Sainte-Justine'; %HOLE NAMING SOURCE a1b2, detector A
    NIRS.Cf.H.p =  job.prjfile;
    %Information on stimuli triggers
    if isfield(SD,'auxtrig')
        NIRS.Dt.fir.pp(1).aux5 = SD.auxtrig;
    end
    %SAVE IN BV FORMAT AND ADD TO AUX INFO 
    if 0 %job.cf1.save_aux
        NIRS.Dt.fir.aux1 = oldNIRS.aux1;
        NIRS.Dt.fir.aux2 = oldNIRS.aux2;
        NIRS.Dt.fir.aux3 = oldNIRS.aux3;
        NIRS.Dt.fir.aux4 = oldNIRS.aux4;
        NIRS.Dt.fir.aux5i = oldNIRS.aux5i;
    end
    %Information on fiducials
    NIRS.Cf.H.F.r.o.mm.p = mat_Fid'; 
    %sources
    %NIRS.Cf.H.S.n
    NIRS.Cf.H.S.N = SD.n_Src;
    NIRS.Cf.H.S.r.o.mm.p = SD.SrcPos';
    %detectors
    %NIRS.Cf.H.D.n
    NIRS.Cf.H.D.N = SD.n_Det;
    NIRS.Cf.H.D.r.o.mm.p = SD.DetPos';
    %channels
    NIRS.Cf.H.C.n = [SD.ChannelLabels; SD.ChannelLabels];
    NIRS.Cf.H.C.N = size(SD.ml,1);
    NIRS.Cf.H.C.id = [1:size(SD.ml,1); SD.ml(:,1:2)'];
    NIRS.Cf.H.C.wl =SD.ml(:,4)';
      %for Files_idx -- Big loop over BOXY files for this subject  
    %DISTANCE GEOMETRIQUE
    for i=1:size(SD.ml,1)
        x1 =  NIRS.Cf.H.S.r.o.mm.p(1,SD.ml(i,1));
        x2 =  NIRS.Cf.H.D.r.o.mm.p(1,SD.ml(i,2));
        y1 =  NIRS.Cf.H.S.r.o.mm.p(2,SD.ml(i,1));
        y2 =  NIRS.Cf.H.D.r.o.mm.p(2,SD.ml(i,2));
        z1 =  NIRS.Cf.H.S.r.o.mm.p(3,SD.ml(i,1));
        z2 =  NIRS.Cf.H.D.r.o.mm.p(3,SD.ml(i,2));
        NIRS.Cf.H.C.gp(i,1) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);%À CALCULER DISTANCE GÉOMETRIQUE
        pos(i,1) = (x2-x1)/2+x1;
        pos(i,2) = (y2-y1)/2+y1;
        pos(i,3) = (z2-z1)/2+z1; 
        pos(i,4) = NIRS.Cf.H.C.gp(i,1);
    end
    
    %idbad = find(NIRS.Cf.H.C.gp<job.distmin|NIRS.Cf.H.C.gp>job.distmax);
    NIRS.Cf.H.C.ok = SD.ok;
    


end
    %write out NIRS in .mat file
    save(fullfile(pathout,'NIRS.mat'),'NIRS');   
   % save([oldNIRS.subj_path,job.config_path.output_path, '\', 'oldNIRS'],'oldNIRS');   
    outNIRSmat = [outNIRSmat; fullfile(pathout,'NIRS.mat')];
out.NIRSmat = outNIRSmat;