function newPrj = create_project_SNIRF(filedata)
    % This function creates a new project object from the SNIRF file
    % information contained in the filedata structure and returns it.
    
    % Template folder TODO
    path = ['.',filesep,'MontageTemplate',filesep];
    
    % Create a project object.
    newPrj = IO_Project_Data;

    % Create a digitization for the subject fiducials (SNIRF data).
    dig_subjectHelmet = Digitization;

        % Digitization of the subject probes and fiducials.
    subjectProbes = struct;

        % Detectors
    for detectorIndex = size(filedata.detectorPosition,1):-1:1                                     % Iterate through all detectors and add them to the probe array.
        subjectProbes(detectorIndex).Type = 400;
        subjectProbes(detectorIndex).Label = char(filedata.detectorLabels(detectorIndex));
        subjectProbes(detectorIndex).Coord = [filedata.detectorPosition(detectorIndex,:)./100, 1]; % Divide by 100, otherwise the helmet would be too big for the fids.
    end                                                                                            % Note however that this is probably an error with the unit detection
                                                                                                   % and show be addressed in the read/write code. Appending a 1 is
        % Sources                                                                                  % necessary to respect the probe format.
    for sourceIndex = size(filedata.sourcePosition,1):-1:1
        sourceDataTemp(sourceIndex).Type = 400;
        sourceDataTemp(sourceIndex).Label = char(filedata.sourceLabels(sourceIndex));
        sourceDataTemp(sourceIndex).Coord = [filedata.sourcePosition(sourceIndex,:)./100, 1];      % ^^^ Same as above. ^^^
    end
    subjectProbes = [subjectProbes, sourceDataTemp];

        % Fiducials
    if isfield(filedata,'landmarkLabels')
        % Find the index of the corresponding landmark. To be recognized,
        % the landmark has to have a standard name (case insensitive).
        nasIndex = find(upper(filedata.landmarkLabels) == "NAS" | upper(filedata.landmarkLabels) == "NASION");
        lpaIndex = find(upper(filedata.landmarkLabels) == "LPA" | upper(filedata.landmarkLabels) == "LEFT PRE-AURICULAR");
        rpaIndex = find(upper(filedata.landmarkLabels) == "RPA" | upper(filedata.landmarkLabels) == "RIGHT PRE-AURICULAR");
        % Verify only one index was found for each landmark.
        assert(length(nasIndex)==1 && length(lpaIndex)==1 && length(rpaIndex)==1, "Error. Missing or duplicate landmark.")
        % Add the fiducial coordinates to the SaveData.
        subjectFiducials = [filedata.landmarkPosition(nasIndex,:); filedata.landmarkPosition(lpaIndex,:); filedata.landmarkPosition(rpaIndex,:)];
    else
        error("No fiducial data found in the SNIRF file.")
    end

        % Add the probes and the fiducials to the subject helmet digitization.
    dig_subjectHelmet = set_vProbes(dig_subjectHelmet, subjectProbes);
    dig_subjectHelmet = set_matFiducials(dig_subjectHelmet, subjectFiducials);
    newPrj = set_Dig_SubjectFids(newPrj, dig_subjectHelmet);
    
    % Create a helmet object from the partial helmet and the subject
    % fiducials.
    helmet = Helmet(get_Dig_SubjectFids(newPrj));
    newPrj = set_Helmet(newPrj, helmet);
    
    % Create a MRI data object by reading the template MRI file.
    fileroot = 'template_fid';
    [oMRI, matVoxels] = load_MRI(path, fileroot);
    
    % Load segmentation data.
    
    fileroot = 'template_cortex';
    cortexSeg = load_segmentation(path, fileroot);
    
    fileroot = 'template_skin';
    skinSeg = load_segmentation(path, fileroot);
    
    % Treat segmentation data.
    
    bSkinSeg = false;
    bCortexSeg = false;
    
    % Add cortex segmentation data.
    
    [nVoxX, nVoxZ, nVoxY ] = size(matVoxels);
    clear matVoxels;
    [nVoxCortexSegX, nVoxCortexSegZ, nVoxCortexSegY ] = size(cortexSeg.matSegmentation);

        %Validation de la taille de la segmentation du cortex
    if( nVoxCortexSegX == nVoxX && nVoxCortexSegZ == nVoxZ && nVoxCortexSegY <= nVoxY )

        %Uniformisation des tailles
        if( nVoxCortexSegY < nVoxY )
                cortexSeg.matSegmentation(:,:,nVoxCortexSegY+1:nVoxY) = zeros(nVoxX,nVoxZ,nVoxY-nVoxCortexSegY, 'uint8');                    

            %Formattage de la matrice de segmentation:
            %->(X,Y,Z): C'est plus intuitif
            tmpSeg = zeros( nVoxX, nVoxY, nVoxZ, 'uint8' );
            disp('Reshaping Segmentation');
            for iY=1:nVoxY
                tmpSeg(:,iY,:) = reshape( cortexSeg.matSegmentation(:,:,iY), nVoxX, 1, nVoxZ );
            end
        disp('Done');

        oMRI = set_CortexSegmentation( oMRI, tmpSeg );
        bCortexSeg = true;

        else
            h = 10;
            warndlg( 'La segmentation n''à pu etre chargée' );
            waitfor(h);
            warndlg( 'Les dimensions XZ du fichier de segmentation de cortex ne correspondent pas aux dimensions du fichier de données' );
            waitfor(h);
        end
    end
    
    % Add skin segmentation data.
    
    [nVoxSkinSegX, nVoxSkinSegZ, nVoxSkinSegY ] = size(skinSeg.matSegmentation);
    
        %Validation de la taille de la segmentation de la peau
    if( nVoxSkinSegX == nVoxX && nVoxSkinSegZ == nVoxZ && nVoxSkinSegY <= nVoxY )

        %Uniformisation des tailles
        if( nVoxSkinSegY < nVoxY )
            skinSeg.matSegmentation(:,:,nVoxSkinSegY+1:nVoxY) = zeros(nVoxX,nVoxZ,nVoxY-nVoxSkinSegY, 'uint8');                    
        end

        %Formattage de la matrice de segmentation:
        %->(X,Y,Z): C'est plus intuitif
        tmpSeg = zeros( nVoxX, nVoxY, nVoxZ, 'uint8' );
        disp('Reshaping Segmentation');
        for iY=1:nVoxY
            tmpSeg(:,iY,:) = reshape( skinSeg.matSegmentation(:,:,iY), nVoxX, 1, nVoxZ );
        end
        disp('Done');

        oMRI = set_SkinSegmentation( oMRI, tmpSeg );
        bSkinSeg = true;

    else
        h = warndlg( 'La segmentation n''à pu etre chargée' );
        waitfor(h);
        h = warndlg( 'Les dimensions XZ du fichier de segmentation de peau ne correspondent pas aux dimensions du fichier de données' );
        waitfor(h);
    end
    
    Helm = get_Helmet(newPrj);
    
    % Esthablish a coregistration matrix.
    oMRI = CoRegisterToHelmet(oMRI, Helm);
    
    % Fit holes on the segmentation.
    Helm = Fit_Helmet_OnMRISegmentations(Helm, oMRI, bSkinSeg, bCortexSeg, false);
    newPrj = set_Helmet(newPrj, Helm);
    
    % Load meshes.
        % Load cortex hight vertex and low vertex.
    fileroot = 'template_cortex';
    FullFileName_CortexHiMesh = [fullfile(path, fileroot), '.SRX'];
    
    [vertex_matrix, faces_matrix, ~, ~] = Read_iMagic_SRX( FullFileName_CortexHiMesh );
    CortexMeshHiVertex = vertex_matrix;
    CortexMeshHiIndex = faces_matrix;
    CortexMeshHiVcolor = zeros(size(CortexMeshHiVertex,1),1);
    
    CortexMeshLowRes = CortexMeshHiVertex;
    CortexMeshLoIndex = CortexMeshHiIndex;
    CortexMeshLoVcolor = CortexMeshHiVcolor;
   
        % Load skin vertex.
    fileroot = 'template_skin';
    FullFileName_SkinMesh = [fullfile(path, fileroot), '.SRX'];
    [vertex_matrix, faces_matrix, ~, ~] = Read_iMagic_SRX( FullFileName_SkinMesh );
    SkinMeshVertex = vertex_matrix;
    SkinMeshIndex = faces_matrix;
    SkinVcolor = zeros(size(SkinMeshVertex,1),1);
    
        % Add meshes to the MRI data.
            % HiRes cortex mesh.
    oMRI = set_CortexMeshHiRes(oMRI, CortexMeshHiVertex,CortexMeshHiIndex);
    oMRI = set_CortexHiResVcolor(oMRI, CortexMeshHiVcolor);
            % LoRes cortex mesh.
    oMRI = set_CortexMeshLowRes(oMRI, CortexMeshLowRes, CortexMeshLoIndex);
    oMRI = set_CortexLowResVcolor(oMRI, CortexMeshLoVcolor);
            % Skin mesh.
    oMRI = set_SkinMesh(oMRI, SkinMeshVertex, SkinMeshIndex);
    oMRI = set_SkinVcolor(oMRI, SkinVcolor);
    
    % Set atlas of MRI.
    fileroot = 'template_Atlas';
    FullFileName_Atlas = [fullfile(path, fileroot), '.img'];
    oMRI = coloratlas(oMRI, FullFileName_Atlas);
    
    % Set MRI data.
    newPrj = set_MRI_Data(newPrj, oMRI);
end

function [mri_data, matVoxels] = load_MRI(path, fileroot)
    % This function loads the MRI data in the .hdr and .vox file from the
    % template folder and creates a MRI data object that it then returns.

    % Get the complete filepath for the .hdr and .vox file.
    FullFileName_Hdr = [fullfile(path, fileroot), '.hdr'];
    FullFileName_Vox = [fullfile(path, fileroot), '.VOX'];
    
    %----------------------------------------------------------------------
    %--------------- SECTION HDR (ASCII) ----------------------------------
    %----------------------------------------------------------------------
        
    fid_Hdr = fopen(FullFileName_Hdr);
    if( fid_Hdr < 2 )
        warndlg( 'Error loading Volume Data File (ERR113)' );
    end

    %XRES
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR114)');
        return;
    end
    XRes=strread(a_line, '%d');

    %YRES
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR115)');
        return;
    end
    YRes=strread(a_line, '%d');

    %ZRES
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR116)');
        return;
    end
    ZRes=strread(a_line, '%d');

    %Bytes per pixels (jumped)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR117)');
        return;
    end
    BytesPerPix=strread(a_line, '%d');

    %XZDim (vox dim inside a slice), in mm
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR118)');
        return;
    end
    XZVoxDim = strread(a_line);
    if( numel(XZVoxDim) == 2 )
        XDim = XZVoxDim(1);
        ZDim = XZVoxDim(2);
    else
        XDim = XZVoxDim; %mm
        ZDim = XZVoxDim; %mm
    end

    %YDim (slice thickness) in mm
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR119)');
        return;
    end
    YDim=strread(a_line);

    %Plane
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR120)');
        return;
    end
    Plane=strread(a_line, '%d');

    %Jump this line - Skin treshold
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR121)');
        return;
    end

    %Jump this line - Cortex treshold
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR122)');
        return;
    end

    %LPA
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR123)');
        return;
    end
    voxFids(2,:)=strread(a_line);

    %RPA
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR124)');
        return;
    end
    voxFids(3,:)=strread(a_line);

    %NAS
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR125)');
        return;
    end
    voxFids(1,:)=strread(a_line);

    %Jump this line (VERTEX)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR126)');
        return;
    end

    %Jump this line (ANTCOMM? - TAILARACH SYSTEM ORIGIN)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR127)');
        return;
    end

    %Jump this line - ATLAS filename associated
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR128)');
        return;
    end

    %Jump this line - MIP filename associated
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR129)');
        return;
    end

    %Jump this line (INI)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR130)');
        return;
    end

    %NOS
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR131)');
        return;
    end
    voxFids(4,:)=strread(a_line);

    if isempty(find(voxFids,1))
        warndlg('No fiducial point marked (ERR132)');
        return;
    end

    %----------------------------------------------------------------------
    %--------------- SECTION VOX (BINARY) ---------------------------------
    %----------------------------------------------------------------------
    fid_Vox=fopen(FullFileName_Vox);

    if( BytesPerPix == 1 )
        sDataType = 'uint8';
        sDataConversion = @uint8;
    elseif( BytesPerPix == 2 )
        sDataType = 'uint16';
        sDataConversion = @uint16;
    else
        warndlg('Error loading Volume Data File: unsupported bbp (ERR132)');
        return;
    end

    %Initialisation
    countTotal = 0;
    MaxPacketSize = 2500000;
    VoxelsToRead = XRes*YRes*ZRes;
    Data1D = zeros(VoxelsToRead, 1, sDataType);

    hWaitbar = waitbar(0,'Loading data...');
    while( countTotal < VoxelsToRead )

        %Lit en format 'double': convertir en uint8/16
        [Packet, count] =  fread(fid_Vox, MaxPacketSize, sDataType);

            Data1D(countTotal+1:countTotal+count) = sDataConversion(Packet);

        countTotal = countTotal+count;
        waitbar(countTotal/VoxelsToRead);
    end
    close(hWaitbar);

    matVoxels = sDataConversion(reshape( Data1D, XRes, ZRes, YRes ));

    clear Data1D % Clear the temporary matrix necessary for the reshape operation.

    fclose(fid_Hdr);
    fclose(fid_Vox);

    % Treat the MRI data and create MRI object.

        %Les dimensions ne sont pas tres intuitives (X,Z,Y)...
    [nVoxX, nVoxZ, nVoxY ] = size(matVoxels);

        %Initialisation
    sMRI.img = zeros( nVoxX, nVoxY, nVoxZ, sDataType );

        %Formattage de la matrice de voxels ->(X,Y,Z): C'est plus intuitif
    disp('Reshaping Voxels');
    for iY=1:nVoxY
        sMRI.img(:,iY,:) = reshape( matVoxels(:,:,iY), nVoxX, 1, nVoxZ );
    end
    disp('Done');

        %Sauvegarde des informations
    sMRI.hdr.dime.bitpix = BytesPerPix;
    sMRI.hdr.dime.pixdim = [ 0, XDim, YDim, ZDim ];
    sMRI.hdr.dime.glmax = max(max(max(matVoxels)));

        %Constructeur objet MRI_Data (prend la struct en parametre)
    mri_data = MRI_Data( sMRI ); 
    
        %Sauvegarde des fiducies: handles.VOXFids = [NAS; LPA; RPA]
    mri_data = set_matFiducials( mri_data, voxFids(1:3,:) );

end

function segmentationStruct = load_segmentation(path, filename)
    % This function loads a segmentation into a struct for further
    % processing.
    
    FullFileName_Hdr = [fullfile(path, filename), '.hdr'];
    FullFileName_Seg = [fullfile(path, filename), '.SEG'];
    
    segmentationStruct = struct;
    
    fid_Hdr = fopen(FullFileName_Hdr);
    if( fid_Hdr < 2 )
        warndlg( 'Error loading Segmentation File (ERR103)' );
    end
        
    %Taille X de la segmentation (Horizontal)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR104)' );
        matSegmentation = [];
        return;
    end
    XSegRes=strread(a_line, '%d' );
    segmentationStruct.XSegRes = XSegRes;
    
    
    %Taille Y de la segmentation (Slice)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR105)' );
        matSegmentation = [];
        return;
    end
    YSegRes=strread(a_line, '%d' );
    segmentationStruct.YSegRes = YSegRes;
    
                
    %Taille Z de la segmentation (Vertical)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR106)' );
        matSegmentation = [];
        return;
    end
    ZSegRes=strread(a_line, '%d' );
    segmentationStruct.ZSegRes = ZSegRes;
    
    
    %Index des images
    vIndex = [];
    a_line=fgets(fid_Hdr);
    while( a_line ~= -1 )
        %Lecture de l'index et correction de 0-indexed a 1-indexed
        vIndex(numel(vIndex)+1) = strread(a_line, '%d' )+1;
        a_line=fgets(fid_Hdr);
    end
    MaxSliceNo = max( vIndex );
    segmentationStruct.vIndexSeg = vIndex;

    %----------------------------------------------------------------------
    %--------------- SECTION SEG (BINARY) ---------------------------------
    %----------------------------------------------------------------------
    fid_Seg=fopen(FullFileName_Seg);
    
    %Initialisation
    countTotal = 0;
    MaxPacketSize = 2500000;
    VoxelsToRead = XSegRes*YSegRes*ZSegRes;
    Data1D = zeros(VoxelsToRead, 1, 'uint8');
        
    hWaitbar = waitbar(0,'Loading data...');
    while( countTotal < VoxelsToRead )
        
        %Lit en format 'double': convertir en uint8/16
        [Packet, count] =  fread(fid_Seg, MaxPacketSize, 'uint8');
        
        Data1D(countTotal+1:countTotal+count) = uint8(Packet);
      
        countTotal = countTotal+count;
        
        waitbar(countTotal/VoxelsToRead);
    end
    close(hWaitbar);
    
    mat_tmp = uint8(reshape( Data1D, XSegRes, ZSegRes, YSegRes ));
    
    segmentationStruct.matSegmentation = zeros( XSegRes, ZSegRes, MaxSliceNo, 'uint8' );
    
    segmentationStruct.matSegmentation( :, :, vIndex ) = mat_tmp;

    fclose(fid_Hdr);
    fclose(fid_Seg);
end