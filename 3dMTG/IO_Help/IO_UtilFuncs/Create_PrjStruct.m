function SaveStruct = Create_PrjStruct( oPrjData )
    
    %initialisation
    SaveStruct = [];
    
    %----------------------------------------------------------------------
    %             ---------- SECTION DIGITIZATIONS ----------
    %----------------------------------------------------------------------
    %Nouvellement ajouté SB
    %Transfer de la digitalisation de casque entier (Celle dont l'origine est a la pointe du "Brainsight Calibration Block", en format LOCATOR) 
    oDig_CompleteHelmet = get_Dig_CompleteHelmet( oPrjData );
    vProbesCH = get_vProbes( oDig_CompleteHelmet );
    for( p=1:numel(vProbesCH) )
        SaveStruct.m_DigCompHelm.v_Probes(p).Type   = vProbesCH(p).Type;
        SaveStruct.m_DigCompHelm.v_Probes(p).Label  = vProbesCH(p).Label;
        SaveStruct.m_DigCompHelm.v_Probes(p).Coord  = vProbesCH(p).Coord;
    end
    SaveStruct.m_DigCompHelm.matFiducials = get_matFiducials( oDig_CompleteHelmet );

    %Nouvellement ajouté SB
    %Transfer de la digitalisation de casque partiel (Fids+X1X2X3X4+possiblement quelques points supplementaires) 
    oDig_SubjectFids = get_Dig_SubjectFids(oPrjData);
    vProbesSF = get_vProbes( oDig_SubjectFids );
    for( p=1:numel(vProbesSF) )
        SaveStruct.m_DigSubjFids.v_Probes(p).Type  = vProbesSF(p).Type;
        SaveStruct.m_DigSubjFids.v_Probes(p).Label = vProbesSF(p).Label;
        SaveStruct.m_DigSubjFids.v_Probes(p).Coord = vProbesSF(p).Coord;
    end
    SaveStruct.m_DigSubjFids.matFiducials = get_matFiducials( oDig_SubjectFids );
    
    %----------------------------------------------------------------------
    %                 ---------- SECTION HELMET ----------
    %----------------------------------------------------------------------
    oHelmet = get_Helmet(oPrjData);
    vHoles = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );
    
    for( p=1:numel(vHoles) )
        SaveStruct.m_Helmet.v_Holes(p).Type      = vHoles(p).Type;
        SaveStruct.m_Helmet.v_Holes(p).Label     = vHoles(p).Label;
        SaveStruct.m_Helmet.v_Holes(p).Coord.x   = vHoles(p).Coord.x;
        SaveStruct.m_Helmet.v_Holes(p).Coord.y   = vHoles(p).Coord.y;
        SaveStruct.m_Helmet.v_Holes(p).Coord.z   = vHoles(p).Coord.z;
        SaveStruct.m_Helmet.v_Holes(p).CanBeDet  = vHoles(p).CanBeDet;
        SaveStruct.m_Helmet.v_Holes(p).CanBeSrc  = vHoles(p).CanBeSrc;
        SaveStruct.m_Helmet.v_Holes(p).RegionsID = vHoles(p).RegionsID;
        SaveStruct.m_Helmet.v_Holes(p).Normal.x  = vHoles(p).Normal.x;
        SaveStruct.m_Helmet.v_Holes(p).Normal.y  = vHoles(p).Normal.y;
        SaveStruct.m_Helmet.v_Holes(p).Normal.z  = vHoles(p).Normal.z;
        
        %Nouvellement ajouté SB
        SaveStruct.m_Helmet.v_Holes(p).Transformation       = vHoles(p).Transformation;    %Matrice de transformation (pre-calcul, pour l'orientation de chaque trou. Voir Draw_Holes.m)
        SaveStruct.m_Helmet.v_Holes(p).SkinDepth            = vHoles(p).SkinDepth;         
        SaveStruct.m_Helmet.v_Holes(p).CortexDepth          = vHoles(p).CortexDepth; 
        SaveStruct.m_Helmet.v_Holes(p).IsFromCompleteHelmet = vHoles(p).IsFromCompleteHelmet;
    end

    SaveStruct.m_Helmet.Mtg_Data.mat_RegionNorms  = sMtg.mat_RegionNorms;
    SaveStruct.m_Helmet.Mtg_Data.mat_RegionBorder = sMtg.mat_RegionBorder;
    SaveStruct.m_Helmet.Mtg_Data.matFiducials     = sMtg.matFiducials;
    SaveStruct.m_Helmet.Mtg_Data.v_HolesMtg       = sMtg.v_HolesMtg;
    SaveStruct.m_Helmet.Mtg_Data.v_HolesEle       = sMtg.v_HolesEle;
    SaveStruct.m_Helmet.Mtg_Data.v_pSrc           = sMtg.v_pSrc;
    SaveStruct.m_Helmet.Mtg_Data.v_pDet           = sMtg.v_pDet;
    SaveStruct.m_Helmet.Mtg_Data.v_pEle           = sMtg.v_pEle;
    SaveStruct.m_Helmet.Mtg_Data.v_HolesMtg_Selected = sMtg.v_HolesMtg_Selected; %JT 1 = selectionne 0 = non selectionne
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.Min_Dist_Mux      = sMtg.Gen_Params.Min_Dist_Mux;
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.Max_Dist_Mux      = sMtg.Gen_Params.Max_Dist_Mux;
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.DistContamination = sMtg.Gen_Params.DistContamination;
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.Nb_Sources_Disp   = sMtg.Gen_Params.Nb_Sources_Disp;
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.Nb_Detect_Disp    = sMtg.Gen_Params.Nb_Detect_Disp;
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.Nb_Longueur_Onde  = sMtg.Gen_Params.Nb_Longueur_Onde;
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.NbBanks           = sMtg.Gen_Params.NbBanks;
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.NbSrcPerBank      = sMtg.Gen_Params.NbSrcPerBank;
    if ~isfield(sMtg.Gen_Params,'AcqSystem')
        sMtg.Gen_Params.AcqSystem = 'ISS';
    end
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.AcqSystem         = sMtg.Gen_Params.AcqSystem;
    
    %Nouvellement ajouté SB
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.UniformSkinDepth     = sMtg.Gen_Params.UniformSkinDepth;
    SaveStruct.m_Helmet.Mtg_Data.Gen_Params.bUseUniformSkinDepth = sMtg.Gen_Params.bUseUniformSkinDepth;
    
    %----------------------------------------------------------------------
    %             ---------- SECTION MRI_DATA ----------
    %----------------------------------------------------------------------
    oMRI_Data = get_MRI_Data( oPrjData );
    vMarkers = get_vMarkers( oMRI_Data );
    [VertexBufferCortexHiRes,  IndexBufferCortexHiRes]  = get_CortexMeshHiRes( oMRI_Data );
    [VertexBufferCortexLowRes, IndexBufferCortexLowRes] = get_CortexMeshLowRes( oMRI_Data );
    [VertexBufferSkin,         IndexBufferSkin]         = get_SkinMesh( oMRI_Data );
    
    Vcolormri = get_CortexHiResVcolor( oMRI_Data );

    %Nouvellement ajouté SB
    for( p=1:numel(vMarkers) )
        SaveStruct.m_MRI_Data.vMarkers(p).Coord = vMarkers(p).Coord;
        SaveStruct.m_MRI_Data.vMarkers(p).Label = vMarkers(p).Label;
    end
    SaveStruct.m_MRI_Data.matVoxels             = get_matVoxels( oMRI_Data );
    SaveStruct.m_MRI_Data.matCortexSegmentation = get_CortexSegmentation( oMRI_Data );
    SaveStruct.m_MRI_Data.matSkinSegmentation   = get_SkinSegmentation( oMRI_Data );
    SaveStruct.m_MRI_Data.matFiducials          = get_matFiducials( oMRI_Data );
    SaveStruct.m_MRI_Data.matCoRegistration     = get_matTransform( oMRI_Data );
    SaveStruct.m_MRI_Data.VoxDim                = get_PhysicalVoxDim( oMRI_Data );
    SaveStruct.m_MRI_Data.VoxDepth              = get_VoxDepth( oMRI_Data );
    SaveStruct.m_MRI_Data.Vox_IMax              = get_IMaxVox( oMRI_Data );
    SaveStruct.m_MRI_Data.Vox_IMin              = get_IMinVox( oMRI_Data );
    SaveStruct.m_MRI_Data.Surfaces.CortexHiRes.VertexBuffer  = VertexBufferCortexHiRes;
    SaveStruct.m_MRI_Data.Surfaces.CortexHiRes.IndexBuffer   = IndexBufferCortexHiRes;
    SaveStruct.m_MRI_Data.Surfaces.Skin.Vcolor               = zeros(size(VertexBufferSkin,1),1);
    SaveStruct.m_MRI_Data.Surfaces.CortexLowRes.VertexBuffer = VertexBufferCortexLowRes;
    SaveStruct.m_MRI_Data.Surfaces.CortexLowRes.IndexBuffer  = IndexBufferCortexLowRes;
    SaveStruct.m_MRI_Data.Surfaces.Skin.VertexBuffer         = VertexBufferSkin;
    SaveStruct.m_MRI_Data.Surfaces.Skin.IndexBuffer          = IndexBufferSkin;
    SaveStruct.m_MRI_Data.Surfaces.CortexHiResVcolor         = Vcolormri;
        