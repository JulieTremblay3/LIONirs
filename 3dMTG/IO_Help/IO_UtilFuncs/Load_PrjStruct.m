function oPrjData = Load_PrjStruct( LoadedStruct, bComputeNeighborhood )
    
    %Initialisation d'un nouveau projet par le constructeur
    oPrjData  = IO_Project_Data();
    
    %La structure prend le nom de la variable "SaveStruct" quand elle est
    %sauvegardee. Remettre le champ SaveStruct dans LoadedStruct.
    if( isfield(LoadedStruct, 'SaveStruct') )
        LoadedStruct = LoadedStruct.SaveStruct;
    end
    
    %Gestion des versions de sauvegarde precedentes (avant les strcutures
    %de sauvegarde)
    if( ~isfield(LoadedStruct, 'm_Helmet') )
        
        if( isfield(LoadedStruct, 'PrjData') )
            if( isobject(LoadedStruct.PrjData) )
                LoadedStruct.m_Helmet = get_Helmet( LoadedStruct.PrjData );
            else
                LoadedStruct.m_Helmet = LoadedStruct.PrjData.m_Helmet;
            end
        else
            errordlg( 'Invalid Project Data File', 'Error' );
            return;
        end
    end
    
    %----------------------------------------------------------------------
    %             ---------- SECTION DIGITIZATIONS ----------
    %----------------------------------------------------------------------
    %Transfer de la digitalisation de casque entier (Celle dont l'origine est a la pointe du "Brainsight Calibration Block", en format LOCATOR) 
    %New Digitization - Complete Helmet (
    oDig_CompHelm = Digitization();
    vProbesCH = get_vProbes( oDig_CompHelm );
    matFiducials = get_matFiducials( oDig_CompHelm );
    EmptyProbe = vProbesCH(1);
    
    %Positions de trous et informations de trous
    if( isfield(LoadedStruct, 'm_DigCompHelm') )
        if( isfield(LoadedStruct.m_DigCompHelm, 'v_Probes') )
            for( p=1:numel(LoadedStruct.m_DigCompHelm.v_Probes) )
                %Initialisation
                vProbesCH(p) = EmptyProbe;
                
                if( isfield( LoadedStruct.m_DigCompHelm.v_Probes(p), 'Type' ) )
                    vProbesCH(p).Type   = LoadedStruct.m_DigCompHelm.v_Probes(p).Type;
                end

                if( isfield( LoadedStruct.m_DigCompHelm.v_Probes(p), 'Label' ) )
                    vProbesCH(p).Label  = LoadedStruct.m_DigCompHelm.v_Probes(p).Label;
                end
                
                if( isfield( LoadedStruct.m_DigCompHelm.v_Probes(p), 'Coord' ) )
                    vProbesCH(p).Coord  = LoadedStruct.m_DigCompHelm.v_Probes(p).Coord;
                end
            end
        end
        
        if( isfield(LoadedStruct.m_DigCompHelm, 'matFiducials') )
            matFiducials = LoadedStruct.m_DigCompHelm.matFiducials;
        end
    end
    oDig_CompHelm = set_vProbes( oDig_CompHelm, vProbesCH );
    oDig_CompHelm = set_matFiducials( oDig_CompHelm, matFiducials );
    oPrjData = set_Dig_CompleteHelmet( oPrjData, oDig_CompHelm );
    
    
    
    %Transfer de la digitalisation de casque partiel (Fids+X1X2X3X4+possiblement quelques points supplementaires) 
    %New Digitization
    oDig_SubjFids = Digitization();
    vProbesSF = get_vProbes( oDig_SubjFids );
    matFiducials = get_matFiducials( oDig_SubjFids );
    EmptyProbe = vProbesSF(1);
    
    %Positions de trous et informations de trous
    if( isfield(LoadedStruct, 'm_DigSubjFids') )
        if( isfield(LoadedStruct.m_DigSubjFids, 'v_Probes') )
            for( p=1:numel(LoadedStruct.m_DigSubjFids.v_Probes) )
                %Initialisation
                vProbesSF(p) = EmptyProbe;
                
                if( isfield( LoadedStruct.m_DigSubjFids.v_Probes(p), 'Type' ) )
                    vProbesSF(p).Type   = LoadedStruct.m_DigSubjFids.v_Probes(p).Type;
                end

                if( isfield( LoadedStruct.m_DigSubjFids.v_Probes(p), 'Label' ) )
                    vProbesSF(p).Label  = LoadedStruct.m_DigSubjFids.v_Probes(p).Label;
                end
                
                if( isfield( LoadedStruct.m_DigSubjFids.v_Probes(p), 'Coord' ) )
                    vProbesSF(p).Coord  = LoadedStruct.m_DigSubjFids.v_Probes(p).Coord;
                end
            end
        end
        
        if( isfield(LoadedStruct.m_DigSubjFids, 'matFiducials') )
            matFiducials = LoadedStruct.m_DigSubjFids.matFiducials;
        end
    end
    oDig_SubjFids = set_vProbes( oDig_SubjFids, vProbesSF );
    oDig_SubjFids = set_matFiducials( oDig_SubjFids, matFiducials );
    oPrjData = set_Dig_SubjectFids( oPrjData, oDig_SubjFids );
    
    %----------------------------------------------------------------------
    %                 ---------- SECTION HELMET ----------
    %----------------------------------------------------------------------
    %New Helmet
    oHelmet = Helmet(); %Constructeur de Helmet
    vHoles = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );
    EmptyHole = vHoles(1);
    
    %Positions de trous et informations de trous
    if( isfield(LoadedStruct.m_Helmet, 'v_Holes') )
        for( p=1:numel(LoadedStruct.m_Helmet.v_Holes) )
            
            %Initialisation
            vHoles(p) = EmptyHole;
            
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'Type' ) )
                vHoles(p).Type      = LoadedStruct.m_Helmet.v_Holes(p).Type;
            end
            
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'Label' ) )
                vHoles(p).Label     = LoadedStruct.m_Helmet.v_Holes(p).Label;
            end
            
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'Coord' ) )
                if( isfield( LoadedStruct.m_Helmet.v_Holes(p).Coord, 'x' ) )
                    vHoles(p).Coord.x   = LoadedStruct.m_Helmet.v_Holes(p).Coord.x;
                end

                if( isfield( LoadedStruct.m_Helmet.v_Holes(p).Coord, 'y' ) )
                    vHoles(p).Coord.y   = LoadedStruct.m_Helmet.v_Holes(p).Coord.y;
                end

                if( isfield( LoadedStruct.m_Helmet.v_Holes(p).Coord, 'z' ) )
                    vHoles(p).Coord.z   = LoadedStruct.m_Helmet.v_Holes(p).Coord.z;
                end
            end
            
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'Normal' ) )
                if( isfield( LoadedStruct.m_Helmet.v_Holes(p).Normal, 'x' ) )
                    vHoles(p).Normal.x  = LoadedStruct.m_Helmet.v_Holes(p).Normal.x;
                end

                if( isfield( LoadedStruct.m_Helmet.v_Holes(p).Normal, 'y' ) )
                    vHoles(p).Normal.y  = LoadedStruct.m_Helmet.v_Holes(p).Normal.y;
                end

                if( isfield( LoadedStruct.m_Helmet.v_Holes(p).Normal, 'z' ) )
                    vHoles(p).Normal.z  = LoadedStruct.m_Helmet.v_Holes(p).Normal.z;
                end
            end     
            
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'CanBeDet' ) )
                vHoles(p).CanBeDet  = LoadedStruct.m_Helmet.v_Holes(p).CanBeDet;
            end
            
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'CanBeSrc' ) )
                vHoles(p).CanBeSrc  = LoadedStruct.m_Helmet.v_Holes(p).CanBeSrc;
            end

            %Nouvellement ajoute SB
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'Transformation' ) )
                vHoles(p).Transformation = LoadedStruct.m_Helmet.v_Holes(p).Transformation;
            end 
            
            %Nouvellement ajoute SB
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'SkinDepth' ) )
                vHoles(p).SkinDepth = LoadedStruct.m_Helmet.v_Holes(p).SkinDepth;
            end 
            
            %Nouvellement ajoute SB
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'CortexDepth' ) )
                vHoles(p).CortexDepth = LoadedStruct.m_Helmet.v_Holes(p).CortexDepth;
            end 
            
            %Nouvellement ajoute SB
            if( isfield( LoadedStruct.m_Helmet.v_Holes(p), 'IsFromCompleteHelmet' ) )
                vHoles(p).IsFromCompleteHelmet = LoadedStruct.m_Helmet.v_Holes(p).IsFromCompleteHelmet;
            end 
        end
    end
    
    %Chargement du montage
    if( isfield(LoadedStruct.m_Helmet, 'Mtg_Data') )
        
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'mat_RegionNorms') )
            sMtg.mat_RegionNorms  = LoadedStruct.m_Helmet.Mtg_Data.mat_RegionNorms;
        end
        
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'mat_RegionBorder') )
            sMtg.mat_RegionBorder = LoadedStruct.m_Helmet.Mtg_Data.mat_RegionBorder;
        end
        
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'matFiducials') )
            sMtg.matFiducials     = LoadedStruct.m_Helmet.Mtg_Data.matFiducials;
        end
        
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'v_HolesMtg') )
            sMtg.v_HolesMtg       = LoadedStruct.m_Helmet.Mtg_Data.v_HolesMtg;
        end
         
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'v_HolesEle') )
            sMtg.v_HolesEle       = LoadedStruct.m_Helmet.Mtg_Data.v_HolesEle;
        end
        
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'v_pSrc') )
            sMtg.v_pSrc           = LoadedStruct.m_Helmet.Mtg_Data.v_pSrc;
        end
        
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'v_pDet') )
            sMtg.v_pDet           = LoadedStruct.m_Helmet.Mtg_Data.v_pDet;
        end
        
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'v_pEle') )
            sMtg.v_pEle           = LoadedStruct.m_Helmet.Mtg_Data.v_pEle;
        end
        
        %Parametres de montage
        if( isfield(LoadedStruct.m_Helmet.Mtg_Data, 'Gen_Params') )
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'Min_Dist_Mux') )
                sMtg.Gen_Params.Min_Dist_Mux      = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.Min_Dist_Mux;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'Max_Dist_Mux') )
                sMtg.Gen_Params.Max_Dist_Mux      = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.Max_Dist_Mux;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'DistContamination') )
                sMtg.Gen_Params.DistContamination = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.DistContamination;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'Nb_Sources_Disp') )
                sMtg.Gen_Params.Nb_Sources_Disp   = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.Nb_Sources_Disp;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'Nb_Detect_Disp') )
                sMtg.Gen_Params.Nb_Detect_Disp    = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.Nb_Detect_Disp;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'Nb_Longueur_Onde') )
                sMtg.Gen_Params.Nb_Longueur_Onde  = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.Nb_Longueur_Onde;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'NbBanks') )
                sMtg.Gen_Params.NbBanks           = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.NbBanks;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'NbSrcPerBank') )
                sMtg.Gen_Params.NbSrcPerBank      = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.NbSrcPerBank;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'AcqSystem') )
                sMtg.Gen_Params.AcqSystem      = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.AcqSystem;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'UniformSkinDepth') )
                sMtg.Gen_Params.UniformSkinDepth  = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.UniformSkinDepth;
            end
            
            if( isfield(LoadedStruct.m_Helmet.Mtg_Data.Gen_Params, 'bUseUniformSkinDepth') )
                sMtg.Gen_Params.bUseUniformSkinDepth = LoadedStruct.m_Helmet.Mtg_Data.Gen_Params.bUseUniformSkinDepth;
            end
        end
    end
    
    %Chargement du nouveau casque 
    oHelmet = set_vHoles( oHelmet, vHoles );
    oHelmet = set_Mtg( oHelmet, sMtg );
    
    if( exist( 'bComputeNeighborhood' ) )
        if( bComputeNeighborhood )
            oHelmet = ComputeTransform( oHelmet );
            oHelmet = ComputeNeighborhood( oHelmet );
        end
    end
    oPrjData = set_Helmet( oPrjData, oHelmet );
    
    %----------------------------------------------------------------------
    %                 ---------- SECTION MRI_DATA ----------
    %----------------------------------------------------------------------
    %New Helmet
    oMRI = MRI_Data(); %Constructeur de MRI_DATA
    vMarkers = get_vMarkers( oMRI );
    EmptyMarker = vMarkers(1);
    
    %Positions de trous et informations de trous
    if( isfield(LoadedStruct, 'm_MRI_Data') )
        if( isfield(LoadedStruct.m_MRI_Data, 'vMarkers') )
            for( p=1:numel(LoadedStruct.m_MRI_Data.vMarkers) )

                %Initialisation
                vMarkers(p) = EmptyMarker;

                if( isfield( LoadedStruct.m_MRI_Data.vMarkers(p), 'Label' ) )
                    vMarkers(p).Label = LoadedStruct.m_MRI_Data.vMarkers(p).Label;
                end
                
                if( isfield( LoadedStruct.m_MRI_Data.vMarkers(p), 'Coord' ) )
                    vMarkers(p).Coord = LoadedStruct.m_MRI_Data.vMarkers(p).Coord;
                end
            end
            oMRI = set_vMarkers(oMRI, vMarkers);
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'matVoxels') )
            oMRI = set_matVoxels( oMRI, LoadedStruct.m_MRI_Data.matVoxels );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'matCortexSegmentation') )
            oMRI = set_CortexSegmentation( oMRI, LoadedStruct.m_MRI_Data.matCortexSegmentation );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'matSkinSegmentation') )
            oMRI = set_SkinSegmentation( oMRI, LoadedStruct.m_MRI_Data.matSkinSegmentation );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'matFiducials') )
            oMRI = set_matFiducials( oMRI, LoadedStruct.m_MRI_Data.matFiducials );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'matCoRegistration') )
            oMRI = set_matTransform( oMRI, LoadedStruct.m_MRI_Data.matCoRegistration );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'VoxDim') )
            oMRI = set_PhysicalVoxDim( oMRI, LoadedStruct.m_MRI_Data.VoxDim );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'VoxDepth') )
            oMRI = set_VoxDepth( oMRI, LoadedStruct.m_MRI_Data.VoxDepth );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'Vox_IMax') )
            oMRI = set_IMaxVox( oMRI, LoadedStruct.m_MRI_Data.Vox_IMax );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'Vox_IMin') )
            oMRI = set_IMinVox( oMRI, LoadedStruct.m_MRI_Data.Vox_IMin );
        end
        
        if( isfield(LoadedStruct.m_MRI_Data, 'Vox_IMin') )
            oMRI = set_IMinVox( oMRI, LoadedStruct.m_MRI_Data.Vox_IMin );
        end
         
        if( isfield(LoadedStruct.m_MRI_Data, 'Surfaces') )
           
            if( isfield(LoadedStruct.m_MRI_Data.Surfaces, 'CortexHiRes') )
                
                if( isfield(LoadedStruct.m_MRI_Data.Surfaces.CortexHiRes, 'VertexBuffer') && ...
                    isfield(LoadedStruct.m_MRI_Data.Surfaces.CortexHiRes, 'IndexBuffer') )
                    oMRI = set_CortexMeshHiRes( oMRI, LoadedStruct.m_MRI_Data.Surfaces.CortexHiRes.VertexBuffer, ...
                                                      LoadedStruct.m_MRI_Data.Surfaces.CortexHiRes.IndexBuffer );
                    if isfield(LoadedStruct.m_MRI_Data.Surfaces,'CortexHiResVcolor')
                         oMRI = set_CortexHiResVcolor(oMRI,LoadedStruct.m_MRI_Data.Surfaces.CortexHiResVcolor);
                    else
                         oMRI = set_CortexHiResVcolor(oMRI,zeros(size(LoadedStruct.m_MRI_Data.Surfaces.CortexHiRes.VertexBuffer,1),1));
                    end
                end
            end
            
            if( isfield(LoadedStruct.m_MRI_Data.Surfaces, 'CortexLowRes') )         
                if( isfield(LoadedStruct.m_MRI_Data.Surfaces.CortexLowRes, 'VertexBuffer') && ...
                    isfield(LoadedStruct.m_MRI_Data.Surfaces.CortexLowRes, 'IndexBuffer') )
                    oMRI = set_CortexMeshLowRes( oMRI, LoadedStruct.m_MRI_Data.Surfaces.CortexLowRes.VertexBuffer, ...
                                                       LoadedStruct.m_MRI_Data.Surfaces.CortexLowRes.IndexBuffer );
                    if isfield(LoadedStruct.m_MRI_Data.Surfaces.CortexLowRes,'Vcolor');
                         oMRI = set_CortexLowResVcolor(oMRI,LoadedStruct.m_MRI_Data.Surfaces.CortexLowRes.Vcolor);
                    else
                         oMRI = set_CortexLowResVcolor(oMRI,zeros(size(LoadedStruct.m_MRI_Data.Surfaces.CortexLowRes.VertexBuffer,1),1));
                    end
                end
            end

            if( isfield(LoadedStruct.m_MRI_Data.Surfaces, 'Skin') )
                
                if( isfield(LoadedStruct.m_MRI_Data.Surfaces.Skin, 'VertexBuffer') && ...
                    isfield(LoadedStruct.m_MRI_Data.Surfaces.Skin, 'IndexBuffer') )
                    oMRI = set_SkinMesh( oMRI, LoadedStruct.m_MRI_Data.Surfaces.Skin.VertexBuffer, ...
                                               LoadedStruct.m_MRI_Data.Surfaces.Skin.IndexBuffer );
                if isfield(LoadedStruct.m_MRI_Data.Surfaces.Skin,'Vcolor')
                     oMRI = set_SkinVcolor(oMRI,LoadedStruct.m_MRI_Data.Surfaces.Skin.Vcolor);
                else
                     oMRI = set_SkinVcolor(oMRI,zeros(size(LoadedStruct.m_MRI_Data.Surfaces.Skin.VertexBuffer,1),1));
                end
                
                end
            end
        end
                    
        oPrjData = set_MRI_Data( oPrjData, oMRI );             
    end
    