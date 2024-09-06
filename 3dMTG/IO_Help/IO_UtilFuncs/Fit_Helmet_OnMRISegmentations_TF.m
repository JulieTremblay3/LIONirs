% Fonction permettant de calculer la profondeur des trous du casque qui
% fait en sorte que les optodes arrivent sur la peau (intersection de la
% segmentation et de la peau).
%
%
function oHelmet = Fit_Helmet_OnMRISegmentations_TF( oHelmet, oMRI_Data, bProcessSkin, bProcessCortex, bDisplay,mat_transform )

    vHoles = get_vHoles(oHelmet);
    
    %Fonction de transfert des coordonnees de trous en voxel-space du MRI
    [HCoords, HNorms ] = CoRegisterHelmetToMRI_TF( oHelmet, oMRI_Data,1,mat_transform);
    
    %Memorisation des matrices de segmentations
    matSkinSegmentation = get_SkinSegmentation( oMRI_Data );
    matCortexSegmentation = get_CortexSegmentation( oMRI_Data );
    [nVoxSkinSegX,nVoxSkinSegY,nVoxSkinSegZ] = size( matSkinSegmentation );
    [nVoxCortSegX,nVoxCortSegY,nVoxCortSegZ] = size( matSkinSegmentation );
    
    vPhysicalVoxDim = get_PhysicalVoxDim( oMRI_Data );
    if( numel(find( vPhysicalVoxDim == vPhysicalVoxDim(1) ) ) < 3 )
        warndlg( 'Voxel size is not cubic: possible errors in distances calculated on skin segmentation (ERR204)');
    end
    %MaxDepth_m = 0.05; %in meters (5 cm)
    MaxDepth_m = 0.10;
    if vPhysicalVoxDim == 0
        MaxDepth_m = 0.10;
        MaxDepth_vox = 1;
    else
        MaxDepth_vox = ceil(MaxDepth_m/mean(vPhysicalVoxDim));
    end
    %Structure de memorisation des profondeurs de trous.
    %vHoles_SkinDepth = zeros( 1,size(HCoords,1) );
    if( exist( 'bDisplay' ) && bDisplay )
        hold on;
        [VertexBuffer, IndexBuffer] = get_CortexMeshHiRes( oMRI_Data );
        patch('Vertices',VertexBuffer,'Faces',IndexBuffer, ...
              'EdgeColor', 'none', 'FaceColor', [0.5,0.5,0.5], ...
              'FaceLighting', 'gouraud', ...
              'SpecularColorReflectance', 0.5, ...
              'DiffuseStrength', 0.5, ...
              'SpecularStrength', 1);
        light( 'Position', [500,500,500] );
        light( 'Position', [500,-500,500] );
        light( 'Position', [-500,0,0] );
    end
    
    %Pour chaque position de trou ce casque (en voxel-space)
    for( pElem = 1:size(HCoords,1) )
        
        DepthSteps = 0:1:MaxDepth_vox;
        HCoordsStepsX = int16(HCoords(pElem,1)-HNorms(pElem,1)*DepthSteps);
        HCoordsStepsY = int16(HCoords(pElem,2)-HNorms(pElem,2)*DepthSteps);
        HCoordsStepsZ = int16(HCoords(pElem,3)-HNorms(pElem,3)*DepthSteps);
        
        %hold on;
        if( exist( 'bDisplay' ) && bDisplay )
            plot3( HCoordsStepsX, HCoordsStepsY,HCoordsStepsZ, 'r' );
        end
        
        %----------------------------------------------------------------------------------------------------------
        %-- SKIN  -------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------
        if( bProcessSkin && ~isempty(matSkinSegmentation) )
            SkinReached = false;
            SkinDepth = 0;
            %Pour au moins chaque voxel dans la trajectoire de la normale
            for( iDepth=1:numel(DepthSteps) )

                %Si a l'interieur de la limite de voxels (-1 et +1 pour cube)
                if(   ( HCoordsStepsX(iDepth) > 1 && HCoordsStepsX(iDepth) < nVoxSkinSegX ) ...
                    &&( HCoordsStepsY(iDepth) > 1 && HCoordsStepsY(iDepth) < nVoxSkinSegY ) ...
                    &&( HCoordsStepsZ(iDepth) > 1 && HCoordsStepsZ(iDepth) < nVoxSkinSegZ ) )

                    %Verifier que le nombre de voxels du cube centre en
                    %HCoordsStepsX/Y/Z faisant partie de la segmentation de la
                    %peau est d'au moins 13 (13/27 des voxels du cube)
                    if( numel( find( matSkinSegmentation(HCoordsStepsX(iDepth)-1:HCoordsStepsX(iDepth)+1, ...
                                                         HCoordsStepsY(iDepth)-1:HCoordsStepsY(iDepth)+1, ...
                                                         HCoordsStepsZ(iDepth)-1:HCoordsStepsZ(iDepth)+1 ) ) ) >= 13 )
                        %Conversion de profondeur vox->m (distance euclidienne
                        %respectee)
                        SkinDepth = ((iDepth-1)/MaxDepth_vox)*MaxDepth_m;
                        SkinReached = true;
                        if( exist( 'bDisplay' ) && bDisplay )
                            plot3( [HCoordsStepsX(iDepth),HCoordsStepsX(iDepth)-1], ...
                                   [HCoordsStepsY(iDepth),HCoordsStepsY(iDepth)-1], ...
                                   [HCoordsStepsZ(iDepth),HCoordsStepsZ(iDepth)-1], ...
                                   'color', [1,1,1], 'linewidth', 3 );
                        end
                        break;
                    end
                end
            end
            vHoles(pElem).SkinDepth = SkinDepth;
        end
        
        %----------------------------------------------------------------------------------------------------------
        %-- CORTEX  -------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------
        if( bProcessCortex && ~isempty(matCortexSegmentation) )
            CortexReached = false;
            CortexDepth = 0;
            %Pour au moins chaque voxel dans la trajectoire de la normale
            for( iDepth=1:numel(DepthSteps) )

                %Si a l'interieur de la limite de voxels (-1 et +1 pour cube)
                if(   ( HCoordsStepsX(iDepth) > 1 && HCoordsStepsX(iDepth) < nVoxCortSegX ) ...
                    &&( HCoordsStepsY(iDepth) > 1 && HCoordsStepsY(iDepth) < nVoxCortSegY ) ...
                    &&( HCoordsStepsZ(iDepth) > 1 && HCoordsStepsZ(iDepth) < nVoxCortSegZ ) )

                    %Verifier que le nombre de voxels du cube centre en
                    %HCoordsStepsX/Y/Z faisant partie de la segmentation de la
                    %peau est d'au moins 13 (13/27 des voxels du cube)
                    if( numel( find( matCortexSegmentation(HCoordsStepsX(iDepth)-1:HCoordsStepsX(iDepth)+1, ...
                                                           HCoordsStepsY(iDepth)-1:HCoordsStepsY(iDepth)+1, ...
                                                           HCoordsStepsZ(iDepth)-1:HCoordsStepsZ(iDepth)+1 ) ) ) >= 13 )
                        %Conversion de profondeur vox->m (distance euclidienne respectee)
                        CortexDepth = ((iDepth-1)/MaxDepth_vox)*MaxDepth_m;
                        CortexReached = true;
                        if( exist( 'bDisplay' ) && bDisplay )
                            plot3( [HCoordsStepsX(iDepth),HCoordsStepsX(iDepth)-1], ...
                                   [HCoordsStepsY(iDepth),HCoordsStepsY(iDepth)-1], ...
                                   [HCoordsStepsZ(iDepth),HCoordsStepsZ(iDepth)-1], ...
                                   'color', [1,1,1], 'linewidth', 3 );
                        end
                        break;
                    end
                end
            end
            vHoles(pElem).CortexDepth = CortexDepth;
        end
        
    end
    
    oHelmet = set_vHoles( oHelmet, vHoles );
    
