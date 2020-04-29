function Draw_MuxChannels( oDispOpt, oInterDlgComm, oHelmet )

    persistent v_hPrevDispItems;
    
    if( ~isempty(v_hPrevDispItems) )
        delete( v_hPrevDispItems( find( v_hPrevDispItems ~= 0 & ishandle(v_hPrevDispItems) ) ) );
    end
        
    sHelmetAxeDisp = get_HelmetAxeDisp( oInterDlgComm );
    
    sMtg = get_Mtg( oHelmet );
    vHoles = get_vHoles( oHelmet );
    dMaxMx=sMtg.Gen_Params.Max_Dist_Mux;
    dMinMx=sMtg.Gen_Params.Min_Dist_Mux;
    dCntMx=sMtg.Gen_Params.DistContamination;
    dRngMx=0.02;
    v_pDetsFnd = sMtg.v_pDet(find(sMtg.v_pDet));
    v_pSrcsFnd = sMtg.v_pSrc(find(sMtg.v_pSrc));
    nMtgDets = length(v_pDetsFnd);
    nMtgSrcs = length(v_pSrcsFnd);
    nMuxPerDet = floor(sMtg.Gen_Params.NbSrcPerBank / sMtg.Gen_Params.Nb_Longueur_Onde );
    
    if( nMtgSrcs < nMuxPerDet )
        nMuxPerDet = nMtgSrcs;
    end

    if( ~nMtgDets || ~nMtgSrcs )
        return;
    end
    
    %Matrice de coordonnees de source du mtg en vecteurs colonnes:
    mat_SrcCoords = zeros(3, nMtgSrcs);
    for( iSrc=1:nMtgSrcs )
         mat_SrcCoords( :, iSrc )   = [ vHoles(v_pSrcsFnd(iSrc)).Coord.x;
                                        vHoles(v_pSrcsFnd(iSrc)).Coord.y;
                                        vHoles(v_pSrcsFnd(iSrc)).Coord.z ];
    end
    
    %** A REVOIR POUR UNE CONFIG QUI NE POSSEDE PAS TOUTES LES PERIODES **
    LineWidth = 1+get_DispOptChecked( oDispOpt, 'App_BoldLines' );
    
    if( get_DispOptChecked( oDispOpt, 'Opt_MuxChannel' ) )
        for( iDet=1:nMtgDets )
            pDet = v_pDetsFnd(iDet);
            
            matDists = [];
            for( iSrc=1:size(mat_SrcCoords,2) )
                matDists(iSrc) = get_DistHoles( oHelmet, v_pSrcsFnd(iSrc), pDet );
            end

            %Matrice: [c1 = pSrcs, c2 = (longueur mux)^2 ], en ordre de distance.
            %mat_Sorted_Srcs = sortrows( [ v_pSrcsFnd', [[[1,1,1]*[matDiff.^2]]'].^0.5 ], 2);
            mat_Sorted_Srcs = sortrows( [ v_pSrcsFnd', matDists' ], 2);
            dClosest = mat_Sorted_Srcs(1,2);
            
            for( iSrc=1:nMuxPerDet )
                pSrc = mat_Sorted_Srcs(iSrc,1);
                dSrc = mat_Sorted_Srcs(iSrc,2);

                if( dSrc >= dCntMx )
                    break;
                end
                
                if( dSrc <= dClosest + dRngMx )
                    vColor = hsv2rgb( [ (0.5-((dSrc-dClosest)/dRngMx)*0.5), 1.0, 1.0 ] );
                else
                    vColor = hsv2rgb( [ 0.0, 1, 1 ] ); %Red
                end

                if( sHelmetAxeDisp.v_bVisible(pDet) && sHelmetAxeDisp.v_bVisible(pSrc) )
                    
                    %Trajet
                    v_hPrevDispItems(numel(v_hPrevDispItems)+1) = plot3([vHoles(pDet).Coord.x-0.02*vHoles(pDet).Normal.x,vHoles(pSrc).Coord.x-vHoles(pDet).SkinDepth*vHoles(pSrc).Normal.x], ...
                                                                        [vHoles(pDet).Coord.y-0.02*vHoles(pDet).Normal.y,vHoles(pSrc).Coord.y-vHoles(pDet).SkinDepth*vHoles(pSrc).Normal.y], ...
                                                                        [vHoles(pDet).Coord.z-0.02*vHoles(pDet).Normal.z,vHoles(pSrc).Coord.z-vHoles(pDet).SkinDepth*vHoles(pSrc).Normal.z], 'Color', vColor, 'LineWidth', LineWidth );
                end
            end
        end
    end