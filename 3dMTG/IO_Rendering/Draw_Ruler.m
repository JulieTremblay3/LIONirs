function Draw_Ruler( oDispOpt, oInterDlgComm, oHelmet )

    persistent v_hPrevDispItems;
    
    %Effacer les anciens elements graphiques avant d'afficher les nouveaux
    if( ~isempty( v_hPrevDispItems ) )
        delete( v_hPrevDispItems( find( ishandle(v_hPrevDispItems) & (v_hPrevDispItems ~= 0) )));
        v_hPrevDispItems = [];
    end
    
    %Ne pas afficher la regle pendant la rotation de camera
    if( get_HelmetAxeDisp_AxeRotating( oInterDlgComm ) )
        return;
    end
    
    sHelmetAxeDisp = get_HelmetAxeDisp(oInterDlgComm);
    p1 = sHelmetAxeDisp.pRulerFirst;
    p2 = sHelmetAxeDisp.pRulerLast;
    vHoles = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );
    
    %Memorisation des preferences d'affichage 
    LineWidth = 2+get_DispOptChecked( oDispOpt, 'App_BoldLines' );
    FontWeight = 'normal';
    if( get_DispOptChecked( oDispOpt, 'App_BoldLabels' ) )
        FontWeight = 'bold';
    end
    vColor = get_ItemColor( oDispOpt, 'Ruler' );
    
        
    %Si le premier point de la regle existe
    if( p1 )
        %Si ce point est visible de la camera
        if( sHelmetAxeDisp.v_bVisible(p1) )
            NormalP1 = [vHoles(p1).Normal.x, vHoles(p1).Normal.y, vHoles(p1).Normal.z];
            CoordP1 = [vHoles(p1).Coord.x, vHoles(p1).Coord.y, vHoles(p1).Coord.z];
            CoordUnderP1 = CoordP1-NormalP1*vHoles(p1).SkinDepth;
            
            v_hPrevDispItems(1) = plot3( CoordP1(1), ...
                                       CoordP1(2), ...
                                       CoordP1(3), ...
                                       'x', 'color', ...
                                       vColor,...
                                       'LineWidth', 2 );

            v_hPrevDispItems(2) = plot3( [CoordUnderP1(1), CoordP1(1)], ...
                                       [CoordUnderP1(2), CoordP1(2)], ...
                                       [CoordUnderP1(3), CoordP1(3)], ...
                                       'color', ... %'x', 'color', ...
                                       vColor,...
                                       'LineWidth', LineWidth );
        end
    end
    
    %Si le premier point de la regle existe
    if( p2 )
        %Si ce point est visible de la camera
        if( sHelmetAxeDisp.v_bVisible(p2) )
            NormalP2 = [vHoles(p2).Normal.x, vHoles(p2).Normal.y, vHoles(p2).Normal.z];
            CoordP2 = [vHoles(p2).Coord.x, vHoles(p2).Coord.y, vHoles(p2).Coord.z];
            CoordUnderP2 = CoordP2-NormalP2*vHoles(p2).SkinDepth;

            v_hPrevDispItems(3) = plot3( CoordP2(1), ...
                                       CoordP2(2), ...
                                       CoordP2(3), ...
                                       'x', 'color', ...
                                       vColor,...
                                       'LineWidth', 2 );

            v_hPrevDispItems(4) = plot3( [CoordUnderP2(1), CoordP2(1)], ...
                                       [CoordUnderP2(2), CoordP2(2)], ...
                                       [CoordUnderP2(3), CoordP2(3)], ...
                                       'color', ...%'x', 'color', ...
                                       vColor,...
                                       'LineWidth', LineWidth );
        end
    end
    
    %Si les 2 points existent et sont visibles
    if( p1 && p2 )
        if( sHelmetAxeDisp.v_bVisible(p1) && sHelmetAxeDisp.v_bVisible(p2) )
            
            v_hPrevDispItems(5) = plot3( [CoordUnderP1(1); CoordUnderP2(1)], ...
                                       [CoordUnderP1(2); CoordUnderP2(2)], ...
                                       [CoordUnderP1(3); CoordUnderP2(3)], ...
                                       '--', ...,
                                       'color', ...
                                       vColor,...
                                       'LineWidth', LineWidth );
                                     
            dist_cm = get_DistHoles( oHelmet, p1, p2 )*100;
            
            v_hPrevDispItems(6) = text( 0.68,0.04, sprintf( '%2.3f cm', dist_cm ), ...
                                      'Units', 'normalized', ...
                                      'Color', vColor,...
                                      'FontSize', 14, 'FontWeight', FontWeight );
        end
    end