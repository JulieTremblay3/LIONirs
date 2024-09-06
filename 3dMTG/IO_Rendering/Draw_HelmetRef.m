function Draw_HelmetRef( oDispOpt, oInterDlgComm, oHelmet  )

    persistent v_hPrevDispItems;
    
    %Nettoyage
    if( ~isempty( v_hPrevDispItems ) )
        delete( v_hPrevDispItems(find( ishandle(v_hPrevDispItems) & (v_hPrevDispItems ~= 0) )));
        v_hPrevDispItems = [];
    end
    
    FontWeight = 'normal';
    if( get_DispOptChecked( oDispOpt, 'App_BoldLabels' ) )
        FontWeight = 'bold';
    end
    vColor = get_ItemColor( oDispOpt, 'Dig_HelmetRef' );
    
    if( get_DispOptChecked( oDispOpt, 'Dig_HelmetRegistrationPoints' ) )
        
        v_Holes = get_vHoles( oHelmet );
        
        v_iToDisplay = find( [v_Holes.Type] == get_PtTypeNo(oHelmet, 'HelmetRef' ) );
                      
        for( iDisp=1:numel(v_iToDisplay) )
            p = v_iToDisplay(iDisp);
            
            v_hPrevDispItems((iDisp-1)*2+1) = plot3( v_Holes(p).Coord.x, v_Holes(p).Coord.y, v_Holes(p).Coord.z, ...
                                                     '.', 'Color', vColor );
            v_hPrevDispItems((iDisp-1)*2+2) = text(  v_Holes(p).Coord.x, v_Holes(p).Coord.y, v_Holes(p).Coord.z, ...
                                                     sprintf( ' %s',v_Holes(p).Label ), 'FontWeight', FontWeight, ...
                                                     'FontSize', 14, 'Color', vColor );
        end
    end