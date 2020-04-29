function Draw_CoordSysAxes( oDispOpt, oInterDlgComm )

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
    vColor = get_ItemColor( oDispOpt, 'Dig_CoordinateSystemAxes' );
    LineWidth = 1+get_DispOptChecked( oDispOpt, 'App_BoldLines' );
    
    if( get_DispOptChecked( oDispOpt, 'Dig_CoordinateSystemAxes' ) )

        v_hPrevDispItems(1) = plot3( [0;0.04], [0;0], [0;0], 'Color', vColor, 'linewidth', LineWidth );
        v_hPrevDispItems(2) = plot3( [0;0], [0;0.04], [0;0], 'Color', vColor, 'linewidth', LineWidth );
        v_hPrevDispItems(3) = plot3( [0;0], [0;0], [0;0.04], 'Color', vColor, 'linewidth', LineWidth );
        
        v_hPrevDispItems(4) = text( [0.04], [0], [0], 'x', 'Color', vColor, ...
                                                       'FontWeight', FontWeight, 'FontSize', 10 );
        v_hPrevDispItems(5) = text( [0], [0.04], [0], 'y', 'Color', vColor, ...
                                                       'FontWeight', FontWeight, 'FontSize', 10 );
        v_hPrevDispItems(6) = text( [0], [0], [0.04], 'z', 'Color', vColor, ...
                                                       'FontWeight', FontWeight, 'FontSize', 10 );
    end