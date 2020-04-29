function Draw_DigFiducials( oDispOpt, oInterDlgComm, oHelmet )

    persistent v_hPrevDispItems;
    
    %Nettoyage
    if( ~isempty( v_hPrevDispItems ) )
        delete( v_hPrevDispItems(find( ishandle(v_hPrevDispItems) & (v_hPrevDispItems ~= 0) )));
        v_hPrevDispItems = [];
    end
    
    sMtg = get_Mtg( oHelmet );
	matFids = sMtg.matFiducials;
    
    FontWeight = 'normal';
    if( get_DispOptChecked( oDispOpt, 'App_BoldLabels' ) )
        FontWeight = 'bold';
    end
    vColor = get_ItemColor( oDispOpt, 'Dig_Fiducials' );
        
    if( get_DispOptChecked( oDispOpt, 'Dig_SubjectFiducials' ) && ~isempty( find(matFids) ) )
        
        v_hPrevDispItems(1) = plot3( matFids(1,1), matFids(1,2), matFids(1,3), ...
                                     '.', 'Color', vColor );
        v_hPrevDispItems(2) = text(  matFids(1,1), matFids(1,2), matFids(1,3), ...
                                     ' NAS', 'FontWeight', FontWeight, ...
                                     'Color', vColor );
        v_hPrevDispItems(3) = plot3( matFids(2,1), matFids(2,2), matFids(2,3), ...
                                     '.', 'Color', vColor );
        v_hPrevDispItems(4) = text(  matFids(2,1), matFids(2,2), matFids(2,3), ...
                                     ' LPA', 'FontWeight', FontWeight, ...
                                     'Color', vColor );
        v_hPrevDispItems(5) = plot3( matFids(3,1), matFids(3,2), matFids(3,3), ...
                                     '.', 'Color', vColor );
        v_hPrevDispItems(6) = text(  matFids(3,1), matFids(3,2), matFids(3,3), ...
                                     ' RPA', 'FontWeight', FontWeight, ...
                                     'Color', vColor );
    end