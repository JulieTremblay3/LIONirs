function Draw_ItemInformation( oDispOpt, oInterDlgComm, oHelmet )
   
    persistent hLastId;
    persistent pLastId;
    sHelmetAxeDisp = get_HelmetAxeDisp(oInterDlgComm);
    
    %Validation du handle de l'item precedemment trace dans le graphique.
    if( isempty(hLastId) )
        hLastId = 0;
        pLastId = 0;
    elseif( ~ishandle( hLastId ) )
        hLastId = 0;
        pLastId = 0;
    end
    
    if( (hLastId ~= 0) && ishandle(hLastId) ) % Must check explicitely that the handle is not zero (handle cannot be used directly as a boolean).
        delete( hLastId );
        hLastId = 0;
        pLastId = 0;
    end

    FontWeight = 'normal';
    if( get_DispOptChecked( oDispOpt, 'App_BoldLabels' ) )
        FontWeight = 'bold';
    end
    vColor = get_ItemColor( oDispOpt, 'Identification' );
    
    if( sHelmetAxeDisp.pIdentification )
        pItem = sHelmetAxeDisp.pIdentification;
        
        vHoles = get_vHoles( oHelmet );
        
         
        if vHoles(pItem).Coord.y < 0
        ItemString = sprintf( '%s\\rightarrow ', vHoles(pItem).Label );
        hLastId = text( vHoles(pItem).Coord.x, vHoles(pItem).Coord.y, vHoles(pItem).Coord.z, ...
                        ItemString, 'FontWeight', FontWeight,'FontSize', 14, ...
                        'color', vColor, 'HorizontalAlignment', 'right' );
        else
        ItemString = sprintf( '\\leftarrow%s ', vHoles(pItem).Label);
        hLastId = text( vHoles(pItem).Coord.x, vHoles(pItem).Coord.y, vHoles(pItem).Coord.z, ...
                 ItemString, 'FontWeight', FontWeight,'FontSize', 14, ...
                 'color', vColor, 'HorizontalAlignment', 'left' );
        end
        pLastId = pItem;
    end