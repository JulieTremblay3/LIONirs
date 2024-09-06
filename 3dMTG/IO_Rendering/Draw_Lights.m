function Draw_Lights( oDispOpt )

    persistent v_hPrevDispItems;
    
    %Effacer les anciens elements graphiques avant d'afficher les nouveaux
    if( ~isempty( v_hPrevDispItems ) )
        delete( v_hPrevDispItems(find( ishandle(v_hPrevDispItems) & (v_hPrevDispItems ~= 0) )));
        v_hPrevDispItems = [];
    end
    
    %Si l'option d'affichage des reflexions est selectionnee
    %(L'option de reflexions n'est en realite que 3 lumieres)
    if( get_DispOptChecked( oDispOpt, 'App_Reflections' ) )
        v_hPrevDispItems(1) = light( 'Position', [500,500,500] );
        v_hPrevDispItems(2) = light( 'Position', [500,-500,500] );
        v_hPrevDispItems(3) = light( 'Position', [-500,0,0] );
    end