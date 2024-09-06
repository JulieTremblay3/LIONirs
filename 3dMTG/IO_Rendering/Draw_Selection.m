function Draw_Selection( oDispOpt, oInterDlgComm, oHelmet )
    
    persistent v_hPrevDispSel;
    persistent v_pPrevDispSel;
    
    v_pSelection = get_HelmetAxeDisp_vSelection(oInterDlgComm);
    vHoles = get_vHoles( oHelmet );
    
    % Uniformisation des tailles de vecteurs
    NbActuItems = length(v_pSelection);
    NbPrevItems = length(v_pPrevDispSel);
    MaxLength = max(NbActuItems, NbPrevItems);
    v_pPrevDispSel(NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    v_hPrevDispSel(NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    v_pSelection(  NbActuItems+1:MaxLength) = 0;  %Padding a zero

    %Sortir s'il n'y a rien a faire
	if( ~MaxLength )
        return;
    end
    
    %Validation des handles
    v_hPrevDispSel = v_hPrevDispSel .* ishandle(v_hPrevDispSel);
    v_pPrevDispSel = v_pPrevDispSel .* (v_hPrevDispSel ~= 0);
    
    %Creation du vecteur de visibilite
    v_bVisible = zeros(1,MaxLength);
    v_Ind = find(v_pSelection);
    HelmetAxeDisp = get_HelmetAxeDisp(oInterDlgComm);
    
    v_bVisible(v_Ind) = HelmetAxeDisp.v_bVisible(v_pSelection(v_Ind));

    %Effacer du graphique les élements affichés précèdemment qui
    %sont soit maintenant non-visibles ou différents des nouveaux elems.
    v_iToErase = find( v_pPrevDispSel & ( ~v_bVisible | ( v_pPrevDispSel ~= v_pSelection ) ) );
    delete( v_hPrevDispSel( v_iToErase ) );
    v_hPrevDispSel( v_iToErase ) = 0;
    v_pPrevDispSel( v_iToErase ) = 0;
    
    vColor = get_ItemColor( oDispOpt, 'Selection' );
    
    
    %Afficher les élements visibles différents des anciens elems.
    v_iToDisplay = find ( v_bVisible & ( v_pPrevDispSel ~= v_pSelection ) );
    v_pPrevDispSel( v_iToDisplay ) = v_pSelection( v_iToDisplay );                                                  
    for( Pos=1:length(v_iToDisplay) )
        v_hPrevDispSel( v_iToDisplay(Pos) ) = plot3(  vHoles(v_pSelection(v_iToDisplay(Pos))).Coord.x, ...
                                                      vHoles(v_pSelection(v_iToDisplay(Pos))).Coord.y, ...
                                                      vHoles(v_pSelection(v_iToDisplay(Pos))).Coord.z, ...
                                                      'x', 'color', vColor, ...
                                                      'LineWidth', 2 );
    end
    