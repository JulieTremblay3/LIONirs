function Draw_HelmetHoles( oDispOpt, oInterDlgComm, oHelmet )

    %Previously displayed items handles
    persistent v_hPrevDispHoles;    %Holes circles
    persistent v_hPrevDispNormals;  %Holes normals
    persistent v_hPrevDispLabels;   %Holes Labels
    
    sHelmetAxeDisp = get_HelmetAxeDisp(oInterDlgComm);
    vHoles = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );
    
    %Initialisation de la matrice de handles d'items affichés précedemment.
    %Cette section n'est reellement utilie qu'au chargement d'un nouveau
    %casque, lorsque le nombre de trous change entre 2 affichage. 
    if ( length(v_hPrevDispHoles) ~= length(vHoles) )
        v_hPrevDispHoles      = zeros( 1, length(vHoles) );
    end
    if ( length(v_hPrevDispNormals) ~= length(vHoles) )
        v_hPrevDispNormals    = zeros( 1, length(vHoles) );
    end
    if ( length(v_hPrevDispLabels) ~= length(vHoles) )
        v_hPrevDispLabels     = zeros( 1, length(vHoles) );
    end
    
    %Validation des handles (mettre a zero les handles invalides, ce qui n'est pas suppose arriver)
    v_hPrevDispHoles      = v_hPrevDispHoles      .* ishandle(v_hPrevDispHoles);
    v_hPrevDispNormals    = v_hPrevDispNormals    .* ishandle(v_hPrevDispNormals);
    v_hPrevDispLabels     = v_hPrevDispLabels     .* ishandle(v_hPrevDispLabels);
    
    %Le vecteur de visibilite prend le meme nombre d'elements que le
    %vecteur de trous
    if( length(v_hPrevDispHoles) ~= length(sHelmetAxeDisp.v_bVisible) )
        sHelmetAxeDisp.v_bVisible = ones(size(v_hPrevDispHoles));
    end
    
    
    if( get_DispOptChecked( oDispOpt, 'Dig_HelmetHoles' ) )
        %Effacer du graphique les elements caches
        v_iHolesToErase = find( v_hPrevDispHoles & ~sHelmetAxeDisp.v_bVisible );
        %Afficher les élements visibles qui ne le sont pas déja.
        v_iHolesToDisplay = find (    [vHoles.Type] == get_PtTypeNo(oHelmet, 'NormalHole' ) ...
                                   & ~v_hPrevDispHoles & sHelmetAxeDisp.v_bVisible ); 
    else
        %Effacer tous du graphique car inactifs
        v_iHolesToErase = find( v_hPrevDispHoles );
        v_iHolesToDisplay = [];
    end
    delete( v_hPrevDispHoles( v_iHolesToErase ) );
    v_hPrevDispHoles( v_iHolesToErase ) = 0;
    
    
    if( get_DispOptChecked( oDispOpt, 'Dig_HelmetHolesNormals' ) )
        %Effacer du graphique les elements caches
        v_iNormalsToErase = find( v_hPrevDispNormals & ~sHelmetAxeDisp.v_bVisible );
        %Afficher les élements visibles qui ne le sont pas déja.
        v_iNormalsToDisplay = find (    [vHoles.Type] == get_PtTypeNo(oHelmet, 'NormalHole' ) ...
                                   & ~v_hPrevDispNormals & sHelmetAxeDisp.v_bVisible );
    else
        %Effacer tous du graphique car inactifs
        v_iNormalsToErase = find( v_hPrevDispNormals );
        v_iNormalsToDisplay = [];
    end
    delete( v_hPrevDispNormals( v_iNormalsToErase ) );
    v_hPrevDispNormals( v_iNormalsToErase ) = 0;
    
    
    if( get_DispOptChecked( oDispOpt, 'Lbl_HelmetHoleId' ) )
        %Effacer du graphique les elements caches ou qui ne respectent pas
        %l'intervalle d'affichage d'identificateurs
        %disp( sprintf( 'numel(v_hPrevDispLabels): %d', numel(v_hPrevDispLabels) ) );
        %disp( sprintf( 'numel(sHelmetAxeDisp.v_bDisplayLabels): %d', numel(sHelmetAxeDisp.v_bDisplayLabels) ) );
        %disp( sprintf( 'numel(sHelmetAxeDisp.v_bVisible): %d', numel(sHelmetAxeDisp.v_bVisible) ) );
        %disp( sprintf( 'numel(vHoles): %d', numel(vHoles) ) );
        v_iLabelsToErase = find( v_hPrevDispLabels & ( ~v_hPrevDispLabels | ~sHelmetAxeDisp.v_bDisplayLabels ) );
        
        %Afficher les élements visibles qui ne le sont pas déja.
        v_iLabelsToDisplay = find (    [vHoles.Type] == get_PtTypeNo(oHelmet, 'NormalHole' ) ...
                                       & ~v_hPrevDispLabels & sHelmetAxeDisp.v_bVisible & sHelmetAxeDisp.v_bDisplayLabels );
    else
        %Effacer tous du graphique car inactifs
        v_iLabelsToErase = find( v_hPrevDispLabels );
        v_iLabelsToDisplay = [];
    end
    delete( v_hPrevDispLabels( v_iLabelsToErase ) );
    v_hPrevDispLabels( v_iLabelsToErase ) = 0;
    
    %Initialisation des elements d'affaichage
    nPoints = 12;
    DiameterSrc = 0.0054;
    DiameterDet = 0.011;
    vCircleSrc = CreateCircle( DiameterSrc, nPoints );
    vCircleDet = CreateCircle( DiameterDet, nPoints );
    vCircleUnknown = CreateCircle( DiameterSrc/2, nPoints );
%     figure
%     plot(vCircleSrc(:,1),vCircleSrc(:,2))
    %Memorisation des options d'affichage
    FontWeight = 'normal';
    if( get_DispOptChecked( oDispOpt, 'App_BoldLabels' ) )
        FontWeight = 'bold';
    end
    LineWidth = 1+get_DispOptChecked( oDispOpt, 'App_BoldLines' );
    vColor_CompleteHelmet = get_ItemColor( oDispOpt, 'Dig_Holes_Complete' );
    vColor_PartialHelmet  = get_ItemColor( oDispOpt, 'Dig_Holes_Partial' );
    vColor_HelmetHoleId   = get_ItemColor( oDispOpt, 'Lbl_HelmetHoleId' );
    
    %Holes
    for( Pos=1:length(v_iHolesToDisplay) )
        
        p = v_iHolesToDisplay(Pos);

        % Application de la matrice de transformation     
        if( vHoles(p).CanBeDet )
            vCircle = vCircleDet*vHoles(p).Transformation;
        elseif( vHoles(p).CanBeSrc )
            vCircle = vCircleSrc*vHoles(p).Transformation;
        else
            vCircle = vCircleUnknown*vHoles(p).Transformation;
        end

        % Determination de la couleur dependemment de la provenance
        if( vHoles(p).IsFromCompleteHelmet )
            vColor = vColor_CompleteHelmet;
        else
            vColor = vColor_PartialHelmet;
        end

        %Affichage
        v_hPrevDispHoles(p) = plot3( vCircle(:,1), ...
                                     vCircle(:,2), ...
                                     vCircle(:,3), ...
                                     'color', vColor, ...
                                     'LineWidth', LineWidth, ...
                                     'Clipping', 'on' );
    end
    
    %Normals
    for( Pos=1:length(v_iNormalsToDisplay) )
        
        p = v_iNormalsToDisplay(Pos);
        
        %Si une normale est definie
        if( vHoles(p).Normal.x || vHoles(p).Normal.y || vHoles(p).Normal.z )
                
            %Prolonger la normale jusqu'au cortex (minimum de 3 cm)
            if( vHoles(p).CortexDepth < 0.03 )
                vNormal = [ 0,   0,  0.00,  1; ...
                            0,   0, -0.03,  1 ];
            else
                vNormal = [ 0,   0,  0.00,                  1; ...
                            0,   0,-vHoles(p).CortexDepth,  1 ];
            end

            if( vHoles(p).IsFromCompleteHelmet )
                vColor = vColor_CompleteHelmet;
            else
                vColor = vColor_PartialHelmet;
            end
            
            vNormal = vNormal*vHoles(p).Transformation;
            v_hPrevDispNormals(p) = plot3( vNormal(:,1), ...
                                           vNormal(:,2), ...
                                           vNormal(:,3), ...
                                           'color', vColor, ...
                                           'LineWidth', LineWidth+1, ...
                                           'Clipping', 'on' );
        end
    end
    
     
    %Labels
    for( Pos=1:length(v_iLabelsToDisplay) )
        
        p = v_iLabelsToDisplay(Pos);
        
        v_hPrevDispLabels(p) =  text(  vHoles(p).Coord.x, ...
                                       vHoles(p).Coord.y, ...
                                       vHoles(p).Coord.z, ...
                                       sprintf( ' %s',vHoles(p).Label ), ...
                                       'FontWeight', 'normal', ...
                                       'FontSize', 9, ...
                                       'Color', vColor_HelmetHoleId, ...
                                       'HorizontalAlignment', 'center');
                                      
                                                 
    end
    if ( get_DispOptChecked( oDispOpt, 'line_holes ' ) ) 
        %pour tout les droit
        pD = [];
        pG = [];
        for p = 1:numel(vHoles)   
            if vHoles(p).Label(1)=='D';
                pD = [pD,p];
            end 
            if vHoles(p).Label(1)=='G';
                pG = [pG,p];
            end
        end
         for ind_lettre = 1:10
                lettre = char(64+ind_lettre);
                 listP = [];
                 Nb = [];
               for i = 1:numel(pD)                           
                p = pD(i);                  
                if strcmp(vHoles(p).Label(2),lettre)
                    listP = [listP, p];
                    Nb = [Nb, str2num(vHoles(p).Label(3:end))];
                end 
                if ~isempty(listP)
                    [result,ind] = sort(Nb);
                    listP = listP(ind);% 
                    %Dessine la ligne
                    ligne.x = [];
                    ligne.y = [];
                    ligne.z = [];
                    for i = 1:numel(listP)
                            ligne.x = [ligne.x, vHoles(listP(i)).Coord.x];
                            ligne.y = [ligne.y, vHoles(listP(i)).Coord.y];
                            ligne.z = [ligne.z, vHoles(listP(i)).Coord.z];
                    end
                    h = plot3(ligne.x,ligne.y,ligne.z,'LineWidth',1);
                    if mod(ind_lettre,2)
                        set(h,'color','m')
                    else
                        set(h,'color','b')
                    end
                end
                
               end
               listP = [];
               Nb = [];
               for i = 1:numel(pG)
                    p = pG(i);
                if strcmp(vHoles(p).Label(2),lettre)
                    listP = [listP, p];
                    Nb = [Nb, str2num(vHoles(p).Label(3:end))];
                end 
                if ~isempty(listP)
                    [result,ind] = sort(Nb);
                    listP = listP(ind);% 
                    %Dessine la ligne
                    ligne.x = [];
                    ligne.y = [];
                    ligne.z = [];
                    for i = 1:numel(listP)
                            ligne.x = [ligne.x, vHoles(listP(i)).Coord.x];
                            ligne.y = [ligne.y, vHoles(listP(i)).Coord.y];
                            ligne.z = [ligne.z, vHoles(listP(i)).Coord.z];
                    end
                    h = plot3(ligne.x,ligne.y,ligne.z,'LineWidth',1);
                    if mod(ind_lettre,2)
                        set(h,'color','m')
                    else
                        set(h,'color','b')
                    end
                end
                
               end
    end
  %  ligne entre les trous qui se suivent
    end