function Draw_Optodes( oDispOpt, oInterDlgComm, oHelmet )

    persistent v_hPrevDispItems; %Vecteur d'handles de chaque elements pour chaque position
    persistent v_pPrevDispItems;%Vecteur de position (position=indice lineaire de trou dans vHoles)

    sMtg          = get_Mtg( oHelmet );
    vHoles        = get_vHoles( oHelmet );
    sHelmetAxeDisp = get_HelmetAxeDisp( oInterDlgComm );
    v_pActuDispItems = [sMtg.v_pDet, sMtg.v_pSrc];

    %Uniformisation des tailles afin de permettre les algorithmes
    %vectoriels. Le padding a zero n'influence pas les resultats.
    NbActuItems = length(v_pActuDispItems); %Nombre de positions actuellement visibles
    NbPrevItems = length(v_pPrevDispItems); %Nombre de positions precedemment visibles (lors du dernier appel la fonction)
    MaxLength = max(NbActuItems, NbPrevItems);
    v_pActuDispItems(NbActuItems+1:MaxLength) = 0;  %Padding a zero
    v_pPrevDispItems(NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    v_hPrevDispItems.CortexIntersection(NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    v_hPrevDispItems.LabelsHoleId(      NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    v_hPrevDispItems.LabelsFibersId(    NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    v_hPrevDispItems.Optode(            NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    %Sortir s'il n'y a rien a faire (rien a effacer, rien a afficher)
	if( ~MaxLength )
        return;
    end
    
    %Validation des handles: 
    v_hPrevDispItems.CortexIntersection = v_hPrevDispItems.CortexIntersection .* ishandle(v_hPrevDispItems.CortexIntersection);
    v_hPrevDispItems.LabelsHoleId       = v_hPrevDispItems.LabelsHoleId       .* ishandle(v_hPrevDispItems.LabelsHoleId);
    v_hPrevDispItems.LabelsFibersId     = v_hPrevDispItems.LabelsFibersId     .* ishandle(v_hPrevDispItems.LabelsFibersId);
    v_hPrevDispItems.Optode             = v_hPrevDispItems.Optode             .* ishandle(v_hPrevDispItems.Optode);
    %Les handles d'items affiches precedemment doivent etre non nuls,
    %autrement on ne peut considerer que les items sont encore dans le
    %graphique. Ils ont peut-etre ete effaces entre 2 appels de la fonction
    %d'affichage, avec la commande 'cla' par exemple. Ceci (ci-bas) permet donc de
    %verifier afin de reafficher les items effaces.
    v_pPrevDispItems = v_pPrevDispItems .* (   (v_hPrevDispItems.CortexIntersection ~= 0) ...
                                             | (v_hPrevDispItems.LabelsHoleId ~= 0) ...
                                             | (v_hPrevDispItems.LabelsFibersId ~= 0) ...
                                             | (v_hPrevDispItems.Optode ~= 0) );
    
    %Creation du vecteur de visibilite (vecteur de booleens pour savoir si
    %les positions sont vus par dessus [visibles], soit par en-dessous
    %[invisibles]. Le calcul de visibilite est fait dans une autre
    %fonction (UpdateHelmetVisibility.m de oInterDlgComm).
    v_bVisible = zeros(1,MaxLength);
    v_Ind = find(v_pActuDispItems);
    v_bVisible(v_Ind) = sHelmetAxeDisp.v_bVisible(v_pActuDispItems(v_Ind));

    %Effacer du graphique les élements affichés précèdemment qui
    %sont soit maintenant non-visibles ou différents des nouveaux elems.
    v_iToErase = find( v_pPrevDispItems & ( ~v_bVisible | ( v_pPrevDispItems ~= v_pActuDispItems ) ) );
    
    if( get_DispOptChecked( oDispOpt, 'Opt_CortexIntersectionMarker' ) )
        delete( v_hPrevDispItems.CortexIntersection( v_iToErase( find(v_hPrevDispItems.CortexIntersection(v_iToErase)) ) ) ); %Delete sur elems non-nuls seulement (le 'find' sert a retirer les elements nuls)
        v_hPrevDispItems.CortexIntersection( v_iToErase ) = 0;
    else
        %Delete de tous les Handles valides (les Handles valides sont des
        %elements encore affiches) si cet option d'affichage est desactive
        delete( v_hPrevDispItems.CortexIntersection( find( v_hPrevDispItems.CortexIntersection ) ) );
        v_hPrevDispItems.CortexIntersection(         find( v_hPrevDispItems.CortexIntersection ) ) = 0;
    end
    
    if( get_DispOptChecked( oDispOpt, 'Lbl_OptodesHoleId' ) )
        delete( v_hPrevDispItems.LabelsHoleId(       v_iToErase( find(v_hPrevDispItems.LabelsHoleId(      v_iToErase)) ) ) ); %Delete sur elems non-nuls seulement (le 'find' sert a retirer les elements nuls)
        v_hPrevDispItems.LabelsHoleId( v_iToErase ) = 0;
    else
        %Delete de tous les Handles valides (les Handles valides sont des
        %elements encore affiches) si cet option d'affichage est desactive
        delete( v_hPrevDispItems.LabelsHoleId( find( v_hPrevDispItems.LabelsHoleId ) ) );
        v_hPrevDispItems.LabelsHoleId(         find( v_hPrevDispItems.LabelsHoleId ) ) = 0;
    end
    
    %Les labels de fibres doivent tous etre redessines car pour une meme 
    %position les fibres peuvent changer entre 2 appels d'affichage. Donc,
    %on les efface toutes.
    delete( v_hPrevDispItems.LabelsFibersId( find( v_hPrevDispItems.LabelsFibersId ) ) );
    v_hPrevDispItems.LabelsFibersId(         find( v_hPrevDispItems.LabelsFibersId ) ) = 0;
    
    %Les optodes sont toujours actives. Pas besoin de verifier les options.
    delete( v_hPrevDispItems.Optode(             v_iToErase( find(v_hPrevDispItems.Optode(            v_iToErase)) ) ) ); %Delete sur elems non-nuls seulement (le 'find' sert a retirer les elements nuls)
    v_hPrevDispItems.Optode(             v_iToErase ) = 0;
    v_pPrevDispItems( v_iToErase ) = 0;
    
    %Determiner les positions a afficher: 
    v_iToDisplay = find ( v_bVisible & ( v_pPrevDispItems ~= v_pActuDispItems ) );
    v_pPrevDispItems( v_iToDisplay ) = v_pActuDispItems( v_iToDisplay );

    %OptIndexBuffer est le meme pour les 2 fonctions, seul le diametre
    %change.
    if( get_DispOptChecked( oDispOpt, 'Opt_FittedOnSkin' ) )
        [OptSrc_VertexBuffer, OptIndexBuffer] = CreateCylinderMesh( 0.0054, 1, 12 );
        [OptDet_VertexBuffer, OptIndexBuffer] = CreateCylinderMesh( 0.0110, 1, 12 );
    else
        [OptSrc_VertexBuffer, OptIndexBuffer] = CreateCircleMesh( 0.0054, 12 );
        [OptDet_VertexBuffer, OptIndexBuffer] = CreateCircleMesh( 0.0110, 12 );
    end
    
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    %Affichage des Cercles ou des cylindres de sources et de detecteurs
    for( Pos=1:length(v_iToDisplay) )
        
        p = v_pActuDispItems(v_iToDisplay(Pos));
       
        %Detecteur
        if( find( sMtg.v_pDet == p ) )
            OptVertexBuffer = OptDet_VertexBuffer;
        %Source
        else
            OptVertexBuffer = OptSrc_VertexBuffer;
        end
         
            %Detecteur
            if( find( sMtg.v_pDet == p ) )
                %vColor = get_ItemColor( oDispOpt, 'Opt_LabelDetHole' );
                switch sMtg.Gen_Params.AcqSystem
                    case 'ISS'
                     vColor = [255/255, 40/255,40/255];
                    case  'IMAGINC' %or NIRx nomenclature
                     vColor = [0,0,255/255]; %Bleu
                     case  'NIRx' %or NIRx nomenclature
                      vColor = [0,0,255/255]; %Bleu
                end

                %Source
            else 
                switch sMtg.Gen_Params.AcqSystem
                    case 'ISS' 
                    if floor(sMtg.v_HolesMtg(p)/1000) > 0 & floor(sMtg.v_HolesMtg(p)/1000) <= 32%a1b2...                    
                        %Fond noir pour voir le blanc              
                        vColor = [80/255,240/255,60/255];
                    elseif floor(sMtg.v_HolesMtg(p)/1000) > 32 &  floor(sMtg.v_HolesMtg(p)/1000) <= 64 %cd bank
                        vColor = [255/255 128/255 0];
                    elseif floor(sMtg.v_HolesMtg(p)/1000) > 64 &  floor(sMtg.v_HolesMtg(p)/1000) <= 96 %ef bank
                        vColor = [169/255 0 240/255];
                    elseif floor(sMtg.v_HolesMtg(p)/1000) > 96 %gh
                        vColor = [255/255 0/255 128/255];
    %                     vColor = get_ItemColor( oDispOpt, 'Opt_LabelSrcHole' );
                    else                   
                        vColor = [0/255 255/255 255/255];%get_ItemColor( oDispOpt, 'Opt_LabelSrcHole' );
                    end
                    case  'IMAGINC' %or NIRx nomenclature
                        vColor = [255/255, 40/255,40/255];
                     case  'NIRx' %or NIRx nomenclature
                        vColor = [255/255, 40/255,40/255];
                end
            end
        
        if( get_DispOptChecked( oDispOpt, 'Opt_FittedOnSkin' ) )
            OptVertexBuffer = [OptVertexBuffer(:,1),OptVertexBuffer(:,2),OptVertexBuffer(:,3).*vHoles(p).SkinDepth,OptVertexBuffer(:,4)];%vHoles(p).SkinDepth;
        end
        
        OptVertexBuffer = OptVertexBuffer*vHoles(p).Transformation;
        
        v_hPrevDispItems.Optode(v_iToDisplay(Pos)) = patch( 'Vertices', [OptVertexBuffer(:,1), OptVertexBuffer(:,2), OptVertexBuffer(:,3)],...
                                                            'Faces', OptIndexBuffer, ...
                                                            'FaceColor', vColor, ...
                                                            'EdgeColor', 'none', ...
                                                            'FaceLighting', 'gouraud' );
    end
    fontsizenb = 14;
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    %Affichage des marqueurs d'intersection de cortex
    if( get_DispOptChecked( oDispOpt, 'Opt_CortexIntersectionMarker' ) )
        
        [MrkVertexBuffer, MrkIndexBuffer] = CreateSphereMesh( 0.0054, 5 );
        
        for( Pos=1:length(v_iToDisplay) )
        
            %Si la profondeur du cortex est non-nulle (autrement ca n'a pas
            %de sens)
            if( vHoles(p).CortexDepth )
                p = v_pActuDispItems(v_iToDisplay(Pos));

            
            %Detecteur
            if( find( sMtg.v_pDet == p ) )
                     vColor = [255/255, 40/255,40/255];
                     %vColor = [0, 0, 0];
            %Source
            else 
                if floor(sMtg.v_HolesMtg(p)/1000) <= 32 %a1b2... 
                    vColor = [255/255 255/255 0];
                elseif floor(sMtg.v_HolesMtg(p)/1000) > 32 &  floor(sMtg.v_HolesMtg(p)/1000) <= 64
                    vColor = [169/255 0 240/255];
                elseif floor(sMtg.v_HolesMtg(p)/1000) > 64 &  floor(sMtg.v_HolesMtg(p)/1000) <= 96
                    vColor = [169/255 0 240/255];
                elseif floor(sMtg.v_HolesMtg(p)/1000) > 96
                    vColor = [255/255 0/255 128/255];
                else
                    vColor = [0/255 255/255 255/255];%get_ItemColor( oDispOpt, 'Opt_LabelSrcHole' );
                end
            
            end

                v_hPrevDispItems.CortexIntersection(v_iToDisplay(Pos)) = patch( 'Vertices',  [MrkVertexBuffer(:,1)+(vHoles(p).Coord.x-vHoles(p).Normal.x*vHoles(p).CortexDepth),...%*vHoles(p).CortexDepth),...
                                                                                              MrkVertexBuffer(:,2)+(vHoles(p).Coord.y-vHoles(p).Normal.y*vHoles(p).CortexDepth),...%*vHoles(p).CortexDepth),....CortexDepth), ...
                                                                                              MrkVertexBuffer(:,3)+(vHoles(p).Coord.z-vHoles(p).Normal.z*vHoles(p).CortexDepth)], ....CortexDepth)], ...%vHoles(p).CortexDepth)], ....CortexDepth)], ...
                                                                                'Faces', MrkIndexBuffer, ...
                                                                                'FaceColor', vColor, ...
                                                                                'EdgeColor', 'none', ...
                                                                                'FaceLighting', 'gouraud' );
            end
        end
    end

    FontWeight = 'normal';
    if( get_DispOptChecked( oDispOpt, 'App_BoldLabels' ) )
        FontWeight = 'bold';
    end
    
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    %Affichage des labels de trous d'optodes
    if( get_DispOptChecked( oDispOpt, 'Lbl_OptodesHoleId' ) )
        for( Pos=1:length(v_iToDisplay) )
            
            p = v_pActuDispItems(v_iToDisplay(Pos));
            
           if get_DispOptChecked( oDispOpt, 'App_InvertedColorScheme' )
               vColor = [1, 1, 1]; % Change la couleur du texte pour le blanc.
           else
               vColor = [0, 0, 0];         
           end
            v_hPrevDispItems.LabelsHoleId(v_iToDisplay(Pos)) = text(  vHoles(p).Coord.x+vHoles(p).Normal.x*0.0025, ...
                                                                      vHoles(p).Coord.y+vHoles(p).Normal.y*0.0025, ...
                                                                      vHoles(p).Coord.z+vHoles(p).Normal.z*0.0025, ...
                                                                      vHoles(p).Label, 'FontSize', fontsizenb, ...
                                                                      'FontWeight', FontWeight, ...
                                                                      'Color', vColor, ...
                                                                      'VerticalAlignment', 'bottom');
        end
    end

    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------------------------------------------------------------------
    %Affichage des labels de fibres d'optodes
    %Les labels de fibres d'optodes doivent TOUS etre reaffiches car
    %pour une meme position qu'un affichage precedenet, le label de fibre
    %peut changer. Par consequent, il faut boucler pour tous les elements
    %de "v_pActuDispItems", et non seulement ceux qui ont change de
    %position.
    if( get_DispOptChecked( oDispOpt, 'Lbl_OptodesFibersId' ) )
        for( iPos=1:length(v_pActuDispItems) )
        
            %La visibilite doit etre verifiee car contrairement au boucles
            %precedentes, "v_iToDisplay" n'est pas utilise.
            if( v_bVisible(iPos) )
                
                p = v_pActuDispItems(iPos);            
                if get_DispOptChecked( oDispOpt, 'App_InvertedColorScheme' ) %Black background white string
                    vColor = [1, 1, 1]; 
                    set(gca,'color','black') % Crée un 'cube' de fond noir autour du model 3D.
                else            
                    vColor = [0, 0, 0];
                    set(gca,'color','white') % Crée un 'cube' de fond blanc autour du model 3D. Permet de voir les labels (noir) quand le background est noir lui aussi.
                end    
                strFiber = get_HoleFiberID( oHelmet, p );      
                if numel(strFiber)==4
                    strFiber = strFiber(1:2);
                elseif numel(strFiber)==6
                   strFiber = strFiber(1:3);   
                elseif numel(strFiber)==5;
                    if str2num(strFiber(2))==9
                        strFiber = strFiber(1:2);
                    else
                        strFiber = strFiber(1:3);
                    end                    
                end

                v_hPrevDispItems.LabelsFibersId(iPos) = text(  vHoles(p).Coord.x+vHoles(p).Normal.x*0.0025, ...
                                                               vHoles(p).Coord.y+vHoles(p).Normal.y*0.0025, ...
                                                               vHoles(p).Coord.z+vHoles(p).Normal.z*0.0025, ...
                                                               strFiber, 'FontSize', fontsizenb, 'Color',  ...
                                                               vColor, ...
                                                               'FontWeight', FontWeight, ...
                                                               'VerticalAlignment', 'top');
            end
        end
    end