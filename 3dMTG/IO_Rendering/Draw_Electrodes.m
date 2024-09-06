function Draw_Electrodes( oDispOpt, oInterDlgComm, oHelmet )

    persistent v_hPrevDispItems; %Vecteur d'handles de chaque elements pour chaque position
    persistent v_pPrevDispItems; %Vecteur de position (position=indice lineaire de trou dans vHoles)

    sMtg          = get_Mtg( oHelmet );
    vHoles        = get_vHoles( oHelmet );
    sHelmetAxeDisp = get_HelmetAxeDisp( oInterDlgComm );
    v_pActuDispItems = sMtg.v_pEle;
    fontsizenb = 14;
    %Uniformisation des tailles afin de permettre les algorithmes
    %vectoriels. Le padding a zero n'influence pas les resultats.
    NbActuItems = length(v_pActuDispItems); %Nombre de positions actuellement visibles
    NbPrevItems = length(v_pPrevDispItems); %Nombre de positions precedemment visibles (lors du dernier appel la fonction)
    MaxLength = max(NbActuItems, NbPrevItems);
    v_pActuDispItems(NbActuItems+1:MaxLength) = 0;  %Padding a zero
    v_pPrevDispItems(NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    v_hPrevDispItems.holes(NbPrevItems+1:MaxLength)= 0;
    v_hPrevDispItems.label(NbPrevItems+1:MaxLength) = 0;  %Padding a zero
    v_hPrevDispItems.LabelsEleId(NbPrevItems+1:MaxLength) = 0;
    %Sortir s'il n'y a rien a faire (rien a effacer, rien a afficher)
	if( ~MaxLength )
        return;
    end
    
    %Validation des handles: 
    v_hPrevDispItems.label     = v_hPrevDispItems.label        .* ishandle(v_hPrevDispItems.label);
    v_hPrevDispItems.holes     = v_hPrevDispItems.holes        .* ishandle(v_hPrevDispItems.holes);
    v_hPrevDispItems.LabelsEleId = v_hPrevDispItems.LabelsEleId .*ishandle(v_hPrevDispItems.LabelsEleId);
    %Les handles d'items affiches precedemment doivent etre non nuls,
    %autrement on ne peut considerer que les items sont encore dans le
    %graphique. Ils ont peut-etre ete effaces entre 2 appels de la fonction
    %d'affichage, avec la commande 'cla' par exemple. Ceci (ci-bas) permet donc de
    %verifier afin de reafficher les items effaces.
    v_pPrevDispItems = v_pPrevDispItems .*((v_hPrevDispItems.label ~= 0)...
                                            |(v_hPrevDispItems.holes~=0)...
                                            |(v_hPrevDispItems.LabelsEleId ~=0));


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
    
    delete( v_hPrevDispItems.holes( v_iToErase( find(v_hPrevDispItems.holes( v_iToErase) ) ) ) ); %Delete sur elems non-nuls seulement (le 'find' sert a retirer les elements nuls)
    delete( v_hPrevDispItems.label( v_iToErase( find(v_hPrevDispItems.label( v_iToErase) ) ) ) );
    delete( v_hPrevDispItems.LabelsEleId( v_iToErase( find(v_hPrevDispItems.LabelsEleId( v_iToErase) ) ) ) );
    v_hPrevDispItems.holes( v_iToErase ) = 0;
    v_hPrevDispItems.label( v_iToErase ) = 0;
    v_hPrevDispItems.LabelsEleId( v_iToErase ) = 0;
    v_pPrevDispItems( v_iToErase ) = 0;
    
    if ( get_DispOptChecked( oDispOpt, 'Ele_PositionOnHelmet' ) )
        delete( v_hPrevDispItems.holes( v_iToErase( find(v_hPrevDispItems.holes( v_iToErase) ) ) ) ); %Delete sur elems non-nuls seulement (le 'find' sert a retirer les elements nuls)
        v_hPrevDispItems.holes( v_iToErase ) = 0;
        delete( v_hPrevDispItems.label( v_iToErase( find(v_hPrevDispItems.label( v_iToErase) ) ) ) ); %Delete sur elems non-nuls seulement (le 'find' sert a retirer les elements nuls)
        v_hPrevDispItems.label( v_iToErase ) = 0;
        delete(v_hPrevDispItems.LabelsEleId( v_iToErase( find(v_hPrevDispItems.LabelsEleId( v_iToErase) ) ) ) ); 
        v_hPrevDispItems.LabelsEleId( v_iToErase ) = 0;
    else
        %Delete de tous les Handles valides (les Handles valides sont des
        %elements encore affiches) si cet option d'affichage est desactive
        delete( v_hPrevDispItems.holes( find( v_hPrevDispItems.holes ) ) );
        v_hPrevDispItems.holes(         find( v_hPrevDispItems.holes ) ) = 0;
        delete( v_hPrevDispItems.label( find(  v_hPrevDispItems.label) ) );
        v_hPrevDispItems.label(         find(  v_hPrevDispItems.label ) ) = 0;
        delete( v_hPrevDispItems.LabelsEleId( find(  v_hPrevDispItems.LabelsEleId) ) );
        v_hPrevDispItems.LabelsEleId(         find(  v_hPrevDispItems.LabelsEleId ) ) = 0;
        %Rien a afficher. Sortir.
        return;
    end
    %Les labels de fibres doivent tous etre redessines car pour une meme 
    %position les fibres peuvent changer entre 2 appels d'affichage. Donc,
    %on les efface toutes.
    delete( v_hPrevDispItems.label( find( v_hPrevDispItems.label ) ) );
    delete( v_hPrevDispItems.LabelsEleId( find(v_hPrevDispItems.LabelsEleId ) ) );
    
    v_hPrevDispItems.label( find( v_hPrevDispItems.label) ) = 0;
    
    
    %Determiner les positions a afficher: 
    v_iToDisplay = find ( v_bVisible & ( v_pPrevDispItems ~= v_pActuDispItems ) );
    v_pPrevDispItems( v_iToDisplay ) = v_pActuDispItems( v_iToDisplay );

    %[EleSrc_VertexBuffer, EleIndexBuffer] = CreateCircleMesh( 0.0054, 12 );
   % [EleDet_VertexBuffer, EleIndexBuffer] = CreateCircleMesh( 0.0110, 12 );


   [EleSrc_VertexBuffer, EleIndexBuffer] = CreateSphereMesh( 0.0110, 16 );
   [EleDet_VertexBuffer, EleIndexBuffer] = CreateSphereMesh( 0.0110, 16 );
    
   
	%vColor = get_ItemColor( oDispOpt, 'Ele_PositionOnHelmet' );
    vColor = [0/255 0/255 255/255];
     vColor = [255/255 255/255 255/255];
    %Redessine les ronds
    for( Pos=1:length(v_iToDisplay) )        
        p = v_pActuDispItems(v_iToDisplay(Pos));       
        if( vHoles(p).CanBeDet ) %Trou de detecteur (grand diametre)
            EleVertexBuffer = EleDet_VertexBuffer;
        else %Trou de source (petit diametre)
            EleVertexBuffer = EleSrc_VertexBuffer;
        end        
        EleVertexBuffer = EleVertexBuffer*vHoles(p).Transformation;        
        v_hPrevDispItems.holes(v_iToDisplay(Pos)) = patch( 'Vertices', [EleVertexBuffer(:,1), EleVertexBuffer(:,2), EleVertexBuffer(:,3)],...
                                                     'Faces', EleIndexBuffer, ...
                                                     'FaceColor', vColor, ...
                                                     'EdgeColor', 'none', ...
                                                     'FaceLighting', 'gouraud' );
    end
    %Redessine les noms électrodes
    for( Pos=1:length(v_pActuDispItems) )
        p = v_pActuDispItems(Pos);
        if p ~=0
        if isempty(sMtg.v_HolesEle{p})
           label_ele = ' ';           
        else
           label_ele = sMtg.v_HolesEle{p};
        end
        if get_DispOptChecked( oDispOpt, 'App_InvertedColorScheme' ) % Change la couleur du texte pour le blanc...
        v_hPrevDispItems.label(Pos)=  text(  vHoles(p).Coord.x+vHoles(p).Normal.x*0.0025, ...     
                                             vHoles(p).Coord.y+vHoles(p).Normal.y*0.0025, ...
                                             vHoles(p).Coord.z+vHoles(p).Normal.z*0.0025, ...
                                             label_ele , 'FontSize',fontsizenb ,'VerticalAlignment', 'bottom','color', 'w');      
        else % ou le noir.
                 v_hPrevDispItems.label(Pos)=  text(  vHoles(p).Coord.x+vHoles(p).Normal.x*0.0025, ...     
                                             vHoles(p).Coord.y+vHoles(p).Normal.y*0.0025, ...
                                             vHoles(p).Coord.z+vHoles(p).Normal.z*0.0025, ...
                                             label_ele , 'FontSize',fontsizenb ,'VerticalAlignment', 'bottom','color', 'k'); 
        end
        end
    end
      if get_DispOptChecked( oDispOpt, 'App_InvertedColorScheme' ) % Verifier si les couleurs sont inversees et adapter la couleur de la zone entourant le model 3D.
        set(gca,'color','black')
      else
        set(gca,'color','white')
      end
    % Redessine les noms trous
%     for ( Pos=1:length(v_pActuDispItems) )
%         p = v_pActuDispItems(Pos);
%          if p ~=0
%           v_hPrevDispItems.LabelsEleId(Pos) = text(  vHoles(p).Coord.x+vHoles(p).Normal.x*0.0025, ...
%                                                                       vHoles(p).Coord.y+vHoles(p).Normal.y*0.0025, ...
%                                                                       vHoles(p).Coord.z+vHoles(p).Normal.z*0.0025, ...
%                                                                       vHoles(p).Label, 'FontSize',fontsizenb , 'VerticalAlignment', 'top'  );
%                                                                        
%         end
%     end
   
