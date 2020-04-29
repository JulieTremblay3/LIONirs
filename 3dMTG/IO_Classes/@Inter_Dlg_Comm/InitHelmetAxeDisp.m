function CurObj = InitHelmetAxeDisp( CurObj, hAxeHelmet, oHelmet )
    
    %Flag de camorbit en cours
    CurObj.HelmetAxeDisp.IsAxeRotating = false;
    
    sMtg = get_Mtg( oHelmet );

    % Handle du graphique de casque (axe handle)
    CurObj.HelmetAxeDisp.axe_handle = hAxeHelmet;

    % matrice contenant la position des curseurs (indices lineaires de matH)
    CurObj.HelmetAxeDisp.v_pSelection = [];     

    % Identificateur Actuel
    CurObj.HelmetAxeDisp.pIdentification = 0;
    CurObj.HelmetAxeDisp.pRulerFirst = 0;  %1er pt du calcul de distance
    CurObj.HelmetAxeDisp.pRulerLast = 0;   %2è  pt du calcul de distance

    %Matrice de visibilite des trous du casque : 
    %Afficher le casque en entier pour initialiser correctement la cible
    %de la caméra.
    CurObj.HelmetAxeDisp.v_bVisible          = ones(size(get_vHoles(oHelmet)));
    CurObj.HelmetAxeDisp.m_bRegionPtsVisible = ones(size(sMtg.mat_RegionBorder,1),size(sMtg.mat_RegionBorder,2));
    
    %Matrice de visibilite des trous du casque : 
    %Afficher le casque en entier pour initialiser correctement la cible
    %de la caméra.
    sMtg = get_Mtg(oHelmet);
    
    
       
