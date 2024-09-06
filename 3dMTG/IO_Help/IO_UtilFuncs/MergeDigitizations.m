% MergeDigitizations
%
% Fonction permettant de combiner 2 digitalisations en les intégrant dans
% un même système de cordonnées. Afin d'effectuer le passage de systèeme de
% coordonnées, les points X1 X2 X3 sont utilisés comme référence afin
% d'effectuer une "CoRegistration". Ainsi les points X1 X2 X3 des deux
% digitalisations se superposent.
%
% À noter: Les coordonnées des positions de la digitalisation "Destination" 
% resteront intactes dans le résultat "Merged". Ce sont les coordonnées de
% la digitalisation "Source" qui seront translatés et rotationés afin
% de matcher ses points X1 X2 X3 sur les  X1 X2 X3 de la "destination".
% 
% À noter: Si une position de la source a le meme label qu'une position de
% destination, alors la position de destination verra son "Label" se faire
% ajouter un préfixe "Old_".
function [oDig_Merged, bSrcX1X2X3Found, bDstX1X2X3Found] = MergeDigitizations( oDig_Source, oDig_Destination )
    
    %Objet Digitization vide
    oDig_Merged = Digitization();
    bSrcX1X2X3Found = false;
    bDstX1X2X3Found = false;
    
    %Construction des deux objets "Helmet" (un pour la source, un pour
    %la destination). Les objets "Helmet" donnent acces a plus de fonctions
    %que les objets "Digitization" et c'est pourquoi on les construit
    %temporairement dans cette fonction.
    oDstHelm = set_vHolesAsFromPartialHelmet(  Helmet(oDig_Destination), 1:numel( get_vProbes(oDig_Destination) ) );
    oSrcHelm = set_vHolesAsFromCompleteHelmet( Helmet(oDig_Source),      1:numel( get_vProbes(oDig_Source)      ) );
    
    vDstHoles = get_vHoles( oDst );
    vSrcHoles = get_vHoles( oHelmet );
    
    %Transfert des coordonnees sous forme de vecteurs de struct.x/.y/.z. 
    %Seuls les points du type NormalHole seront transferés.
    v_sSrcCoords =  [vHolesSrc(find([vHolesH.Type] == get_PtTypeNo('NormalHole'))).Coord ];
    v_sSrcNormals = [vHolesSrc(find([vHolesH.Type] == get_PtTypeNo('NormalHole'))).Normal];
    
    %Recherche des references X1,X2,X3 dans la source
    iSrcX1 = find(strcmp({vSrcHoles.Label}, 'X1'));
    iSrcX2 = find(strcmp({vSrcHoles.Label}, 'X2'));
    iSrcX3 = find(strcmp({vSrcHoles.Label}, 'X3'));
    
    %Recherche des references X1,X2,X3 dans la destination
    iDstX1 = find(strcmp({vDstHoles.Label}, 'X1'));
    iDstX2 = find(strcmp({vDstHoles.Label}, 'X2'));
    iDstX3 = find(strcmp({vDstHoles.Label}, 'X3'));
    
    bSrcX1X2X3Found = ( ~isempty(iSrcX1) && ~isempty(iSrcX2) && ~isempty(iSrcX3) );
    bDstX1X2X3Found = ( ~isempty(iDstX1) && ~isempty(iDstX2) && ~isempty(iDstX3) );
    
    %Source X1 X2 X3 ref points not found, or nothing to transfer.
    if( ~bSrcX1X2X3Found || isempty(v_sSrcCoords) )
        disp( 'Empty or incorrect Source Digization (X1 X2 X3 missing?)' );
        oDig_Merged = oDig_Destination;
        return;
    end
    
    %Destination X1 X2 X3 ref points not found.
    if( ~bDstX1X2X3Found )
        %if( numel(vDstHoles) == 1 ) 
        %    disp( 'Destination Digitization Empty. Using Source Digitization' );
        %    oDig_Merged = oDig_Source;
        %    return;
        %else
            disp( 'Destination Digization reference points missing(X1 X2 X3)' );
            oDig_Merged = oDig_Destination;
            return;
        %end
    end
    
    matSrcCoords  = [v_sSrcCoords.x ;v_sSrcCoords.y ;v_sSrcCoords.z ]';
    matSrcNormals = [v_sSrcNormals.x;v_sSrcNormals.y;v_sSrcNormals.z]';

    matSrcRefs = [ vSrcHoles(iSrcX1).Coord.x, vSrcHoles(iSrcX1).Coord.y, vSrcHoles(iSrcX1).Coord.z ; ...
                   vSrcHoles(iSrcX2).Coord.x, vSrcHoles(iSrcX2).Coord.y, vSrcHoles(iSrcX2).Coord.z ; ...
                   vSrcHoles(iSrcX3).Coord.x, vSrcHoles(iSrcX3).Coord.y, vSrcHoles(iSrcX3).Coord.z ];

    matDstRefs = [ vDstHoles(iDstX1).Coord.x, vDstHoles(iDstX1).Coord.y, vDstHoles(iDstX1).Coord.z ; ...
                   vDstHoles(iDstX2).Coord.x, vDstHoles(iDstX2).Coord.y, vDstHoles(iDstX2).Coord.z ; ...
                   vDstHoles(iDstX3).Coord.x, vDstHoles(iDstX3).Coord.y, vDstHoles(iDstX3).Coord.z ];

    matCoordsH(:,4) = 1;
    matNormalsH(:,4) = 1;

    %Transformation a 6 Degres de liberte (Rigid body transform)
	%(6DOF = rotation xyz et translation xyz)
    %Application de la matrice de coregistration pour les positions
	matSrcCoordsTransformed = matSrcCoords*Create_6DOF_CoRegistration_Matrix( matSrcRefs, matDstRefs, [1,1,1], [1,1,1] );
    
    %Transformation a 3 Degres de liberte (Affine transform)
	%(3DOF = rotation xyz) <-Bref, pas de translation pour les normales
    %Application de la matrice de coregistration pour les normales
    matSrcNormalsTransformed = matSrcNormals*Create_3DOF_CoRegistration_Matrix( matSrcRefs, matDstRefs );
    
    HolesToBeTransfered = vHolesH(find([vHolesH.Type] == 400));
    
    %Boucle de copie des trous transformes du casque complet
    for( iElem = 1:size(matCoordsTransformed,1) )
        p=numel(vHolesS)+1;
        if(    ~strcmp(HolesToBeTransfered(iElem).Label, 'X1') && ~strcmp(HolesToBeTransfered(iElem).Label, 'X2') ...
            && ~strcmp(HolesToBeTransfered(iElem).Label, 'X3') && ~strcmp(HolesToBeTransfered(iElem).Label, 'X4') )
            vHolesS(p) = HolesToBeTransfered(iElem);
            vHolesS(p).Coord.x = matCoordsTransformed(iElem,1);
            vHolesS(p).Coord.y = matCoordsTransformed(iElem,2);
            vHolesS(p).Coord.z = matCoordsTransformed(iElem,3);
            vHolesS(p).Normal.x = matNormalsTransformed(iElem,1);
            vHolesS(p).Normal.y = matNormalsTransformed(iElem,2);
            vHolesS(p).Normal.z = matNormalsTransformed(iElem,3);
        end
    end
    
    oHelmetOnSubject = set_vHoles( oSubject, vHolesS );

    %Calcul de la transformation de chaque point (matrice utilisee pour
    %l'affichage)
    oHelmetOnSubject = ComputeNeighborhood( ComputeTransform( oHelmetOnSubject ) );
%     
%     vProbes_Sub = get_vProbes( oDig_SubjectFiducials );
%     vProbes_Helm = get_vProbes( oDig_CompleteHelmet );
%     
%     matCoordsHelm = reshape([vProbes_Helm.Coord]', 3, numel([vProbes_Helm.Coord])/3)';
%     matCoordsHelm = matCoordsHelm( find( [vProbes_Helm.Type] == get_PtTypeNo(Helmet(), 'NormalHole' ) ), :);
%     
%     oHelmetOnSubject = Helmet;
% 
%     iX1H = find(strcmp({vProbes_Helm.Label}, 'X1'));
%     iX2H = find(strcmp({vProbes_Helm.Label}, 'X2'));
%     iX3H = find(strcmp({vProbes_Helm.Label}, 'X3'));
% 
%     %Recherche des references du casque partiel du sujet
%     iX1S = find(strcmp({vProbes_Sub.Label}, 'X1'));
%     iX2S = find(strcmp({vProbes_Sub.Label}, 'X2'));
%     iX3S = find(strcmp({vProbes_Sub.Label}, 'X3'));
% 
%     matHelmetRefs = [ vProbes_Helm(iX1H).Coord ; ...
%                       vProbes_Helm(iX2H).Coord ; ...
%                       vProbes_Helm(iX3H).Coord ];
% 
%     matSubjectRefs = [ vProbes_Sub(iX1S).Coord ; ...
%                        vProbes_Sub(iX2S).Coord ; ...
%                        vProbes_Sub(iX3S).Coord ];
% 
% 	matCoRegistration = Create_CoRegistration_Matrix( matHelmetRefs, matSubjectRefs, [1,1,1], [1,1,1] );
%     
%     matCoordsHelm(:,4) = 1;
    

