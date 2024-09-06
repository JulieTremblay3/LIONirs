% RegisterDigitization
%
% Fonction permettant de transférer une digitalisation "Source" dans le
% système de coordonnées d'une deuxième digitalisation dite de
% "Destination". Afin d'effectuer le passage de systèeme de
% coordonnées, les points X1 X2 X3 sont utilisés comme référence afin
% d'effectuer une "Registration". Ainsi, la source est translatée et 
% rotatée de manière à ce que ses points X1 X2 X3 superposent les points X1
% X2 X3 de la destination. La trasformation est affine ("Rigid-Body
% Transform") donc ne provoque pas de déformations. La relation entre les
% positions est conservée avec exactitude.
%
% Note: TOUS les éléments de la Digitalisation source sont transférés de
% système de coordonnees, y compris les fiducies et les X1 X2 X3.
function [oDigRegisteredSrc, bSrcX1X2X3Found, bDstX1X2X3Found] = RegisterDigitization( oDig_Source, oDig_Destination )
    
    %Objet Digitization vide
    oDigRegisteredSrc = Digitization();
    bSrcX1X2X3Found = false;
    bDstX1X2X3Found = false;
    
    vDstProbes = get_vProbes( oDig_Destination );
    vSrcProbes = get_vProbes( oDig_Source );
    
    %Transfert des coordonnees sous forme de matrice de vecteurs. 
    %matSrcCoords = vSrcProbes(:).Coord;
    matSrcCoords = reshape([vSrcProbes(:).Coord], 4, numel(vSrcProbes))';
    
    %Fiducies (Na;Lpa;Rpa)
    matSrcFids = get_matFiducials( oDig_Source );
    matSrcFids(:,4) = 1; %4è coord pour manipulation matricielle
    
    %Recherche des references X1,X2,X3 dans la source
    iSrcX1 = find(strcmp({vSrcProbes.Label}, 'X1'));
    iSrcX2 = find(strcmp({vSrcProbes.Label}, 'X2'));
    iSrcX3 = find(strcmp({vSrcProbes.Label}, 'X3'));
    
    %Recherche des references X1,X2,X3 dans la destination
    iDstX1 = find(strcmp({vDstProbes.Label}, 'X1'));
    iDstX2 = find(strcmp({vDstProbes.Label}, 'X2'));
    iDstX3 = find(strcmp({vDstProbes.Label}, 'X3'));
    
    bSrcX1X2X3Found = ( ~isempty(iSrcX1) && ~isempty(iSrcX2) && ~isempty(iSrcX3) );
    bDstX1X2X3Found = ( ~isempty(iDstX1) && ~isempty(iDstX2) && ~isempty(iDstX3) );
    
    %Source X1 X2 X3 ref points not found, or nothing to transfer.
    if( ~bSrcX1X2X3Found || isempty(matSrcCoords) )
        disp( 'Empty or incorrect Source Digitization (X1 X2 X3 missing?)' );
        return;
    end
    
    %Destination X1 X2 X3 ref points not found.
    if( ~bDstX1X2X3Found )
        disp( 'Destination Digitization reference points missing(X1 X2 X3)' );
        return;
    end
    
    matSrcRefs = [ vSrcProbes(iSrcX1).Coord ; ...
                   vSrcProbes(iSrcX2).Coord ; ...
                   vSrcProbes(iSrcX3).Coord ];

    matDstRefs = [ vDstProbes(iDstX1).Coord; ...
                   vDstProbes(iDstX2).Coord; ...
                   vDstProbes(iDstX3).Coord ];

    %Transformation a 6 Degres de liberte (Rigid body transform)
	%(6DOF = rotation xyz et translation xyz)
    matRegistration = Create_6DOF_CoRegistration_Matrix( matSrcRefs, matDstRefs, [1,1,1], [1,1,1] );
	
    %Application de la matrice de recalage sur les positions
    matSrcCoordsRegistered = matSrcCoords*matRegistration;
    
    %Application de la matrice de recalage sur les fiducies
    matFidsRegistered = matSrcFids*matRegistration;
    
    %Boucle de copie des trous transformes du casque complet
    for( iProbe = 1:numel(vSrcProbes) )
        %Transfert des information de la position numérisée (Type, Label, etc.)
        vRegisteredProbes(iProbe) = vSrcProbes(iProbe);
        %Transfert des coordonnées
        vRegisteredProbes(iProbe).Coord = matSrcCoordsRegistered(iProbe,:);
    end
        
    %Transfert dans l'objet de sortie
    oDigRegisteredSrc = set_vProbes( oDigRegisteredSrc, vRegisteredProbes );
    
    %Transfert dans l'objet de sortie (4è coord ignorée: n'est plus utile)
    oDigRegisteredSrc = set_matFiducials( oDigRegisteredSrc, matFidsRegistered(:,1:3) );
    
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
    

