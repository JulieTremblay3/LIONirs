% FitHelmetOnSubject
%
% Intrants: 
%    oDig_SubjectFiducials ->Digitalisation faite sur le sujet. 
%                            2 possibilités: 
%                                1)Digitalisation contenant les
%                                  positions directement digitalisés sur la
%                                  tête du sujet. Dans ce cas, le recalage
%                                  du casque complet sur cette digitalisation
%                                  n'est pas absolument nécéssaire.
%                                2)Digitalisation ne contenant que les
%                                  fiducies et les références X1,X2,X3.
%                                  Dans ce cas, le casque complet devra
%                                  être recalé sur cette digitalisation.
%
%    oDig_CompleteHelmet   ->Digitalisation complete du casque faite sur
%                            le support de calibration. Cette digitalisation
%                            ne comporte pas de points de fiducies valides,
%                            donc elle doit etre recalée dans le système de
%                            coordonnées du sujet (oDig_SubjectFiducials)
%                            afin d'avoir une quelconque validité.
% Extrants: 
%    oHelmetOnSubject      ->Casque dont les digitalisations sont recalés
%                            dans le système de coordonnées du sujet.
%
% Commentaires: Fonction permettant de recaler la digitalisation du casque 
%               complet dans le systeme de coordonnees du sujet. Les points
%               de reference X1 X2 X3 sont utilises pour y parvenir. Voir 
%               fonction "RegisterDigitization" pour plus de details.
function oHelmetOnSubject = FitHelmetOnSubject( oDig_SubjectFiducials, oDig_CompleteHelmet )
    
    %Détermination de la source et de la destination 
    %(en terme de syst. de coord.)
    oDig_Source = oDig_CompleteHelmet;          %Source = Complete Helmet
    oDig_Destination = oDig_SubjectFiducials;   %Destination = Partial Helmet (Subject Fiducials)
    [oDig_RegisteredDest, bSrcX1X2X3Found, bDstX1X2X3Found] = RegisterDigitization( oDig_Source, oDig_Destination );
    
    %Gestion des erreurs: si la registration nèa pas fonctionnée
    %Source X1 X2 X3 ref points not found, or nothing to transfer from Source.
    if( ~bSrcX1X2X3Found || isempty(oDig_Source) )
        if( ~bSrcX1X2X3Found )
            disp( 'FitHelmetOnSubject::Incomplete Complete Helmet Digitization (X1 X2 X3 missing)' );
        end
        if( isempty(oDig_Source) )
            disp( 'FitHelmetOnSubject::Empty Complete Helmet Digitization (vProbes empty)' );
        end
        disp( '    ->Partial Helmet only will be Loaded' );
        oHelmetOnSubject = set_vHolesAsFromPartialHelmet( Helmet(oDig_SubjectFiducials), 1:numel( get_vProbes(oDig_SubjectFiducials) ) );
        
        return;
    end
    
    %Destination X1 X2 X3 ref points not found, or nothing to transfer from Source.
    if( ~bDstX1X2X3Found || isempty(oDig_Destination) )
        disp( 'FitHelmetOnSubject::Empty or incomplete Partial Helmet Digitization (X1 X2 X3 missing?)' );
        bok = warndlg_modal( 'No valid suject digitization found. Consequently, this project should only be used to visualize the complete helmet digitization.' );
        
        %Note: Helmet is not really on subject because it as not been
        %sucessufully registered. However, it is still practical to allow 
        %it's visualisation.
        oHelmetOnSubject = set_vHolesAsFromCompleteHelmet( Helmet(oDig_CompleteHelmet), 1:numel( get_vProbes(oDig_CompleteHelmet) ) );
        return;
    end
        
    %Objet Helmet créé à partir des coordonnées digitalisés sur le sujet
    %La fonction "set_vHolesAsFromPartialHelmet" crée un Helmet dont chaque
    %trou est identifié comme faisant partie de la digitalisation
    %partielle.
    %Les normales ici sont déja calculées par le constructeur paramétrique
    %de Helmet.
    oSubjectHelm = set_vHolesAsFromPartialHelmet( Helmet(oDig_SubjectFiducials), 1:numel( get_vProbes(oDig_SubjectFiducials) ) );
    
    %Objet Helmet créé à partir des coordonnées digitalisés sur le support
    %de calibration (coordonnées du casque déja recalées sur le sujet).
    %La fonction "set_vHolesAsFromPartialHelmet" crée un Helmet dont chaque
    %trou est identifié comme faisant partie de la digitalisation
    %complete du casque.
    %Les normales ici sont déja calculées par le constructeur paramétrique
    %de Helmet.
    oRegHelm = set_vHolesAsFromCompleteHelmet( Helmet(oDig_RegisteredDest), 1:numel( get_vProbes(oDig_RegisteredDest) ) );
    
    %Casque partiel:
    %Vecteur de trous de casque digitalisés (chaque elem contient son Label,
    %sa Normale, sa Position, son type, et les autres flags).
    vHolesSub = get_vHoles( oSubjectHelm );
    
    %Casque complet:
    %Vecteur de trous de casque digitalisés (chaque elem contient son Label,
    %sa Normale, sa Position, son type, et les autres flags).
    vHolesReg = get_vHoles( oRegHelm );
    
    %Trous de casque devant être transférés (tous les trous normaux qui ne
    %sont pas des trous de référence ou de test)
    vHolesToBeTransfered = vHolesReg( find([vHolesReg.Type] == get_PtTypeNo(oRegHelm, 'NormalHole')) );
    
    %Boucle de copie des trous transformes du casque complet
    for( iElem = 1:numel(vHolesToBeTransfered) )
        p=numel(vHolesSub)+1;
        %Verifier que ce nèest pas une trou X1, X2, X3, ou X4
        if(    ~strcmp(vHolesToBeTransfered(iElem).Label, 'X1') && ~strcmp(vHolesToBeTransfered(iElem).Label, 'X2') ...
            && ~strcmp(vHolesToBeTransfered(iElem).Label, 'X3') && ~strcmp(vHolesToBeTransfered(iElem).Label, 'X4') )
            vHolesSub(p) = vHolesToBeTransfered(iElem);
        end
    end
    
    oHelmetOnSubject = set_vHoles( oSubjectHelm, vHolesSub );

    %Calcul de la transformation de chaque point (matrice utilisee pour
    %l'affichage)
    oHelmetOnSubject = ComputeNeighborhood( ComputeTransform( oHelmetOnSubject ) );
    
    
    
    
%     
%     
%     
%     
%     
%     vHolesS = get_vHoles( oSubject );
%     vHolesH = get_vHoles( oHelmet );
%     
%     %Transfert des coordonnees sous forme de matrices
%     CoordsH = [vHolesH(find([vHolesH.Type] == 400)).Coord];
%     NormalsH = [vHolesH(find([vHolesH.Type] == 400)).Normal];
%     
%     if( isempty(CoordsH) )
%         disp( 'Coregistration points not found in Complete Digitization. Subject Digitization only will be used.' );
%         oHelmetOnSubject = oSubject;
%         return;
%     end
% 
%     matCoordsH = [CoordsH.x;CoordsH.y;CoordsH.z]';
%     matNormalsH = [NormalsH.x;NormalsH.y;NormalsH.z]';
% 
%     %Recherche des references du casque complet
%     iX1H = find(strcmp({vHolesH.Label}, 'X1'));
%     iX2H = find(strcmp({vHolesH.Label}, 'X2'));
%     iX3H = find(strcmp({vHolesH.Label}, 'X3'));
%     
%     if( isempty(iX1H) ||isempty(iX2H) || isempty(iX3H) )
%         disp( 'Coregistration points not found in Complete Digitization. Subject Digitization only will be used.' );
%         oHelmetOnSubject = oSubject;
%         return;
%     end
% 
%     %Recherche des references du casque partiel du sujet
%     iX1S = find(strcmp({vHolesS.Label}, 'X1'));
%     iX2S = find(strcmp({vHolesS.Label}, 'X2'));
%     iX3S = find(strcmp({vHolesS.Label}, 'X3'));
%     
%     if( isempty(iX1S) ||isempty(iX2S) || isempty(iX3S) )
%         if( numel(vHolesS) == 1 ) 
%             disp( 'Helmet Reference points not found in Subject Digitization. Complete Helmet Digitization only will be used.' );
%             oHelmetOnSubject = oHelmet;
%             return;
%         else
%             warndlg_modal( 'Helmet Reference points (X1 X2 X3 X4) not found in Subject Digitization. Complete Helmet cannot be coregistered so it won''t be used.' );
%             oHelmetOnSubject = oSubject;
%             return;
%         end
%     end
% 
%     matHelmetRefs = [ vHolesH(iX1H).Coord.x, vHolesH(iX1H).Coord.y, vHolesH(iX1H).Coord.z ; ...
%                       vHolesH(iX2H).Coord.x, vHolesH(iX2H).Coord.y, vHolesH(iX2H).Coord.z ; ...
%                       vHolesH(iX3H).Coord.x, vHolesH(iX3H).Coord.y, vHolesH(iX3H).Coord.z ];
% 
%     matSubjectRefs = [ vHolesS(iX1S).Coord.x, vHolesS(iX1S).Coord.y, vHolesS(iX1S).Coord.z ; ...
%                        vHolesS(iX2S).Coord.x, vHolesS(iX2S).Coord.y, vHolesS(iX2S).Coord.z ; ...
%                        vHolesS(iX3S).Coord.x, vHolesS(iX3S).Coord.y, vHolesS(iX3S).Coord.z ];
% 
% 
%     matCoordsH(:,4) = 1;
%     matNormalsH(:,4) = 1;
% 
%     %Transformation a 6 Degres de liberte (Rigid body transform)
% 	%(6DOF = rotation xyz et translation xyz)
%     %Application de la matrice de coregistration pour les positions
% 	matCoordsTransformed = matCoordsH*Create_6DOF_CoRegistration_Matrix( matHelmetRefs, matSubjectRefs, [1,1,1], [1,1,1] );
%     
%     %Transformation a 3 Degres de liberte (Affine transform)
% 	%(3DOF = rotation xyz) <-Bref, pas de translation pour les normales
%     %Application de la matrice de coregistration pour les normales
%     matNormalsTransformed = matNormalsH*Create_3DOF_CoRegistration_Matrix( matHelmetRefs, matSubjectRefs );
%     
%     HolesToBeTransfered = vHolesH(find([vHolesH.Type] == 400));
%     
%     %Boucle de copie des trous transformes du casque complet
%     for( iElem = 1:size(matCoordsTransformed,1) )
%         p=numel(vHolesS)+1;
%         if(    ~strcmp(HolesToBeTransfered(iElem).Label, 'X1') && ~strcmp(HolesToBeTransfered(iElem).Label, 'X2') ...
%             && ~strcmp(HolesToBeTransfered(iElem).Label, 'X3') && ~strcmp(HolesToBeTransfered(iElem).Label, 'X4') )
%             vHolesS(p) = HolesToBeTransfered(iElem);
%             vHolesS(p).Coord.x = matCoordsTransformed(iElem,1);
%             vHolesS(p).Coord.y = matCoordsTransformed(iElem,2);
%             vHolesS(p).Coord.z = matCoordsTransformed(iElem,3);
%             vHolesS(p).Normal.x = matNormalsTransformed(iElem,1);
%             vHolesS(p).Normal.y = matNormalsTransformed(iElem,2);
%             vHolesS(p).Normal.z = matNormalsTransformed(iElem,3);
%         end
%     end
%     
%     oHelmetOnSubject = set_vHoles( oSubject, vHolesS );
% 
%     %Calcul de la transformation de chaque point (matrice utilisee pour
%     %l'affichage)
%     oHelmetOnSubject = ComputeNeighborhood( ComputeTransform( oHelmetOnSubject ) );
% %     
% %     vProbes_Sub = get_vProbes( oDig_SubjectFiducials );
% %     vProbes_Helm = get_vProbes( oDig_CompleteHelmet );
% %     
% %     matCoordsHelm = reshape([vProbes_Helm.Coord]', 3, numel([vProbes_Helm.Coord])/3)';
% %     matCoordsHelm = matCoordsHelm( find( [vProbes_Helm.Type] == get_PtTypeNo(Helmet(), 'NormalHole' ) ), :);
% %     
% %     oHelmetOnSubject = Helmet;
% % 
% %     iX1H = find(strcmp({vProbes_Helm.Label}, 'X1'));
% %     iX2H = find(strcmp({vProbes_Helm.Label}, 'X2'));
% %     iX3H = find(strcmp({vProbes_Helm.Label}, 'X3'));
% % 
% %     %Recherche des references du casque partiel du sujet
% %     iX1S = find(strcmp({vProbes_Sub.Label}, 'X1'));
% %     iX2S = find(strcmp({vProbes_Sub.Label}, 'X2'));
% %     iX3S = find(strcmp({vProbes_Sub.Label}, 'X3'));
% % 
% %     matHelmetRefs = [ vProbes_Helm(iX1H).Coord ; ...
% %                       vProbes_Helm(iX2H).Coord ; ...
% %                       vProbes_Helm(iX3H).Coord ];
% % 
% %     matSubjectRefs = [ vProbes_Sub(iX1S).Coord ; ...
% %                        vProbes_Sub(iX2S).Coord ; ...
% %                        vProbes_Sub(iX3S).Coord ];
% % 
% % 	matCoRegistration = Create_CoRegistration_Matrix( matHelmetRefs, matSubjectRefs, [1,1,1], [1,1,1] );
% %     
% %     matCoordsHelm(:,4) = 1;
%     
