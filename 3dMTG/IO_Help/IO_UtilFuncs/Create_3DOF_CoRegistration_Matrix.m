function matCoReg = Create_3DOF_CoRegistration_Matrix( matSourceRefs, matDestinationRefs, bIsSrcSystemRH, bIsDestSystemRH )
%
% CRÉATION D'UNE MATRICE DE TRANSFORMATION À 3 DEGRÉS DE LIBERTÉ
% (rotations seulement, pas de translation comme 6DOF)
%
% CETTE TRANSFORMATION (3DOF) EST LINEAIRE, CONTRAIREMENT À 6DOF
% (ELLE PEUT PAR CONSEQUENT ETRE UTILISEE POUR TRANSFORMER DES NORMALES)
%
% INTRANTS:                  _                _
%           matSourceRefs : |  s1_x s1_y s1_z  |
%                           |  s2_x s2_y s2_z  |
%                           |_ s3_x s3_y s3_z _|
%                                 _                _
%           matDestinationRefs : |  d1_x d1_y d1_z  |
%                                |  d2_x d2_y d2_z  |
%                                |_ d3_x d3_y d3_z _|
%
%           bIsSrcSystemRH : Flag systeme source Right-Hand? (true/false)
%
%           bIsDestSystemRH : Flag systeme destination Right-Hand? (true/false)
%
% EXTRANTS:
%                                   
%           matCoReg : Matrice de transformation 4x4 avec les composantes
%                      de translation dans le bas.

    %*******************************************************************
    % ETAPE I : Mise a l'echelle selon la taille physique des voxels
    %*******************************************************************
    
    %Mise a l'echelle de la matrice referentielle source 
    matSourceRefs = [ matSourceRefs(:,1), ...
                      matSourceRefs(:,2), ...
                      matSourceRefs(:,3)  ];
                  
    %Mise a l'echelle de la matrice referentielle de destination 
    matDestinationRefs = [ matDestinationRefs(:,1), ...
                           matDestinationRefs(:,2), ...
                           matDestinationRefs(:,3)  ];
                       
    %*******************************************************************
    % ETAPE II : Construction des matrices de changement de base (BS,BD)
    %*******************************************************************
    
    %**************************
    % Base des donnees sources
    %**************************
    
    %  SourcePt2-> o      (origine)
    %                         0-------y
    %                         |          o <-SourcePt3
    %                         |
    %                         |
    %                         x
    %                         
    %                         o <-SourcePt1
    %
    % NOTE: POUR DE MEILLEURS RESULTATS:
    %                         SourcePt1=NAS, SourcePt2=LPA, SourcePt3=RPA
    SourcePt1 = [ matSourceRefs(1,1), matSourceRefs(1,2), matSourceRefs(1,3) ]; %Pt1=[x1,y1,z1]
    SourcePt2 = [ matSourceRefs(2,1), matSourceRefs(2,2), matSourceRefs(2,3) ]; %Pt2=[x2,y2,z2]
    SourcePt3 = [ matSourceRefs(3,1), matSourceRefs(3,2), matSourceRefs(3,3) ]; %Pt3=[x3,y3,z3]
    Origine_BS = [mean([SourcePt2(1), SourcePt3(1)]), ... %x
                  mean([SourcePt2(2), SourcePt3(2)]), ... %y
                  mean([SourcePt2(3), SourcePt3(3)])];    %z
    vSrc_a = SourcePt1-Origine_BS;
    vSrc_b = SourcePt3-Origine_BS;
    vSrc_c = SourcePt2-Origine_BS;
    
    %Vecteur u: 1er vecteur de la base orthonormee commune
    vSrc_u = vSrc_a / (vSrc_a*vSrc_a')^(1/2);
    
    %Vecteur v: 2e vecteur de la base orthonormee commune
    PuC = (vSrc_c*vSrc_u')*vSrc_u; %Projection de vSrc_c sur vSrc_u
    vSrc_v = (vSrc_c-PuC) / ((vSrc_c-PuC)*(vSrc_c-PuC)')^(1/2);
    
    %Vecteur w: 3e vecteur de la base orthonormee commune
    %Produit scalaire 'Right Hand'
    vSrc_w = [ vSrc_u(2)*vSrc_v(3)-vSrc_u(3)*vSrc_v(2), -(vSrc_u(1)*vSrc_v(3)-vSrc_u(3)*vSrc_v(1)), vSrc_u(1)*vSrc_v(2)-vSrc_u(2)*vSrc_v(1) ];
    if( exist( 'bIsSrcSystemRH' ) )
    	if( ~bIsSrcSystemRH )
            %Correction du produit scalaire si 'Left Hand'
            vSrc_w = -vSrc_w;
        end
    end
    
    %Base du systeme source (Base constituee des pts de reference)
    BS = [ vSrc_u(1), vSrc_v(1), vSrc_w(1); ...
           vSrc_u(2), vSrc_v(2), vSrc_w(2); ...
           vSrc_u(3), vSrc_v(3), vSrc_w(3) ];
       
    %*********************************
    % Base des donnees de destination
    %*********************************
    
    %    DestPt2-> o      (origine)   
    %                         0-------y
    %                         |          o <-DestPt3
    %                         |
    %                         |
    %                         x
    %                         
    %                         o <-DestPt1
    %
    % NOTE: POUR DE MEILLEURS RESULTATS:
    %                         SourcePt1=NAS, SourcePt2=LPA, SourcePt3=RPA
    DestPt1 = [ matDestinationRefs(1,1), matDestinationRefs(1,2), matDestinationRefs(1,3) ]; %Pt1=[x1,y1,z1]
    DestPt2 = [ matDestinationRefs(2,1), matDestinationRefs(2,2), matDestinationRefs(2,3) ]; %Pt2=[x2,y2,z2]
    DestPt3 = [ matDestinationRefs(3,1), matDestinationRefs(3,2), matDestinationRefs(3,3) ]; %Pt3=[x3,y3,z3]
    Origine_BD = [mean([DestPt2(1), DestPt3(1)]), ... %x
                  mean([DestPt2(2), DestPt3(2)]), ... %y
                  mean([DestPt2(3), DestPt3(3)])];    %z
    vDest_a = DestPt1-Origine_BD;
    vDest_b = DestPt3-Origine_BD;
    vDest_c = DestPt2-Origine_BD;
    
    %Vecteur u: 1er vecteur de la base orthonormee commune
    vDest_u = vDest_a / (vDest_a*vDest_a')^(1/2);
   
    %Vecteur v: 2e vecteur de la base orthonormee commune
    PuC = (vDest_c*vDest_u')*vDest_u; %Projection de vDest_c sur vDest_u
    vDest_v = (vDest_c-PuC) / ((vDest_c-PuC)*(vDest_c-PuC)')^(1/2);
    
    %Vecteur w: 3e vecteur de la base orthonormee commune
    %Produit scalaire 'Right Hand'
    vDest_w = [ vDest_u(2)*vDest_v(3)-vDest_u(3)*vDest_v(2), -(vDest_u(1)*vDest_v(3)-vDest_u(3)*vDest_v(1)), vDest_u(1)*vDest_v(2)-vDest_u(2)*vDest_v(1) ];
    if( exist( 'bIsDestSystemRH' ) )
    	if( ~bIsDestSystemRH )
            %Correction du produit scalaire si 'Left Hand'
            vDest_w = -vDest_w;
        end
    end
    
    %Base du systeme de destination (Base constituee des pts de reference)
    BD = [ vDest_u(1), vDest_v(1), vDest_w(1); ...
           vDest_u(2), vDest_v(2), vDest_w(2); ...
           vDest_u(3), vDest_v(3), vDest_w(3) ];
       
    %*******************************************************************
    % ETAPE III : Composition de transformation (translations, changements de bases, homotheties, ...)
    %*******************************************************************

%     Echelle_S = (vSrc_a*vSrc_a')^0.5+(vSrc_b*vSrc_b')^0.5+(vSrc_c*vSrc_c')^0.5
%     Echelle_D = (vDest_a*vDest_a')^0.5+(vDest_b*vDest_b')^0.5+(vDest_c*vDest_c')^0.5
%     FacteurAgrandissement = Echelle_D/Echelle_S;
%     
%     disp('**************************************************************');
%     disp('* Systeme Source:')
%     disp( sprintf( '* Physical Dimensions = [ %.5f %.5f %.5f ]', vSourcePhysicalDim(1), vSourcePhysicalDim(2), vSourcePhysicalDim(3) ) );
%     disp( sprintf( '* Dist LPA-RPA = %.4f cm', 100*((matSourceRefs(2,:)-matSourceRefs(3,:))*(matSourceRefs(2,:)-matSourceRefs(3,:))')^0.5 ) );
%     disp('*');
%     disp('* Systeme Destination:')
%     disp( sprintf( '* Physical Dimensions = [ %.5f %.5f %.5f ]', vDestPhysicalDim(1), vDestPhysicalDim(2), vDestPhysicalDim(3) ) );
%     disp( sprintf( '* Dist LPA-RPA = %.4f cm', 100*((matDestinationRefs(2,:)-matDestinationRefs(3,:))*(matDestinationRefs(2,:)-matDestinationRefs(3,:))')^0.5 ) );
%     disp('*');
%     disp( sprintf('* FacteurAgrandissement: %.4f', FacteurAgrandissement ) );
%     disp('**************************************************************');
%     
%     disp('* Systeme Source:')
%     matSourceRefs
%     disp('* Systeme Destination:')
%     matDestinationRefs
    
    %Inversion de Matrice (matrice orthonormee transposee)
    iBD = BD';

    %Matrice de decomposition en composantes de la base source
    T1 = [ BS(1,:), 0;...
           BS(2,:), 0;...
           BS(3,:), 0;...
           0, 0, 0, 1 ];
    
    %Matrice de recomposition a partir des composantes de la base de destination
    T2 = [ iBD(1,:), 0;...
           iBD(2,:), 0;...
           iBD(3,:), 0;...
           0, 0, 0, 1 ];
    
    %Composition de transformation
    matCoReg = T1*T2;
    