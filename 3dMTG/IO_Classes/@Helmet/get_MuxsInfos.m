% mat_Pts = [ [x1, x2, x3, x4];
%             [y1, y2, y3, y4];
%             [z1, z2, z3, z4] ]
%
% mat_Norms = [ [x1, x2, x3, x4];
%               [y1, y2, y3, y4];
%               [z1, z2, z3, z4] ]
%
% mat_pMuxCouples = [ pDet1, pSrc1;
%                     pDet1, pSrc2;
%                     pDet2, pSrc1;
%                     pDet2, pSrc2;
%                     pDet2, pSrc3 ]
function [ v_DetAreaCovered ]= get_MuxsInfos( Helm, mat_pMuxCouples )

    v_DetAreaCovered = [];
    
    vHoles = get_vHoles( Helm );
    sMtg = get_Mtg( Helm );
    
    mat_MP1   = zeros(3, size(mat_pMuxCouples,1) ); %Matrices de vecteurs-colonnes
    mat_MP2   = zeros(3, size(mat_pMuxCouples,1) );
    mat_Pts   = zeros(3, size(mat_pMuxCouples,1) );
    mat_Norms = zeros(3, size(mat_pMuxCouples,1) );
    %mat_fMP1  = zeros(2, size(mat_pMuxCouples,1) ); 
    %mat_fMP2  = zeros(2, size(mat_pMuxCouples,1) );
    
    %Mise en forme matricielle
    for( iMux=1:size(mat_pMuxCouples,1) )
        
        pDet = mat_pMuxCouples(iMux,1);
        pSrc = mat_pMuxCouples(iMux,2);
        mat_MP1(:,iMux) = [ vHoles(pDet).Coord.x; vHoles(pDet).Coord.y; vHoles(pDet).Coord.z ];
        mat_MP2(:,iMux) = [ vHoles(pSrc).Coord.x; vHoles(pSrc).Coord.y; vHoles(pSrc).Coord.z ];
        
        %Vecteur de normales de muxs.
        mat_Pts(:,iMux) = (mat_MP1(:,iMux)+mat_MP2(:,iMux))/2;
        mat_Norms(:,iMux) = [ vHoles(pDet).Normal.x + vHoles(pSrc).Normal.x;
                              vHoles(pDet).Normal.y + vHoles(pSrc).Normal.y;
                              vHoles(pDet).Normal.z + vHoles(pSrc).Normal.z; ];
        mat_Norms(:,iMux) = mat_Norms(:,iMux)/( mat_Norms(:,iMux)'*mat_Norms(:,iMux) )^(0.5);
    end
    
    %Pour chaque region d'interet    
    for( r=1:size(sMtg.mat_RegionNorms,1) )
        
        v_iPtsFound = find( sMtg.mat_RegionNorms(r,:,1) | sMtg.mat_RegionNorms(r,:,2) | sMtg.mat_RegionNorms(r,:,3) );
        
        if( ~isempty(v_iPtsFound) )
            
            disp('Contruction de matrice de projection');
            
            %Utiliser la normale du premier point (elles sont toutes
            %pareilles pour la region)
            vn   = [ sMtg.mat_RegionNorms(r,v_iPtsFound(1),1);
                     sMtg.mat_RegionNorms(r,v_iPtsFound(1),2);
                     sMtg.mat_RegionNorms(r,v_iPtsFound(1),3) ];
            vz = [ 0; 0; 1; ];

            % Produit vectoriel de n et z
            A = vn(2)*vz(3)-vn(3)*vz(2);
            B = vn(3)*vz(1)-vn(1)*vz(3);
            C = vn(1)*vz(2)-vn(2)*vz(1);

            tetha = 2*atan2( find_distance( vn(1)-vz(1), vn(2)-vz(2), vn(3)-vz(3), 0,0,0 )/2, ...
                              find_distance( vn(1)+vz(1), vn(2)+vz(2), vn(3)+vz(3), 0,0,0 )/2 );

            %Matrice de rotation autour du vecteur v=(A,B,C)
            matRot = [ (1-A^2)*cos(tetha)+A^2,          A*B*(1-cos(tetha))-C*sin(tetha), A*C*(1-cos(tetha))+B*sin(tetha); ...
                       A*B*(1-cos(tetha))+C*sin(tetha), (1-B^2)*cos(tetha)+B^2,          B*C*(1-cos(tetha))-A*sin(tetha); ...
                       A*C*(1-cos(tetha))-B*sin(tetha), B*C*(1-cos(tetha))+A*sin(tetha), (1-C^2)*cos(tetha)+C^2 ];
                   
            matProjPXY = [ 1, 0, 0;
                           0, 1, 0;
                           0, 0, 0 ];
            
            %disp('Application de la transformation');
            matFlatten = matProjPXY*matRot;
            
            v_iPtsFound(length(v_iPtsFound)+1) = v_iPtsFound(1);
            
            mat_vXYZ = [ sMtg.mat_RegionBorder(r,v_iPtsFound,1); ...
                         sMtg.mat_RegionBorder(r,v_iPtsFound,2); ...
                         sMtg.mat_RegionBorder(r,v_iPtsFound,3)  ];
            
            mat_vXYZ = matFlatten*mat_vXYZ;
            
            mat_fRiP1 = [ mat_vXYZ(1,1:size(mat_vXYZ,2)-1); ...
                          mat_vXYZ(2,1:size(mat_vXYZ,2)-1) ]';
            
            mat_fRiP2 = [ mat_vXYZ(1,2:size(mat_vXYZ,2)  ); ...
                          mat_vXYZ(2,2:size(mat_vXYZ,2)  ) ]';
            
            mat_fMP1 = matFlatten*mat_MP1;
            mat_fMP2 = matFlatten*mat_MP2;
            
            plot3( mat_vXYZ(1,:), mat_vXYZ(2,:), mat_vXYZ(3,:), 'y--' );
            
            %Semble ok
            %for( iTest=1:size(mat_fMP1,2) )
            %    plot3( [ mat_fMP1(1,iTest),mat_fMP2(1,iTest) ], ...
            %           [ mat_fMP1(2,iTest),mat_fMP2(2,iTest) ], ...
            %           [ mat_fMP1(3,iTest),mat_fMP2(3,iTest) ], ...
            %           'r-' );
            %end
            %break;
            
            
            
        end
    end
                     
            
            
%                     
%             %Pour chaque point:
%             for( iPt=1:length(mat_Pts) )
%                 
%                 %Matrice-colonne de vecteurs normaux
%                 v_nPt =  mat_Norms(:,iPt);
%         
%         
%         
%     %initialisation de la matrice de presence dans la region d'interet
%     v_RegId_Elemt = zeros( 1, length(mat_Norms) );
%     
% 	%Pour chaque region d'interet    
%     for( r=1:size(sMtg.mat_RegionNorms,1) )
%         
%         v_iPtsFound = find( sMtg.mat_RegionNorms(r,:,1) | sMtg.mat_RegionNorms(r,:,2) | sMtg.mat_RegionNorms(r,:,3) );
%         
%         Color =  [r*0.25,1-r*0.25,1-r*0.25];
%         if( ~isempty(v_iPtsFound) )
%             
%             %disp('Contruction de matrice de projection');
%             
%             %Utiliser la normale du premier point (elles sont toutes
%             %pareilles
%             vn   = [ sMtg.mat_RegionNorms(r,v_iPtsFound(1),1);
%                      sMtg.mat_RegionNorms(r,v_iPtsFound(1),2);
%                      sMtg.mat_RegionNorms(r,v_iPtsFound(1),3) ];
%             vz = [ 0; 0; 1; ];
% 
%             % Produit vectoriel de n et z
%             A = vn(2)*vz(3)-vn(3)*vz(2);
%             B = vn(3)*vz(1)-vn(1)*vz(3);
%             C = vn(1)*vz(2)-vn(2)*vz(1);
% 
%             tetha = 2*atan2( find_distance( vn(1)-vz(1), vn(2)-vz(2), vn(3)-vz(3), 0,0,0 )/2, ...
%                               find_distance( vn(1)+vz(1), vn(2)+vz(2), vn(3)+vz(3), 0,0,0 )/2 );
% 
%             %Matrice de rotation autour du vecteur v=(A,B,C)
%             matRot = [ (1-A^2)*cos(tetha)+A^2,          A*B*(1-cos(tetha))-C*sin(tetha), A*C*(1-cos(tetha))+B*sin(tetha); ...
%                        A*B*(1-cos(tetha))+C*sin(tetha), (1-B^2)*cos(tetha)+B^2,          B*C*(1-cos(tetha))-A*sin(tetha); ...
%                        A*C*(1-cos(tetha))-B*sin(tetha), B*C*(1-cos(tetha))+A*sin(tetha), (1-C^2)*cos(tetha)+C^2 ];
%                    
%             matProjPXY = [ 1, 0, 0;
%                            0, 1, 0;
%                            0, 0, 0 ];
%             
%             %disp('Application de la transformation');
%             matFlatten = matProjPXY*matRot;
%             
%             v_iPtsFound(length(v_iPtsFound)+1) = v_iPtsFound(1);
%             
%             mat_vXYZ = [ sMtg.mat_RegionBorder(r,v_iPtsFound,1); ...
%                          sMtg.mat_RegionBorder(r,v_iPtsFound,2); ...
%                          sMtg.mat_RegionBorder(r,v_iPtsFound,3)  ];
%             
%             mat_vXYZ = matFlatten*mat_vXYZ;
%             
%             vP1 = [ mat_vXYZ(1,1:size(mat_vXYZ,2)-1); ...
%                     mat_vXYZ(2,1:size(mat_vXYZ,2)-1) ]';
%             
%             vP2 = [ mat_vXYZ(1,2:size(mat_vXYZ,2)  ); ...
%                     mat_vXYZ(2,2:size(mat_vXYZ,2)  ) ]';
%             
%             %plot3( mat_vXYZ(1,:), mat_vXYZ(2,:), mat_vXYZ(3,:), 'y--' );
%             
%             %Ca a l'air ben bon ca ^
%             %*************************************************************
%             
%             
%             %Pour chaque point:
%             for( iPt=1:length(mat_Pts) )
%                 
%                 %Matrice-colonne de vecteurs normaux
%                 v_nPt =  mat_Norms(:,iPt);
% 
%                          
%                 if( v_nPt'*vn > 0 )
%                     
%                     v_cPt =  mat_Pts(:,iPt);
%                     v_cPt = matFlatten*v_cPt;
%                     
%                     %*************************************************************
%                     %Ici, nous avons la projection de la region et du point
%                     %*************************************************************
%                     Found = find(   (  ( vP1(:,2) < v_cPt(2) & vP2(:,2) > v_cPt(2) ) ...
%                                       |( vP1(:,2) > v_cPt(2) & vP2(:,2) < v_cPt(2) ) ) ...
%                                       &( vP1(:,1) < v_cPt(1) | vP2(:,1) < v_cPt(1) ) );
%                     
%                     %Recherche de croisements invalides
%                     for( jLine=1:length(Found) )
%                         Pos = Found(jLine);
%                         %plot3( [vP1(Pos,1),vP2(Pos,1)], [vP1(Pos,2),vP2(Pos,2)],[0,0], 'r' );
%                     
%                         if( ~(vP1(Pos,1) < v_cPt(1) && vP2(Pos,1) < v_cPt(1)) )
%                             
%                             DeltaX_P1P2 = abs( vP2(Pos,1)-vP1(Pos,1) );
%                             DeltaY_P1P2 = abs( vP2(Pos,2)-vP1(Pos,2) );
%                             
%                             if( vP1(Pos,2) < vP2(Pos,2) )
%                                 LowerPt = vP1(Pos,:);
%                                 HigherPt = vP2(Pos,:);
%                             else
%                                 LowerPt = vP2(Pos,:);
%                                 HigherPt = vP1(Pos,:);
%                             end
%                             
%                             DeltaY_Cur_Fr = abs( v_cPt(2)-LowerPt(2) );
%                             X_Intersect = LowerPt+(HigherPt-LowerPt)*(DeltaY_Cur_Fr/DeltaY_P1P2);
%                             if( X_Intersect(1) > v_cPt(1) )
%                                 %Invalidation des croisements bidon
%                                 %plot3( [vP1(Pos,1),vP2(Pos,1)], [vP1(Pos,2),vP2(Pos,2)],[0,0], 'b' );
%                                 Found(jLine) = 0;
%                             end
%                         end
%                     end
%                     
%                     %Retrait des croisements invalides
%                     Found = find(Found);
%                     
%                     if( mod( length(Found), 2 ) )
%                         %Inside
%                         v_RegId_Elemt(iPt) = v_RegId_Elemt(iPt) + 2^(iPt-1);
%                         %plot3( v_cPt(1), v_cPt(2), v_cPt(3), 'r*' );
%                     end
%                 end
%             end
%         end
%     end


