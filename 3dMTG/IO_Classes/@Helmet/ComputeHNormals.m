%**************************************************************************
%   Fonction permettant de calculer les normales des trous du casque
%**************************************************************************
function oHelmet = ComputeHNormals( oHelmet, PlaneRange )

%Etablir le voisinage
    if( ~oHelmet.NeighborhoodOk )
        oHelmet = ComputeNeighborhood( oHelmet );
    end
    
    %Acces a la matrice de trous
    v_Holes = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );
    
    %Centre de la tete: approximativement la moyenne(LPA,RPA)
    vHead_Center = [ mean([sMtg.matFiducials(2,1), sMtg.matFiducials(3,1)]), ...
                     mean([sMtg.matFiducials(2,2), sMtg.matFiducials(3,2)]), ...
                     mean([sMtg.matFiducials(2,3), sMtg.matFiducials(3,3)]) ];
                
    %Mettre les coordonnees sous forme de matrice de vecteurs
    %Ici, toutes les coordonnees sont translatees afin que le point (0,0,0)
    %se situe au centre de la tete. Cette translation ne change rien a
    %l'ANGLE des normales, mais peut possiblement changer la DIRECTION des
    %normales. Pourquoi? Parce que la methode d'approximation de plan par
    %les moindres carrees utilisee donne une normale dont la direction est
    %determinee par la position de l'origine: toutes les normales ont donc
    %une direction qui s'eloigne du centre. Bref, on veut que toutes les
    %normales du casque s'eloignent de la tete, et en translatant toutes
    %les positions afin de mettre l'origine au centre de la tete, on
    %s'assure de la validite de la direction des normale. Toutefois, je me
    %repete, cette translation de l'origine ne change pas l'ANGLE des
    %normales mais seulement la direction (+/-) des normales.
    %
    %Cette etape vient du fait qu'il soit possible de digitaliser les
    %casques avec le 'Calibration Block' sous Brainsight et, dans un tel
    %cas, l'origine n'est pas dans le casque mais un peu en avant du
    %casque. Donc, cette etape permet d'utiliser les digitalisations faites
    %en format LOCATOR avec le CALIBRATION BLOCK.
    for( Pos=1:numel(oHelmet.v_Holes) );
        v_Coords_xyz(Pos,:) = [v_Holes(Pos).Coord.x, v_Holes(Pos).Coord.y, v_Holes(Pos).Coord.z]-vHead_Center;
    end
    
 
    %Pour tous les trous
    for( p=1:numel(oHelmet.v_Holes) );
        
        %Trou existant
        if( v_Holes(p).Type == get_PtTypeNo(oHelmet, 'NormalHole' ) )

            oHelmet.v_Holes(p).Normal.x = 0;
            oHelmet.v_Holes(p).Normal.y = 0;
            oHelmet.v_Holes(p).Normal.z = 0;

            %Trouver les voisins inclus dans un perimetre de 0.05 m (5 cm)
            PosNeighbors = find( [v_Holes(p).Neighbors.v_Near.Dist] < 0.05 );

            if( ~isempty( PosNeighbors ) )
 
                %Convertir l'indice de voisin (Pos) en indice de trou (p)
                S.type = '()';
                S.subs = {PosNeighbors};
                v_pNei = subsref([v_Holes(p).Neighbors.v_Near.p],S);
                v_pNei = v_pNei( find(v_pNei) );
                if ( numel(v_pNei) > 2 )& numel(v_Holes) > 100
                    
                    %matrice de vecteurs de positions de voisins
                    v_NeiCoords = v_Coords_xyz(v_pNei,:);

                    %ALGORITME D'APPROXIMATION DE PLAN PAR LES MOINDRES
                    %CARRES (forme standard)
                    sumxy = v_NeiCoords(:,1)'*v_NeiCoords(:,2); %produit scalaire x avec y = sum(x*y);
                    sumxz = v_NeiCoords(:,1)'*v_NeiCoords(:,3); %produit scalaire x avec y = sum(x*y);
                    sumyz = v_NeiCoords(:,2)'*v_NeiCoords(:,3); %produit scalaire x avec y = sum(x*y);
                    sumxx = v_NeiCoords(:,1)'*v_NeiCoords(:,1); %produit scalaire x avec x = sum(x^2);
                    sumyy = v_NeiCoords(:,2)'*v_NeiCoords(:,2); %produit scalaire y avec y = sum(y^2);
                    sumzz = v_NeiCoords(:,3)'*v_NeiCoords(:,3); %produit scalaire z avec z = sum(z^2);

                    A = [ sumxx, sumxy, sumxz;
                          sumxy, sumyy, sumyz;
                          sumxz, sumyz, sumzz ];

                    b = [ -sum( v_NeiCoords(:,1) );  % Sommation des xi
                          -sum( v_NeiCoords(:,2) );  % Sommation des yi
                          -sum( v_NeiCoords(:,3) ) ];% Sommation des zi

                    n = inv(A)*b;

                    d = 1/(n'*n)^(1/2);

                    oHelmet.v_Holes(p).Normal.x = -d*n(1);
                    oHelmet.v_Holes(p).Normal.y = -d*n(2);
                    oHelmet.v_Holes(p).Normal.z = -d*n(3);
                    %FIN ALGORITME D'APPROXIMATION DE PLAN PAR LES MOINDRES
                    %CARRES (forme standard)
                elseif  sum(abs(v_Coords_xyz(:,3)))>0 %3d plane  go to the center of the head (0,0,0)
                    norm = sqrt(v_Holes(p).Coord.x^2+v_Holes(p).Coord.y^2+v_Holes(p).Coord.z^2);                    
                    oHelmet.v_Holes(p).Normal.x = v_Holes(p).Coord.x/norm; 
                    oHelmet.v_Holes(p).Normal.y = v_Holes(p).Coord.y/norm;
                    oHelmet.v_Holes(p).Normal.z = v_Holes(p).Coord.z/norm;
                else %2d
                    oHelmet.v_Holes(p).Normal.x = 0.001;
                    oHelmet.v_Holes(p).Normal.y = 0.001;
                    oHelmet.v_Holes(p).Normal.z = 0.1;
                end
            end
        end % for
    end %for 
    1;
    %oHelmet = FilterNormals( oHelmet );

    