%**************************************************************************
% 
% Classe : Helmet            Méthode : ComputeNeighborhood()
%
% Intrants : Helm         -> Objet Helmet
%
% Extrants : Objet Helmet modifié
%
% Commentaires : Méthode permettant de bâtir le graphe de voisinage 
%                des trous (individuellement)
%
%**************************************************************************
function Helm = ComputeNeighborhood( Helm )

    SRC_DET_IO = get_PtTypeNo(Helm, 'NormalHole' );
    
    v_Holes = get_vHoles( Helm ); %Matrice d'information sur les trous
    
    %Pour le waitbar
    NbHolesTotal = numel( find( [v_Holes.Type] == SRC_DET_IO ) );
    if( ~NbHolesTotal )
        return;
    end
    
    NbHolesCompleted = 0;
    hWaitbar = waitbar(0,'Loading data...');
    vToBeSorted = zeros(NbHolesTotal,2);
    
    % Preallocation
    v_Near(NbHolesTotal).p = 0;
    v_Near(NbHolesTotal).Dist = 0;
    for( i = 1:numel(v_Holes) )
        v_Near(i).p = 0;
        v_Near(i).Dist = 0;
    end
    
    %Pour chaque trou du casque...
    for( pc = 1:numel(v_Holes) )
        if( v_Holes(pc).Type == SRC_DET_IO )

            %Position du trou courant
            cPos = [v_Holes(pc).Coord.x, ...
                    v_Holes(pc).Coord.y, ...
                    v_Holes(pc).Coord.z];

            NbVoisins = 0;
            
            %Pour chaque voision possible (soit le casque au complet)
            for( pv = 1:numel(v_Holes) )

                if( v_Holes(pv).Type == SRC_DET_IO )

                    %Position du trou possiblement voisin 
                    %(recherche de voisin)
                    rPos = [v_Holes(pv).Coord.x, ...
                            v_Holes(pv).Coord.y, ...
                            v_Holes(pv).Coord.z];

                    %Distance des 2 trous (courant et recherché)
                    P1P2 = rPos-cPos;
                    Dist = (P1P2*P1P2')^0.5; %racine de Produit scalaire
                    
                    %Voisin
                    NbVoisins = NbVoisins+1;
                    vToBeSorted(NbVoisins,1) = Dist;
                    vToBeSorted(NbVoisins,2) = pv;

                end %if
            end %for

            % Tri des voisins trouvés par ordre de distance
            v_Sorted = sortrows(vToBeSorted,1);

            % Reinitialisation pour un prochain tri
            vToBeSorted = zeros(NbHolesTotal,2);

            % Preallocation
            v_Holes(pc).Neighbors.v_Near = v_Near;
%             v_Holes(pc).Neighbors.v_Near(NbHolesTotal).p = 0;
%             v_Holes(pc).Neighbors.v_Near(NbHolesTotal).Dist = 0;

            % Transfert dans le vecteur de voisinage final
            for( n=1:NbHolesTotal )
                v_Holes(pc).Neighbors.v_Near(n).Dist = v_Sorted(n,1);
                v_Holes(pc).Neighbors.v_Near(n).p    = v_Sorted(n,2);
            end
            
%             if( mod( pc, 100 ) == 0 )
%                 disp( '-------------------------------------------------------' );
%                     
%                 for( n=1:NbHolesTotal )
%                     disp( sprintf( 'Dist=%2.3f  pNeigh=%d',v_Holes(pc).Neighbors.v_Near(n).Dist*100, v_Holes(pc).Neighbors.v_Near(n).p ) );
%                 end
%             end

            %Sauvegarde du nombre de voisins pour ce trou
            v_Holes(pc).Neighbors.Nb = NbHolesTotal;

            NbHolesCompleted = NbHolesCompleted+1;
            waitbar(NbHolesCompleted/NbHolesTotal);

        end %( matH(lin,col).RegionFlag ~= RIEN )
    end %for
    
    Helm.NeighborhoodOk = true;
    
    close(hWaitbar);
    
    Helm = set_vHoles( Helm, v_Holes );
    
    