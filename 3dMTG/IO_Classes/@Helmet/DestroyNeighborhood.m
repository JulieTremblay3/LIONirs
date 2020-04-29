%**************************************************************************
% 
% Classe : Helmet            Méthode : DestroyNeighborhood()
%
% Intrants : Helm         -> Objet Helmet
%
% Extrants : Objet Helmet modifié
%
% Commentaires : Méthode permettant de bâtir le graphe de voisinage 
%                des trous (individuellement)
%
%**************************************************************************
function Helm = DestroyNeighborhood( Helm )

    SRC_DET_IO = 400; %Source ou detecteur
    
    v_Holes = get_vHoles( Helm ); %Matrice d'information sur les trous

            
    %Pour chaque trou du casque...
    for( pc = 1:numel(v_Holes) )
        if( v_Holes(pc).Type == SRC_DET_IO )

            v_Empty(1).p=0;
            v_Empty(1).Dist=0;
            
            v_Holes(pc).Neighbors.v_Near = v_Empty;

            %Sauvegarde du nombre de voisins pour ce trou
            v_Holes(pc).Neighbors.Nb = 0;

        end %( matH(lin,col).RegionFlag ~= RIEN )
    end %for
    
    
    Helm = set_vHoles( Helm, v_Holes );
    
    