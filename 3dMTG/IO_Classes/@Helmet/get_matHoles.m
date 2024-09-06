%**************************************************************************
% FONCTION : HELMET/get_matHoles
%
% INTRANTS : obj          -> objet Helmet
%
% EXTRANTS : matrice d'infos sur les trous du casque. 
%            Dim1 : Suffixe numérique de l'idnetificateur de trou*
%            Dim2 : Préfixe alphabétique de l'idnetificateur de trou*
%
% *L'identificateur de trou est selon le format suiv -> ex: 'GB13'
%            
% UTILITE : Fonction permettant d'accéder à la matrice de
%           d'infos sur les trous du casque. 
%**************************************************************************
function matHoles = get_matHoles(obj)

    for( i=1:size(obj.mat_ilHoles,1) )
        for( j=1:size(obj.mat_ilHoles,2) )
            il = obj.mat_ilHoles(i,j);
            if(il)
                matHoles(i,j).Coord = obj.v_Holes(il).Coord;
                matHoles(i,j).RegionFlag = obj.v_Holes(il).Type == 400;
                matHoles(i,j).IsDetect = obj.v_Holes(il).CanBeDet;
            else
                matHoles(i,j).Coord.x = 0.0;
                matHoles(i,j).Coord.y = 0.0;
                matHoles(i,j).Coord.z = 0.0;
                matHoles(i,j).RegionFlag = 0;
                matHoles(i,j).IsDetect = 0;
            end
        end
    end