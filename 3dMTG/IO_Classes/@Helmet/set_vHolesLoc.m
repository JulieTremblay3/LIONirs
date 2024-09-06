%**************************************************************************
% FONCTION : HELMET/set_vHolesLoc 
%
% INTRANTS : obj          -> objet Helmet
%            new_vHoles   -> nouveau vecteur de coordonnées de trous de 
%                            casque.
%
% EXTRANTS : Objet modifié.
%
% UTILITE : Fonction permettant d'attribuer une nouvelle
%           matrice de coordonnées de trous de casque.
%**************************************************************************
function obj = set_vHolesLoc(obj, new_vHoles)
    
    %Reinitialisation
    obj = Clear_vHoles(obj);
    obj.v_HolesReady = false;
    
    for( p=1:numel(new_vHoles) )
        obj.v_Holes(p).Coord.x = new_vHoles(p).x;
        obj.v_Holes(p).Coord.y = new_vHoles(p).y;
        obj.v_Holes(p).Coord.z = new_vHoles(p).z;
    end
    
    %Should compute norms, neighbors, etc. here.
