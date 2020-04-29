%**************************************************************************
% FONCTION : HELMET/get_vHolesLoc
%
% INTRANTS : obj          -> objet Helmet
%
% EXTRANTS : vecteur de coordonnées de trous du casque.
%            
%**************************************************************************
function vHoles = get_vHolesLoc(obj)

    vHoles(1,1).x = 0;
    vHoles(1,1).y = 0;
    vHoles(1,1).z = 0;
    
    for( p=1:numel(obj.v_Holes) )
        vHoles(p).x = obj.v_Holes(p).Coord.x;
        vHoles(p).y = obj.v_Holes(p).Coord.y;
        vHoles(p).z = obj.v_Holes(p).Coord.z;
    end
    