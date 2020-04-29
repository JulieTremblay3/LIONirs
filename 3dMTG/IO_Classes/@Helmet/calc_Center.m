%**************************************************************************
% FONCTION : HELMET/calc_Center
%
% INTRANTS : obj          -> objet Helmet
%
% EXTRANTS : casque modifie
%
%**************************************************************************
function obj = calc_Center(obj)
    
    Coords = [obj.v_Holes.Coord];
    
    minX = min( [Coords.x] );
    maxX = max( [Coords.x] );
    minY = min( [Coords.y] );
    maxY = max( [Coords.y] );
    minZ = min( [Coords.z] );
    maxZ = max( [Coords.z] );
    
    obj.v_Center = [ minX+maxX, minY+maxY, minZ+maxZ ] / 2;