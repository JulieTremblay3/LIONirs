% Permet d'ajuster le rayon des sphères qui affiche les trous
%
% Intrants: oDispOpt  : Objet IO_DisplayOptions
%           strField  : diametre.
%
% Extrants: oDispOpt : Objet modifie
%
%**************************************************************************
function oDispOpt = set_SphereDiameter( oDispOpt, strField )
  
   oDispOpt.SphereDiameter = strField;