% Permet de choisir les zones de l'atlas à afficher
%
% Intrants: oDispOpt  : Objet IO_DisplayOptions
%           strField  : nom de l'image ajoute.
%
% Extrants: oDispOpt : Objet modifie
%
%**************************************************************************
function oDispOpt = set_DispAtlasZone( oDispOpt, strField )
  
    oDispOpt.AtlasZone = strField;