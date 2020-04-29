% Permet de connaitre l'etat d'un option d'affichage
%
% Intrants: oDispOpt  : Objet IO_DisplayOptions
%           strField  : string de l'option d'affichage
%
% Extrants: bIsActive : 0 = deasactive
%                       1 = active
%
function strField = get_DispAtlasZone( oDispOpt )

   strField = oDispOpt.AtlasZone ;
    
