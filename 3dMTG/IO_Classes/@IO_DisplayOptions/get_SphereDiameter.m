% Permet de connaitre l'etat d'un option d'affichage
%
% Intrants: oDispOpt  : Objet IO_DisplayOptions
%           strField  : diametre
%
% Extrants: bIsActive : 0 = deasactive
%                       1 = active
%
function strField = get_SphereDiameter( oDispOpt )

   strField = oDispOpt.SphereDiameter ;
    
