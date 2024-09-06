%**************************************************************************
% FONCTION : HELMET/get_NbLines
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Nombre de lignes de la matrice de trous.
%            
% UTILITE : Fonction retournant le nombre trous maximal de la matrice 
%           matHolesLoc pour toutes ses colonnes. 
%**************************************************************************
function NbLines = get_NbLines(obj)
    NbLines = size(obj.mat_ilHoles,1);