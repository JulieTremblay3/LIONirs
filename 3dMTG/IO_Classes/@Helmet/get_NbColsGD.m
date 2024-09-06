%**************************************************************************
% FONCTION : HELMET/get_NbColsGD
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Nombre de colonnes de trous du casque.
%            
% UTILITE : Fonction retournant le nombre de colonnes de trous d'un casque.
%           Il y a une colonne de trous qui correspond � chaque pr�fixe
%           alphab�tique diff�rent. Ex: GB10, GB11, GB12 sont sur la m�me
%           colonne de trous alors que DA10, DB10 et GE11 sont tous sur 
%           des colonnes diff�rentes cons�quemment � leur pr�fixe 
%           alphab�tique d'identificateur diff�rent.
%**************************************************************************
function NbColsGD = get_NbColsGD(obj)
    NbColsGD = size(obj.mat_ilHoles,2);
    