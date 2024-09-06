%**************************************************************************
% FONCTION : HELMET/get_NbColsGD
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Nombre de colonnes de trous du casque.
%            
% UTILITE : Fonction retournant le nombre de colonnes de trous d'un casque.
%           Il y a une colonne de trous qui correspond à chaque préfixe
%           alphabétique différent. Ex: GB10, GB11, GB12 sont sur la même
%           colonne de trous alors que DA10, DB10 et GE11 sont tous sur 
%           des colonnes différentes conséquemment à leur préfixe 
%           alphabétique d'identificateur différent.
%**************************************************************************
function NbColsGD = get_NbColsGD(obj)
    NbColsGD = size(obj.mat_ilHoles,2);
    