%**************************************************************************
% FONCTION : HELMET/set_matHoles 
%
% INTRANTS : obj          -> objet Helmet
%            new_vHoles -> nouvelle matrice d'information sur les trous
%
% EXTRANTS : Objet Helmet modifié.
%
% UTILITE : Fonction permettant d'attribuer un nouveau
%           vecteur d'infos de trous de casque.
%**************************************************************************
function obj = set_vHoles(obj, new_vHoles)
   
    obj.v_Holes = new_vHoles;
    
