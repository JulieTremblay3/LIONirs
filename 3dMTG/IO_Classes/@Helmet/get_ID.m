%**************************************************************************
% FONCTION : HELMET/get_matHolesLoc
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Identificateur numérique du casque (numéro, ex: 0000002)
%            
% UTILITE : Fonction retournant le numéro d'identification d'un casque.
%**************************************************************************
function ID_Casque = get_ID(obj)
    ID_Casque = obj.ID;