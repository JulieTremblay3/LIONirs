%**************************************************************************
% FONCTION : HELMET/get_matHolesLoc
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Identificateur num�rique du casque (num�ro, ex: 0000002)
%            
% UTILITE : Fonction retournant le num�ro d'identification d'un casque.
%**************************************************************************
function ID_Casque = get_ID(obj)
    ID_Casque = obj.ID;