%**************************************************************************
% FONCTION : HELMET/set_ID
%
% INTRANTS : obj          -> objet Helmet
%            ID_Helmet    -> Identificateur numérique du casque (Ex:000001)
%
% EXTRANTS : Objet Helmet modifié.
%
% UTILITE : Permet d'établir le numéro d'identification du casque.
%**************************************************************************
function obj = set_Mtg(obj, ID_Helmet )

    obj.Mtg_Data = ID_Helmet;