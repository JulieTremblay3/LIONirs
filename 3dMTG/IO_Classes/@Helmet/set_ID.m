%**************************************************************************
% FONCTION : HELMET/set_ID
%
% INTRANTS : obj          -> objet Helmet
%            ID_Helmet    -> Identificateur num�rique du casque (Ex:000001)
%
% EXTRANTS : Objet Helmet modifi�.
%
% UTILITE : Permet d'�tablir le num�ro d'identification du casque.
%**************************************************************************
function obj = set_Mtg(obj, ID_Helmet )

    obj.Mtg_Data = ID_Helmet;