%**************************************************************************
% FONCTION : HELMET/get_Mtg_PeriodMultiplicity
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Structure d'information sur le montage du casque
%            
% UTILITE : Donne le nombre de fibres (ou de regrouppement de fibres par 
%           periode)qui s'illuminent en meme temps.
%**************************************************************************
function nMultiplicity = get_Mtg_PeriodMultiplicity(obj)
    nMultiplicity = obj.Mtg_Data.Gen_Params.NbBanks;