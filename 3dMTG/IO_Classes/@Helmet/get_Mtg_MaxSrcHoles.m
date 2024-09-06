%**************************************************************************
% FONCTION : HELMET/get_Mtg_MaxSrcHoles
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Structure d'information sur le montage du casque
%            
% UTILITE : Donne le nombre de fibres (ou de regrouppement de fibres par 
%           periode)qui s'illuminent en meme temps.
%**************************************************************************
function nMaxHoles = get_Mtg_MaxSrcHoles(obj)
    nMaxHoles = floor( obj.Mtg_Data.Gen_Params.Nb_Sources_Disp ...
                     / obj.Mtg_Data.Gen_Params.Nb_Longueur_Onde );