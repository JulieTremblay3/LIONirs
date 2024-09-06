%**************************************************************************
% FONCTION : HELMET/get_Mtg_NbSrcPeriods
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Structure d'information sur le montage du casque
%            
% UTILITE : Donne le nombre maximal de sources aux alentours d'un meme
%           detecteur (dist<distContamination).
%**************************************************************************
function nPeriods = get_Mtg_NbSrcPeriods(obj)

    nPeriods = floor( obj.Mtg_Data.Gen_Params.NbSrcPerBank ...
                      /obj.Mtg_Data.Gen_Params.Nb_Longueur_Onde);
    