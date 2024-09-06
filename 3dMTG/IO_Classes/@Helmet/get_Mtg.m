%**************************************************************************
% FONCTION : HELMET/get_Mtg
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Structure d'information sur le montage du casque
%            
% UTILITE : Fonction retournant la strcuture d'information sur le montage
%           a faire sur le casque. N.B: Chaque montage n'est valide que
%           pour un casque donné.
%**************************************************************************
function Mtg_Data = get_Mtg(obj)
    Mtg_Data = obj.Mtg_Data;
    Mtg_Data.v_pSrc = Mtg_Data.v_pSrc( find(Mtg_Data.v_pSrc) );
    Mtg_Data.v_pDet = Mtg_Data.v_pDet( find(Mtg_Data.v_pDet) );
    