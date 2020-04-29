%**************************************************************************
% FONCTION : HELMET/get_Mtg_SrcHolePeriod
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : Structure d'information sur le montage du casque
%            
% UTILITE : Retourne le numero de periode (1-NbPeriods). 
%           Retourne 0 si aucune source n'est attribuee au trou.
%**************************************************************************
function PeriodNo = get_Mtg_SrcHolePeriod( obj, pSrc )
    

    SrcCode = ( obj.Mtg_Data.v_HolesMtg 
    

    PeriodNo = floor( obj.Mtg_Data.Gen_Params.Nb_Sources_Disp ...
                     / obj.Mtg_Data.Gen_Params.Nb_Longueur_Onde );
                 
                 
    obj.Mtg_Data.v_HolesMtg(1)             
                 obj.Mtg_Data.v_HolesMtg(1) = 0; % 0XX0XX (src1src2) ou XX000000 det1
    