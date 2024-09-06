%**************************************************************************
% FONCTION : HELMET/get_MtgDet_SrcInRangeOfGoodData
%
% INTRANTS : obj     -> objet Helmet
%            pDet    -> indice lineaire du trou
%
% EXTRANTS : vecteur de sources a proximite
%**************************************************************************
function v_pSrcsInRange = get_MtgDet_SrcInRangeOfGoodData( obj, pDet )
    
    v_pSrcsInRange = [];
    sMtg = obj.Mtg_Data;
    for( iSrc=1:numel(sMtg.v_pSrc) )
        pSrc = sMtg.v_pSrc(iSrc);
        if(    get_DistHoles( obj, pSrc, pDet ) > sMtg.Gen_Params.Min_Dist_Mux ...
            && get_DistHoles( obj, pSrc, pDet ) < sMtg.Gen_Params.Max_Dist_Mux)
            v_pSrcsInRange = [ v_pSrcsInRange, pSrc ];
        end
    end
    