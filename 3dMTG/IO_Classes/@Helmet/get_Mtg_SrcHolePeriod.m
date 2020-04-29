%**************************************************************************
% FONCTION : HELMET/get_Mtg_SrcHolePeriod
%
% INTRANTS : obj     -> objet Helmet
%            pSrc    -> indice lineaire du trou
%
% EXTRANTS : Periode d'illumination de source.
%            
% UTILITE : Retourne le numero de periode (1-NbPeriods). 
%           Retourne 0 si aucune source n'est attribuee au trou.
%**************************************************************************
function PeriodNo = get_Mtg_SrcHolePeriod( obj, pSrc )
    
    if( pSrc > numel(obj.Mtg_Data.v_HolesMtg) )
        disp( sprintf( 'Error: get_Mtg_SrcHolePeriod(), pSrc:%d, numel(vHolesMtg):%d', ...
             pSrc, numel(obj.Mtg_Data.v_HolesMtg) ) );
         PeriodNo = -1;
    else
        PeriodNo = get_Mtg_SrcFibersPeriod( obj, obj.Mtg_Data.v_HolesMtg(pSrc) );
    end