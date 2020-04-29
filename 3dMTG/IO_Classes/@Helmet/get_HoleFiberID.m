%**************************************************************************
% FONCTION : HELMET/get_HoleFiberID
%
% INTRANTS : obj          -> objet Helmet
%                         -> Indice lineaire du 1er trou
%                         -> Indice lineaire du 2e  trou
%
% EXTRANTS : Distance separant les deux trous specifies. 
%
%**************************************************************************
function strFiber = get_HoleFiberID( obj, p )

    if( p > numel(obj.Mtg_Data.v_HolesMtg) )
        disp( 'get_HoleFiberID: subscript value exceeds numel(v_HolesMtg)' ); 
        disp( sprintf( '    p:%d     numel:%d', p, numel(obj.Mtg_Data.v_HolesMtg) ) );
        strFiber = 'ERR';
        return;
    elseif( ~obj.Mtg_Data.v_HolesMtg(p) )
        strFiber = '';
        return;
    end
     
    strFiber = get_Mtg_nFib2strFib( obj, obj.Mtg_Data.v_HolesMtg(p) );