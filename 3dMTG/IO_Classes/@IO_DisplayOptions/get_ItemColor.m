%function [rVal,gVal,bVal] = get_ItemColor( oObj, strField )
function vColor = get_ItemColor( oDispOpt, strField )
    
    matColor = getfield( oDispOpt.ColorSchemes, strField );
    
    if( oDispOpt.v_bItemChk( oDispOpt.Subscripts.App_InvertedColorScheme ) )
        vColor = matColor(2,:);
    else
        vColor = matColor(1,:);
    end