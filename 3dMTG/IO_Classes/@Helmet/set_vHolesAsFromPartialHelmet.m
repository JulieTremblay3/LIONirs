function oHelmet = set_vHolesAsFromPartialHelmet( oHelmet, v_p )
    v_p = v_p( find( (v_p > 0) & (v_p <= numel(oHelmet.v_Holes)) ) );
    
    for( iElem=1:numel(v_p) )
        oHelmet.v_Holes(v_p(iElem)).IsFromCompleteHelmet = false;
    end
    
    