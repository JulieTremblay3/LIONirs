function oPrj = set_Dig_CompleteHelmet( oPrj, oDig )
    
    if( isa( oDig, 'Digitization' ) )
        oPrj.m_oDig_CompleteHelmet = oDig;
    else
        disp('set_Dig_CompleteHelmet:Error: not a Digitization object');
    end
