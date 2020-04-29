function oPrj = set_Dig_SubjectFids( oPrj, oDig )
    
    if( isa( oDig, 'Digitization' ) )
        oPrj.m_oDig_SubjectFids = oDig;
    else
        disp('set_Dig_CompleteHelmet:Error: not a Digitization object');
    end
