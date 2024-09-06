function strFib = ISSIDnumber2ImagincHubert(nFib)
%CONVERTI LES NUMEROS DE SOURCE ISS
% a1b2 018001 -> E1
% a3b4 020003 -> E2
% a5b6 022005 -> E3
%1st Src (XXX001-XXX016)
        Src1_Code = (mod(nFib, 1000)+1)/2;
        
        if( Src1_Code == 0 )
            strFib = sprintf( '' );
        else
            Src1_E = char('E');
            Src1_Nb = char(int2str(Src1_Code));      
            strFib = sprintf( '%c%c%c', Src1_E, Src1_Nb);
        end