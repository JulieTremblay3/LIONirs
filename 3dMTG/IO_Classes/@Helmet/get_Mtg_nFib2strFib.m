%**************************************************************************
% FONCTION : HELMET/get_Mtg_nFib2strFib
%
% INTRANTS : obj -> objet Helmet
%
% EXTRANTS : String d'identification de la fibre optique propre au system
% ISS ou IMAGINC
%            
% UTILITE : Retourne la string correspondant aux # de fibres recu en param
%**************************************************************************
function strFib = get_Mtg_nFib2strFib( obj, nFib )
    if ~isfield(obj.Mtg_Data.Gen_Params,'AcqSystem')
        obj.Mtg_Data.Gen_Params.AcqSystem = 'ISS'
    end
         
    switch obj.Mtg_Data.Gen_Params.AcqSystem
        case 'ISS' 
    %%%%%%%%% NOMENCLATURE ISS Utilisé avec 128 source et 16 detecteur
    if( nFib >= 1000000 )        
        %Detectors
        Det_Code = nFib/1000000;
        if( Det_Code == 0 )
            strFib = sprintf( '' );
        else 
            strFib = sprintf( '%c', char( int16('A')-1+Det_Code) );
        end
    else
        %Sources        
        %1st Src (XXX001-XXX016)
        Src1_Code = mod(nFib, 1000);
        
        if( Src1_Code == 0 )
            strFib = sprintf( '' );
        else
            Src1_Alpha = char( int16('a')+floor((Src1_Code-1)/obj.Mtg_Data.Gen_Params.NbSrcPerBank) );
            Src1_Number = mod((Src1_Code-1),obj.Mtg_Data.Gen_Params.NbSrcPerBank)+1;
            strFib = sprintf( '%c%d', Src1_Alpha, Src1_Number );

            %2nd Src (001XXX-016XXX)
            if( obj.Mtg_Data.Gen_Params.Nb_Longueur_Onde > 1 )
                Src2_Code = floor(nFib/1000);
                Src2_Alpha = char( int16('a')+floor((Src2_Code-1)/obj.Mtg_Data.Gen_Params.NbSrcPerBank) );
                Src2_Number = mod((Src2_Code-1),obj.Mtg_Data.Gen_Params.NbSrcPerBank)+1;
                strFib = sprintf( '%s%c%d', strFib, Src2_Alpha, Src2_Number );
            end
        end
    end     
    %%%%%%%%FIN NOMENCLATURE ISS
        case 'IMAGINC'
            if( nFib >= 1000000 )        
                %Detectors
                Det_Code = nFib/1000000;
                if( Det_Code == 0 )
                    strFib = sprintf( '' );
                else 
                    strFib = sprintf( 'D%d', Det_Code );
                end
            else
            %Sources        
            %E1 a1b2 iss number 018001, 
            %E2 a3b4 iss number 020003
            %E3 a5b6 iss number 022005
            strFib = ISSIDnumber2Imaginc(nFib);
            end
        case 'NIRx'
            if( nFib >= 1000000 )        
                %Detectors
                Det_Code = nFib/1000000;
                if( Det_Code == 0 )
                    strFib = sprintf( '' );
                else 
                    strFib = sprintf( 'D%d', Det_Code );
                end
            else
            %Sources        
            %E1 a1b2 iss number 018001, 
            %E2 a3b4 iss number 020003
            %E3 a5b6 iss number 022005
            strFib = ISSIDnumber2Imaginc(nFib);
            end
        case 'IMAGINCHUBERT'
            if( nFib >= 1000000 )        
                %Detectors
                Det_Code = nFib/1000000;
                if( Det_Code == 0 )
                    strFib = sprintf( '' );
                else 
                    strFib = sprintf( 'D%d', Det_Code );
                end
            else
            %Sources        
            %E1 a1b2 iss number 001001
            strFib = ISSIDnumber2ImagincHubert(nFib);
            end
        otherwise
            strFib = 'ND'
    end
