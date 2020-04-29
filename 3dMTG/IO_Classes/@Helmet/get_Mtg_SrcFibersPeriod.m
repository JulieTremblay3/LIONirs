%**************************************************************************
% FONCTION : HELMET/get_Mtg_SrcFibersPeriod
%
% INTRANTS : obj     -> objet Helmet
%            SrcCode -> Code de source format 0XX0XX
%
% EXTRANTS : Periode d'illumination de source.
%            
% UTILITE : Retourne le numero de periode (1-NbPeriods). 
%           Retourne 0 si aucune source n'est attribuee au trou.
%**************************************************************************
function PeriodNo = get_Mtg_SrcFibersPeriod( obj, SrcCode )
    
    PeriodNo = 0;
    
    if( SrcCode ) 
        if( obj.Mtg_Data.Gen_Params.Nb_Longueur_Onde == 1 )
            PeriodNo = mod(SrcCode-1,obj.Mtg_Data.Gen_Params.NbSrcPerBank)+1;
        elseif( obj.Mtg_Data.Gen_Params.Nb_Longueur_Onde == 2 )
            UpperCode = floor(SrcCode/1000);
            LowerCode = mod(SrcCode,1000);
            
            UpperPeriodNo = floor(mod(UpperCode-1,obj.Mtg_Data.Gen_Params.NbSrcPerBank)/2)+1;
            LowerPeriodNo = floor(mod(LowerCode-1,obj.Mtg_Data.Gen_Params.NbSrcPerBank)/2)+1;
            
            %Verification, au cas ou une des deux sources serait
            %non-attribuee
            if( UpperCode )
                PeriodNo = UpperPeriodNo;
            elseif( LowerCode )
                PeriodNo = LowerPeriodNo;
            end
            
            %Verification: les deux periodes sont-elles les memes?
            if( UpperPeriodNo ~= LowerPeriodNo )
                PeriodNo = -1;
                disp( sprintf( 'Error in get_Mtg_SrcFibersPeriod(), SrcCode:%d, UpperPeriodNo:%d, LowerPeriodNo:%d', ...
                      SrcCode, UpperPeriodNo, LowerPeriodNo ) );
            end
        end
    end