function bContaminated = get_MtgSrc_Contaminated( Helm )
    
    %disp( 'get_MtgSrc_Contaminated()' );
    
    sMtg = get_Mtg( Helm );
    MaxPeriods = get_Mtg_NbSrcPeriods( Helm );
    
    bContaminated = true;
    
    %Verifier pour tous les detecteurs
    %disp( 'for( iDet=1:sMtg.NbDet )' );
    for( iDet=1:numel(sMtg.v_pDet) )
        pDet = sMtg.v_pDet(iDet);
        %disp( sprintf( '    iDet=%d (%s)', iDet, get_HoleFiberID( Helm, pDet ) ) );
        
        v_pSrcInRange = get_MtgDet_SrcInRangeOfContamination( Helm, pDet );
        
        if( numel(v_pSrcInRange) > MaxPeriods )
            return;
        end
        
        PeriodsUsed = zeros( MaxPeriods );
    
        %Verifier pour toutes les source en vue
        %disp( '    for( iSrc=1:numel(v_pSrcInRange) )' );
        for( iSrc=1:numel(v_pSrcInRange) )
            pSrc = v_pSrcInRange(iSrc);
            %disp( sprintf( '        iSrc=%d (%s)', iSrc, get_HoleFiberID( Helm, pSrc ) ) );
            
            Period = get_Mtg_SrcHolePeriod( Helm, pSrc );
            %disp( sprintf( '        %s : Period:%d', get_HoleFiberID(
            %Helm, pSrc ), Period ) );

            %Source toujours non-attribuee
            if( ~Period )
                %(DO NOTHING)
                
            %Erreur inconnue: Periode hors-bornes
            elseif( Period < 1 || Period > MaxPeriods )
                disp( sprintf( 'VerifiyContamination() pSrc:%d, pDet:%d, Period:%d', ...
                               pSrc, pDet, Period ) );
                return;
                
            %Contamination si la periode d'illumination est occupee
            elseif( PeriodsUsed(Period) )
                return;
            
            %Marquer la periode (1ere fois)
            else
                PeriodsUsed(Period) = true;
            end
        end
    end
    
    %Rendu ici, on s'est assure que le montage ne comporte pas de
    %contamination
    bContaminated = false;