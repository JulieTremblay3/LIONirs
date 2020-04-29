function [oHelmet, oDigSource, oDigDestination] = TransferDigitizationItem( oHelmet, oDigSource, oDigDestination, v_pItemsFromSource );
    
    vHoles = get_vHoles( oHelmet );
    
    %Cellules de Labels des selections
    vCellLabels = {vHoles(v_pItemsFromSource).Label};
    
    %Booleens de provenance
    v_bFromCompleteHelmet = [vHoles(v_pItemsFromSource).IsFromCompleteHelmet];
    
    %Verifier que la provenance des elements est homogene
    if(   ( numel(find(v_bFromCompleteHelmet)) ~= numel(v_bFromCompleteHelmet) )...
        && ~isempty(find(v_bFromCompleteHelmet)) )
        warndlg_modal( 'Invalid selection: different items provenance' );
        return;
    end
    
    %Flag: quel est la source?
    if( find(v_bFromCompleteHelmet) )
        bSourceIsCompleteHelmet = true;
    else
        bSourceIsCompleteHelmet = false;
    end
    
    [oDigRegisteredSrc, bSrcX1X2X3Found, bDstX1X2X3Found] = RegisterDigitization( oDigSource, oDigDestination );
    
    vSrcProbes    = get_vProbes( oDigSource );
    vRegSrcProbes = get_vProbes( oDigRegisteredSrc );
    vDstProbes    = get_vProbes( oDigDestination );
    
    v_pSrcProbesTransfered = zeros(1,numel(vCellLabels));
    
    %Pour chaque Label d'item selectionne
    for( iElem = 1:numel(vCellLabels) )
        
        pSrcProbe = find( strcmp({vSrcProbes.Label}, vCellLabels{iElem} ) );
        pDstProbe = find( strcmp({vDstProbes.Label}, vCellLabels{iElem} ) );
        
        if( numel(pSrcProbe)==1 && pSrcProbe )
            
            %Si l'element existe deja dans la destination
            if( numel(pDstProbe)==1 && pDstProbe )
                %Remplacer l'element
                vDstProbes(pDstProbe) = vRegSrcProbes(pSrcProbe);                
            else
                %Ajouter l'element
                vDstProbes(numel(vDstProbes)+1) = vRegSrcProbes(pSrcProbe);
            end
            
            v_pSrcProbesTransfered(iElem) = pSrcProbe;
        end
        
        %Ne rien faire: plusieurs element ont le meme Label, ce qui ne
        %devrait pas arriver.
        if( numel(pSrcProbe)>1 || numel(pDstProbe)>1 )
            disp( 'TransferDigitizationItem::Cannot work with Digitizations containing multiple elements with the same Label' );
        end
    end
    
    %Effacer les elements transferes
    v_pSrcProbesRemaining = [1:numel(vSrcProbes)];
    v_pSrcProbesRemaining( v_pSrcProbesTransfered ) = 0;
    vSrcProbes = vSrcProbes( find(v_pSrcProbesRemaining) );
    
    %Enregistrer les changements dans les digitalisations.
    oDigSource = set_vProbes( oDigSource, vSrcProbes );
    oDigDestination = set_vProbes( oDigDestination, vDstProbes );
    
    %Reconstruire un casque: effectuer la separation appropriee et
    %recalculer les normales.
    if( bSourceIsCompleteHelmet )
        oHelmet = FitHelmetOnSubject( oDigDestination, oDigSource );
    else
        oHelmet = FitHelmetOnSubject( oDigSource, oDigDestination );
    end
    
    