function [Helm, bOk] = Set_Mtg_OptimisedSrcFibersCfg( Helm )
            
    vHoles = get_vHoles(Helm);
    sMtg = get_Mtg(Helm);
    
    bOk = false;
    
    %Fonction de grouppement de sources par paire (si 2 longueur d'ondes)
    %PhysicalSrcCombinations(:,1) = 830 nm
    %PhysicalSrcCombinations(:,2) = 690 nm
    [PhysicalSrcGroupSyncNo,PhysicalSrcCombinations] = get_PhysicalSrcGroups( Helm );

    %Fonction permettant de construire les structures d'information: Ces
    %structures contiennent l'information sur les interferences src-det a
    %proximite. Les 
    [SrcInfo, DetInfo, PeriodsDisp] = Built_SrcDetStructures( Helm );
    
    %Donne le nombre de periodes disponibles pour l'illumination des
    %sources (fibres optiques combines)
    NbPeriods = get_Mtg_NbSrcPeriods(Helm); %Nombre de Moments d'impulsions (périodes) libres
    
 	%**********************
    % ETAPE 2 - ALGO
    %**********************
    Finished = false;
    ImpossibleMtg = false;
    pContaminationDet = 0;
    while( ~Finished && ~ImpossibleMtg && numel(sMtg.v_pDet) )
        
        %******************************************************************
        % ETAPE 2.1: Recherche du détecteur soumis au plus grand nombre de
        %            contraintes. 
        %******************************************************************
        %Recherche du min{ NbAvailablePeriods-RemainingSrcs }
        iDet = 0;
        MinDiffFound = NbPeriods+1;
        for( i=1:numel(sMtg.v_pDet) )
            
            %S'il reste des sources a placer pour ce detecteur
            if( DetInfo(i).RemainingSrcs )
                
                %Evaluer le 'niveau' de possibilités
                ActualDiff = DetInfo(i).NbAvailablePeriods - DetInfo(i).RemainingSrcs;

                if( ActualDiff < MinDiffFound )

                    if( ActualDiff < 0 )
                        ImpossibleMtg = true;
                        pContaminationDet = sMtg.v_pDet(i);
                    end

                    MinDiffFound = ActualDiff;
                    iDet = i; %Most constrained detector - begin with this one
                end
            end
        end
        
        if( ~iDet )
            Finished = true;
            break;
        end
            

        %******************************************************************
        % ETAPE 2.2: Pour le détecteur iDet, recherche de la source ayant 
        %            la période possible la moins sollicitée et qui, parmi
        %            les sources ayant cette période, est celle qui possède
        %            le moins de périodes possibles. (est-ce correct ???)
        %******************************************************************  
        
        NbSrcInRange = length(DetInfo(iDet).v_iSrcInRange);
        
        %Recherche de la source possédant la période possible qui revient le moins souvent   
        PeriodsRepeat = zeros(1,NbPeriods);

        %Sommation des repetitions de periodes
        for( i=1:NbSrcInRange )
            iSrc = DetInfo(iDet).v_iSrcInRange(i);
            
            %Pour les sources qui n'ont toujours pas de periodes attribues
            if( ~SrcInfo(iSrc).Period )
                PeriodsRepeat = PeriodsRepeat+SrcInfo(iSrc).v_AvailablePeriods;
            end
        end
        
        %voir: find() nonzeros() pour eliminer la boucle
        %Recherche de la periode non-nulle revenant le moins souvent
        nMinPeriodRepeat = NbSrcInRange + 1;
        TargetPeriod = 0;
        for( i=1:NbPeriods )
            if( PeriodsRepeat(i) )
                if( PeriodsRepeat(i) < nMinPeriodRepeat )
                    nMinPeriodRepeat = PeriodsRepeat(i);
                    TargetPeriod = i;
                end
            end
        end
        
        if( ~TargetPeriod )
            ImpossibleMtg = true;
            break;
        end
        
        %Rechercher, parmi les source possédant la période recherchée,
        %celle dont les possibilités de périodes sont le plus
        %limitées.
        MinPossibilities = NbPeriods+1;
        iSrcFound = 0;
        for( i=1:NbSrcInRange )
            iSrc = DetInfo(iDet).v_iSrcInRange(i);
            
            %Si la source possède la période recherchée
            if( SrcInfo(iSrc).v_AvailablePeriods(TargetPeriod) )
                
                %Vérifier si cette source possède le minimum de
                %possibilités.
                if( length( find( SrcInfo(iSrc).v_AvailablePeriods > 0 ) ) < MinPossibilities )
                    MinPossibilities = length( find( SrcInfo(iSrc).v_AvailablePeriods > 0 ) );
                    iSrcFound = iSrc;
                end
            end
        end
        
        %******************************************************************
        % ETAPE 2.3: Nous connaissons maintenant le trou de source
        %            (iSrcFound) ainsi que la période d'impulsion cible
        %            (TargetPeriod). Il reste maintenant à attribuer 
        %            les numéros de sources physiques (fibres optiques) 
        %            dans ce trou et à mettre à jour la structure.
        %******************************************************************
        
        SrcInfo(iSrcFound).Period = TargetPeriod;
        SrcInfo(iSrcFound).v_AvailablePeriods = zeros(1,NbPeriods);
        SrcInfo(iSrcFound).v_ConflictingPeriods = ones(1,NbPeriods);
        
        for( i=1:length( SrcInfo(iSrcFound).v_iDetInRange ) )
            iCloseDet = SrcInfo(iSrcFound).v_iDetInRange(i);
            
            DetInfo(iCloseDet).v_UsedPeriods(TargetPeriod) = true;
            DetInfo(iCloseDet).NbAvailablePeriods = DetInfo(iCloseDet).NbAvailablePeriods-1;
            DetInfo(iCloseDet).RemainingSrcs = DetInfo(iCloseDet).RemainingSrcs-1;
        end
        
        
        %*** A AJOUTER : LORSQUE LES BANQUES SONT EPUISEES, RETIRER LES
        %                PERIODES QUI NE SONT PLUS DISPONIBLES. ->Fait
        %Soit:
        PeriodsDisp(TargetPeriod) = PeriodsDisp(TargetPeriod)-1;
        
        %Si une periode s'est completement videe
        if( ~PeriodsDisp(TargetPeriod) )
            
            for( i=1:numel(sMtg.v_pDet) )

                if( isempty( find( find( DetInfo(i).v_UsedPeriods ) == TargetPeriod ) ) )
                    DetInfo(i).v_UsedPeriods(TargetPeriod) = true;
                    DetInfo(i).NbAvailablePeriods = DetInfo(i).NbAvailablePeriods-1;
                end
            end
        end
        
        %Mise a jour des sources qui n'ont pas encore de periode attribuee
        for( i=1:numel(sMtg.v_pSrc) )
            
            if( ~SrcInfo(i).Period )
                NbDetInRange = length( SrcInfo(i).v_iDetInRange );

                for( j=1:NbDetInRange )
                    iCloseDet = SrcInfo(i).v_iDetInRange(j);
                    SrcInfo(i).v_ConflictingPeriods = SrcInfo(i).v_ConflictingPeriods | DetInfo(iCloseDet).v_UsedPeriods | ( PeriodsDisp == 0 );
                end

                SrcInfo(i).v_AvailablePeriods = ~SrcInfo(i).v_ConflictingPeriods;
            end
        end
    end
    
    %******************************************************************
    % ETAPE 3: Attribution des sources physiques 
    %******************************************************************
    
    if( ImpossibleMtg )
        disp( sprintf( 'Contamination Detected: cannot have more \n than %d sources in range of detectors.', NbPeriods) );
        %h = warndlg( sprintf( 'Contamination Detected: cannot have more \n than %d sources in range of detectors.', NbPeriods) );
        %h = errordlg( sprintf( 'Contamination Detected: cannot have more \n than %d sources in range of detectors.', NbPeriods) ,'Warning', 'on');
        %uiwait(h);
        return;
    end  
    
	PeriodUsage = zeros( NbPeriods );
    StrDebugPer = '';
    StrDebug_v_pSrc = 'v_pSrc = [';
    
    %Boucle d'attribution de sources qui sont "in detector range"
    for( i=1:numel(sMtg.v_pSrc) )
        p = sMtg.v_pSrc(i);

        %Sauter les banques dont la periode est utilisee
        %(nPeriod = Period + NbPeriode*PeriodesUtilisees)
        if( SrcInfo(i).Period )
            CombiNo = SrcInfo(i).Period+PeriodUsage(SrcInfo(i).Period)*NbPeriods;
            sMtg.v_HolesMtg(p) = PhysicalSrcCombinations(CombiNo,1) + PhysicalSrcCombinations(CombiNo,2)*1000;
            PeriodUsage(SrcInfo(i).Period) = PeriodUsage(SrcInfo(i).Period)+1;
        end
    end
    
    %Boucle d'attribution de sources qui sont "out of detector range"
    for( i=1:numel(sMtg.v_pSrc) )
        p = sMtg.v_pSrc(i);

        %Sauter les banques dont la periode est utilisee
        if( ~SrcInfo(i).Period )
            CombiNo = 0;
            TargetPeriod = 0;
            
            %Recherche d'une periode libre (N'importe quelle)
            for( iSeekPeriod=1:NbPeriods ) 
                if( PeriodUsage(iSeekPeriod) < 4)%get_Mtg_PeriodMultiplicity(Helm) )
                    TargetPeriod = iSeekPeriod;
                    break;
                end
            end
            
            %Il ne reste plus de fibres disponibles.
            if( ~TargetPeriod )
                disp( 'Cannot add more sources: Max amount reached' );
                %h = errordlg( sprintf( 'Cannot add more sources: Max amount reached' ) ,'Warning' );
                %uiwait(h);
                return;
            end
            
            CombiNo = TargetPeriod+PeriodUsage(TargetPeriod)*NbPeriods;
            sMtg.v_HolesMtg(p) = PhysicalSrcCombinations(CombiNo,1) + PhysicalSrcCombinations(CombiNo,2)*1000;
            PeriodUsage(TargetPeriod) = PeriodUsage(TargetPeriod)+1;
        end
    end
        disp( StrDebugPer );
        %disp( sprintf( '%s ]', StrDebug_v_pSrc ) );
 

    Helm = set_Mtg(Helm, sMtg);
    bOk = true;
    
%--------------------------------------------------------------------------
% fonction Built_SrcDetStructures
%
% Extrants:
%
% DetInfo(1..NbDet).v_iSrcInRange(0..?)
%                  .v_UsedPeriods(1..NbPeriods)
%                  .NbAvailablePeriods
%                  .RemainingSrcs
%
% SrcInfo(1..NbSrc).v_iDetInRange(0..?)
%                  .v_ConflictingPeriods (1..NbPeriods)
%                  .v_AvailablePeriods(1..NbPeriods)
%                  .Period
% 
% PeriodsDisp(1..NbPeriods) : booleens de disponibilite
%
%--------------------------------------------------------------------------
function [SrcInfo, DetInfo, PeriodsDisp] = Built_SrcDetStructures( Helm )
        
    SrcInfo = [];
    DetInfo = [];
    PeriodsDisp = [];
    
    sMtg = get_Mtg( Helm );
    DistCont = sMtg.Gen_Params.DistContamination;
    NbPeriods = get_Mtg_NbSrcPeriods(Helm);
    
    %Preparation du tableau des periodes 
    for( NoPeriod=1:NbPeriods )
        PeriodsDisp(NoPeriod) = get_Mtg_PeriodMultiplicity(Helm);
    end  
    
    for( i=1:numel(sMtg.v_pDet) )
        %Initialisation de la structure de Detecteurs
        DetInfo(i).v_iSrcInRange = [];  %Liste des sources à proximité de ce détecteur
        DetInfo(i).v_UsedPeriods = zeros(1,NbPeriods);  %Moments d'impulsions (périodes) utilisées
        DetInfo(i).NbAvailablePeriods = NbPeriods; %Nombre de Moments d'impulsions (périodes) libres
        DetInfo(i).RemainingSrcs = 0;   %Nombre de sources restant a placer
    end
    
    %Initialisation du graphe de proximite des sources et des detecteurs.
    for( i=1:numel(sMtg.v_pSrc) )
        
        SrcInfo(i).v_iDetInRange = []; %Liste des détecteurs à proximité de cette source
        SrcInfo(i).v_ConflictingPeriods = zeros(1,NbPeriods); %Union des moments d'impulsions utilisés par les détecteurs à proximité.
        SrcInfo(i).v_AvailablePeriods = ones(1,NbPeriods); %Liste des périodes disponibles (qui ne sont pas conflictuelles).
        SrcInfo(i).Period = 0;         %Moment d'impulsion attribué (période)
        
        pSrc = sMtg.v_pSrc(i);
        
        for( j=1:numel(sMtg.v_pDet) )
            
            pDet = sMtg.v_pDet(j);
            
            %Verifier si la source peut etre percue par le detecteur
            if( get_DistHoles( Helm, pSrc, pDet ) < DistCont )
                
             	NbDetsInRange = length(SrcInfo(i).v_iDetInRange)+1;
                
                SrcInfo(i).v_iDetInRange(NbDetsInRange) = j;
                
                SourcesAroundThisDet = length(DetInfo(j).v_iSrcInRange)+1;
                
                %Ajout de la source dans la liste de src a proximite du det
                DetInfo(j).v_iSrcInRange(SourcesAroundThisDet) = i;

                %Mise a jour du nombre de sources a placer pour ce
                %detecteur.
                DetInfo(j).RemainingSrcs = SourcesAroundThisDet;
            end
        end
    end
    
    
    