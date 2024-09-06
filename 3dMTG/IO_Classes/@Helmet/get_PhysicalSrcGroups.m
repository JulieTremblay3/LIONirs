function [PhysicalSrcGroupSyncNo,PhysicalSrcCombinations] = get_PhysicalSrcGroups( oHelmet )
%**************************************************************************
% Regrouppement optimal de sources selon le nombre de longueurs d'ondes
% par trou d'emission.
%
%**************************************************************************
    sMtg = get_Mtg(oHelmet);
    
    %**********************************************************************
    % Attribution automatique optimisé de sources (multiplexées)
    %**********************************************************************
    
    %Nombre de trou maximal pouvant recevoir des sources (depend du nombre
    %de sources par trou)
    NbBanks = sMtg.Gen_Params.NbBanks;
    NbSrcPerBank = sMtg.Gen_Params.NbSrcPerBank;
    NbLambda = sMtg.Gen_Params.Nb_Longueur_Onde;
    MaxSrcHoles = floor((NbBanks*NbSrcPerBank)/NbLambda);
    
    %Liste de combinaisons de sources
    PhysicalSrcCombinations = zeros( MaxSrcHoles, NbLambda );
    PhysicalSrcGroupSyncNo = zeros( MaxSrcHoles );
    
	%Liste de sources physiques utilisees
    UsedPhysicalSrc = zeros( NbBanks, NbSrcPerBank);
    
    %Pour un nombre de src/bnk impair, ne pas utiliser les dernieres srcs
    MaxUsableSrcPerBank = NbSrcPerBank-mod(NbSrcPerBank,2);
    
    %Boucle generatrice de combinaisons de sources
    LastSrcNo = 0;
    SyncGroup = 1;
    for( h=1:MaxSrcHoles )
        
        LastBank = 0;
        HoleSrcNo = [];
        
        for( l=1:NbLambda )
            
            ValidSrcFound = 0;
            
            b = LastBank+1;
            while( b<=NbBanks && ~ValidSrcFound )
                
                %Regrouppement par paquets synchronisés:
                s = 1;
                
                %Groupe impair, ex: a1b2, a3b4, ... (# 2è src suit # 1ère)
                if( mod(SyncGroup,2) || NbLambda == 1 )
                    if( LastSrcNo == MaxUsableSrcPerBank )
                        SyncGroup = SyncGroup+1;
                        LastSrcNo = 0;
                    else
                        s = LastSrcNo+1;
                    end
                %Groupe pair, ex: a2b1, a4b3
                elseif( l == 1 && LastSrcNo == MaxUsableSrcPerBank-NbLambda+1 )
                    SyncGroup = SyncGroup+1;
                    LastSrcNo = 0;
                end

                while( s<=MaxUsableSrcPerBank && ~ValidSrcFound )
                    if( ~UsedPhysicalSrc(b,s) )
                        
                        %Verifier qu'aucune source ne s'allume en meme
                        %temps pour un meme trou
                        SyncOK = true;
                        for( ls=1:length(HoleSrcNo) )
                            if( HoleSrcNo(ls) == s )
                                SyncOK = false;
                                break;
                            end
                        end
                        
                        %Source disp et synchronisation temporelle ok
                        if( SyncOK )
                            UsedPhysicalSrc(b,s) = true;
                            ValidSrcFound = (b-1)*NbSrcPerBank+s;
                            LastBank = b;
                            LastSrcNo = s;
                            HoleSrcNo(length(HoleSrcNo)+1) = s;
                        end
                    end
                    s = s+1;
                end
                b = b+1;
            end
            
            PhysicalSrcCombinations(h,l) = ValidSrcFound;
            PhysicalSrcGroupSyncNo(h) = SyncGroup;
            
        end
    end