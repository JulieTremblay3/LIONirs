function [Helm, bOk] = Set_Mtg_DefaultSrcFibersCfg( Helm )
            
    global PrjData;
    vHoles = get_vHoles(Helm);
    sMtg = get_Mtg(Helm);
    
    bOk = false;
    
    [PhysicalSrcGroupSyncNo,PhysicalSrcCombinations] = get_PhysicalSrcGroups( Helm );
    
    NbAttributedFibers = 0;
    NbFibersRemaining = numel(find(sMtg.v_pSrc));
    
%     if( ~numel( Helm.mat_ilHoles ) )
%         disp( 'DefaultSrcFibersCfg failed:  Empty. Using optimised algorithm instead' );
%         return;
%     end
%     
%     NbColsGD = size( Helm.mat_ilHoles, 2 );
%     NbLines = size( Helm.mat_ilHoles, 1 );
%     
%     %Attribution des fibres de src selon l'ordre correspondant a la map 2D
%     for( Col = 1:NbColsGD )
%         for( Lin = 1:NbLines )
% 
%             %Indice lineaire de passage entre la carte 2D et le vecteur
%             p = Helm.mat_ilHoles(Lin,Col);
% 
%             if( p && Helm.v_Holes(p).Type == 400 && numel(find(sMtg.v_pSrc == p)) )
%                 SrcCode = 0;
% 
%                 if( sMtg.Gen_Params.Nb_Longueur_Onde == 1 )
%                     SrcCode = PhysicalSrcCombinations(NbAttributedFibers+1,1);
%                 elseif( sMtg.Gen_Params.Nb_Longueur_Onde == 2 )
%                     SrcCode = PhysicalSrcCombinations(NbAttributedFibers+1,1)+ ...
%                               PhysicalSrcCombinations(NbAttributedFibers+1,2)*1000;
%                 end
% 
%                 if( SrcCode )
%                     %disp( SrcCode );
%                     sMtg.v_HolesMtg(p) = SrcCode;
%                     NbAttributedFibers = NbAttributedFibers+1;
%                     NbFibersRemaining = NbFibersRemaining-1;
%                 end
%             end
%         end
%     end
%     
    %S'il reste des sources a distribuer (si des trous ne faisaient pas
    % partie de la map 2D), il faut placer les fibres restantes dans ces 
    % trous.
    if( NbFibersRemaining )
        for( iSrc=1:numel(find(sMtg.v_pSrc)) )
            pSrc = sMtg.v_pSrc(iSrc);
            
            if( pSrc )
                if( ~sMtg.v_HolesMtg(pSrc) )
                    
                    SrcCode = 0;
                    if( sMtg.Gen_Params.Nb_Longueur_Onde == 1 )
                        SrcCode = PhysicalSrcCombinations(NbAttributedFibers+1,1);
                    elseif( sMtg.Gen_Params.Nb_Longueur_Onde == 2 )
                        SrcCode = PhysicalSrcCombinations(NbAttributedFibers+1,1)+ ...
                                  PhysicalSrcCombinations(NbAttributedFibers+1,2)*1000;
                    end
                    
                    if( SrcCode )
                        sMtg.v_HolesMtg(pSrc) = SrcCode;
                        NbAttributedFibers = NbAttributedFibers+1;
                        NbFibersRemaining = NbFibersRemaining-1;
                    end
                end
            end
        end
    end
    
    %Si toutes les sources ont ete associes a des fibres, tout est beau
    if( ~NbFibersRemaining )
        bOk = true;
    end
    
    Helm = set_Mtg(Helm, sMtg);