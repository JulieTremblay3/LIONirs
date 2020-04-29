%**************************************************************************
% FONCTION : HELMET/Label2p
%
% INTRANTS : obj          -> objet Helmet
%                         -> strHole
%
% EXTRANTS : Retourne l'indice lineraire a partir du label
%
%**************************************************************************
function p = Label2p( oHelmet, strHole )

    p = 0;
    
    for( pItem=1:numel(oHelmet.v_Holes) );
        
        nCar = length(oHelmet.v_Holes(pItem).Label);
        
        if( strncmp( strHole, oHelmet.v_Holes(pItem).Label, nCar ) )
            p = pItem;
        end
    end