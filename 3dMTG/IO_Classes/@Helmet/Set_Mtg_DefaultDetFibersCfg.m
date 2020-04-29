function [Helm, bOk] = Set_Mtg_DefaultDetFibersCfg( Helm )
            
    global PrjData;
    vHoles = get_vHoles(Helm);
    sMtg = get_Mtg(Helm);
    bOk = false;
    
    if( numel(sMtg.v_pDet) > sMtg.Gen_Params.Nb_Detect_Disp )
        return;
    end
    
    for( iDet=1:numel(sMtg.v_pDet) )
        pDet = sMtg.v_pDet( iDet );
        sMtg.v_HolesMtg( pDet ) = iDet*1000000;
    end
    
    Helm = set_Mtg(Helm, sMtg);
    bOk = true;