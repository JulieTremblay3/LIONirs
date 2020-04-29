function [Helm, bOk] = SetDetectorsFiber_Withoutoptimising(Helm,pDet)
 
sMtg = get_Mtg( Helm );          
i = numel(sMtg.v_pDet);   
sMtg.v_HolesMtg(pDet) = i*1000000;
bOk = true;
Helm = set_Mtg(Helm, sMtg);