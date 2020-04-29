function [Helm, bOk] = SetSourcesFiber_Withoutoptimising(Helm,pSrc)
 
[PhysicalSrcGroupSyncNo,PhysicalSrcCombinations] = get_PhysicalSrcGroups( Helm );
sMtg = get_Mtg( Helm );          
i = numel(sMtg.v_pSrc);   
sMtg.v_HolesMtg(pSrc) = PhysicalSrcCombinations(i,1) + PhysicalSrcCombinations(i,2)*1000;
bOk = true;
Helm = set_Mtg(Helm, sMtg);