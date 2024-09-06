%**************************************************************************
% FONCTION : HELMET/get_DistHoles
%
% INTRANTS : obj          -> objet Helmet
%                         -> Indice lineaire du 1er trou
%                         -> Indice lineaire du 2e  trou
%
% EXTRANTS : Distance separant les deux trous specifies. 
%
%**************************************************************************
function dist = get_DistHoles( obj, pH1, pH2 )
    CoordH1 = [ obj.v_Holes(pH1).Coord.x, 
                obj.v_Holes(pH1).Coord.y, 
                obj.v_Holes(pH1).Coord.z ];
    CoordH2 = [ obj.v_Holes(pH2).Coord.x, 
                obj.v_Holes(pH2).Coord.y, 
                obj.v_Holes(pH2).Coord.z ];
    NormH1 =  [ obj.v_Holes(pH1).Normal.x, 
                obj.v_Holes(pH1).Normal.y, 
                obj.v_Holes(pH1).Normal.z ];
    NormH2 =  [ obj.v_Holes(pH2).Normal.x, 
                obj.v_Holes(pH2).Normal.y, 
                obj.v_Holes(pH2).Normal.z ];
    P1 = CoordH1-NormH1*obj.v_Holes(pH1).SkinDepth;
    P2 = CoordH2-NormH2*obj.v_Holes(pH2).SkinDepth;
         
    P1P2 = P2-P1;
	dist = (P1P2'*P1P2)^0.5;