function pClosest = get_SrsDetClicked( Helm, CurrentLine, camerapos )
% Descrition : Choisit sur le casque la source ou le detecteur le plus pres
% Entrée : Helm : projet casque
% Current line : position du curseur
% Camerapos: Position de la caméra
% Si la position de la camera est à gauche on prend les sources à gauche 
%**************************************************************************
    mPt = CurrentLine;
    sMtg = get_Mtg( Helm);
    vHoles = get_vHoles( Helm );
    Closest = 10000;
    pClosest = 0;
    vSrsDet = [sMtg.v_pDet,sMtg.v_pSrc];    
   CurDisttot=  [];
   CampDisttot = [];
    for( ind_d= 1:length( vSrsDet))
        p = vSrsDet(ind_d);
            %1)Projection orthogonale
            %
            %               ->          ->
            %               proj        v1 
            %mPt(1,:)Eye--------o--------------------mPt(2,:)MouseTgt 
            %           \       :
            %            \      :
            %          -> \     :
            %          v2  \    :
            %               \   :
            %                \  :
            %                 \ :
            %                 Elem

            
         
        
                v1 = mPt(1,:)-mPt(2,:);
                v2 = [vHoles(p).Coord.x-mPt(1,1), vHoles(p).Coord.y-mPt(1,2), vHoles(p).Coord.z-mPt(1,3)];
                proj = ((v2*v1')/(v1*v1'))*v1;
                vDist = v2-proj;
                CurDist =  (vDist*vDist')^0.5;
                %find_distance( vDist(1), vDist(2), vDist(3), 0, 0, 0 );
                CurDisttot = [CurDisttot, CurDist];
                
                %Distance de la camera
                Campos = [vHoles(p).Coord.x-camerapos(1), vHoles(p).Coord.y-camerapos(2),vHoles(p).Coord.z-camerapos(3)];
                CampDist = (Campos*Campos')^0.5;
                CampDisttot = [CampDisttot, CampDist];
                if( CurDist < Closest )
                    Closest = CurDist;
                    pClosest = p;
                end
    end       
distance = CampDisttot - mean(CampDisttot);
for ind = 1:numel(CurDisttot) 
   [C,I] = min(CurDisttot);
    if  distance(I)< 0
     pClosest = vSrsDet(I);
     break
    else
      CurDisttot(I) = max(CurDisttot);
    end        
end
