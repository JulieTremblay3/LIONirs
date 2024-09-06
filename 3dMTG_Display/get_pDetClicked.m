function pClosest = get_pDetClicked( Helm, CurrentLine,detlst )
%Trouve le point du detecteur clique le plus proche dans la liste detlst
   %Droite du click de souris en 3D
    mPt = CurrentLine;

    %**************************************************************************
    % Recherche la source ou le détecteur le plus proche du vecteur de click :
    %**************************************************************************
    sMtg = get_Mtg( Helm);
    vHoles = get_vHoles( Helm );

    Closest = 10000;
    pClosest = 0;
   
    for( ind_d= 1:length( sMtg.v_pDet))
        det_montage = sMtg.v_HolesMtg(sMtg.v_pDet(ind_d))/1000000 ;
        if ~isempty(find(detlst==det_montage))
        p = sMtg.v_pDet(ind_d);
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

                if( CurDist < Closest )
                    Closest = CurDist;
                    pClosest = p;
                end
        end
    end
    1