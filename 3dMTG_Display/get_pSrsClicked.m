function pClosest = get_pSrsClicked( Helm, CurrentLine,srslst)   
%Trouve le point de la source clique la plus proche dans la liste srslst
    %Droite du click de souris en 3D
    mPt = CurrentLine;

    for ind = 1:length(srslst)
        srslst(ind) = SDpairs2Srs_n(srslst(ind),1);
    end
    
    %**************************************************************************
    % Recherche la source ou le détecteur le plus proche du vecteur de click :
    %**************************************************************************
    sMtg = get_Mtg( Helm);
    vHoles = get_vHoles( Helm );

    Closest = 10000;
    pClosest = 0;
    
    for( ind_s = 1:length(sMtg.v_pSrc))
        p = sMtg.v_pSrc(ind_s);
        srs_montage = sMtg.v_HolesMtg(sMtg.v_pSrc(ind_s));
        if srs_montage>0
         
        if ~isempty(find(srslst==srs_montage))
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
    end