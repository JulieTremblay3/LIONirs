function pClosest = get_pItemClicked( oInterDlgComm, oHelmet, CurrentLine )


    %Structure d'information sur le graphique du casque
    sHelmetAxeDisp = get_HelmetAxeDisp(oInterDlgComm);

    %Droite du click de souris en 3D
    mPt = CurrentLine;

    %**************************************************************************
    % Recherche de l'element le plus proche du vecteur de click :
    %**************************************************************************

    vHoles = get_vHoles( oHelmet );
    Closest = 10000;
    pClosest = 0;

    for( p=1:numel(vHoles))
         if( vHoles(p).Type == 400 )
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

            if( CurDist < Closest && sHelmetAxeDisp.v_bVisible(p) )
                Closest = CurDist;
                pClosest = p;
            end
        end
    end
    