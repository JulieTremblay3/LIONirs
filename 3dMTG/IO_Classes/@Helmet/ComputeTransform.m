function oHelmet = ComputeTransform( oHelmet )

    vHoles = oHelmet.v_Holes;
    
    for( p = 1:numel(vHoles) )
        
        vn = [ vHoles(p).Normal.x; vHoles(p).Normal.y; vHoles(p).Normal.z ];
        vz = [ 0; 0; 1; ];

        %Ajustement de l'orientation
        T = AlignementVecteurs( vz, vn );
        
        %Ajustement de la position
        T(4,1) = vHoles(p).Coord.x;
        T(4,2) = vHoles(p).Coord.y;
        T(4,3) = vHoles(p).Coord.z;
        
        vHoles(p).Transformation = T;
    end
    
    oHelmet = set_vHoles( oHelmet, vHoles );
    