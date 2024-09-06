function oDig = set_vProbes( oDig, vProbes )

    if( isfield( vProbes, 'Coord' ) && isfield( vProbes, 'Label' ) && isfield( vProbes, 'Type' ) )
        oDig.vProbes = vProbes;
    else
        disp( 'function set_vProbes : Error loading digitization probes' );
    end
    
