function Display( oDig )
    disp(' ');
    disp([inputname(1),' = ']);
    disp( ' ' );
    disp( 'DIGITIZATION:' );
    disp( ' ' );
    disp( sprintf('       vProbes = struct containing %d element(s)', numel(oDig.vProbes) ) );
    disp( ' ' );