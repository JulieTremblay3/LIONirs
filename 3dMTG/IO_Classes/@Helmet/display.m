function display(obj)
%HELMET/DISPLAY command prompt call
disp(' ');
disp([inputname(1),' = ']);
disp( ' ' );
disp( 'HELMET:' );
disp( ' ' );
disp( sprintf('            ID = %d ', obj.ID ) );
disp( sprintf('       v_Holes = %d struct', size(obj.v_Holes,2)  ) );
disp( sprintf('     Mtg:NbDet = %d', numel( find(obj.Mtg_Data.v_pDet) ) ) );  
disp( obj.Mtg_Data.v_pDet );
disp( ' ' );