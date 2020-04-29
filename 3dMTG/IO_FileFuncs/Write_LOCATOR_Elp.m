function Write_LOCATOR_Elp( oDig_CompleteHelmet, PathFileName )

    if( ~PathFileName ) 
        return;
    end
    
    v_Probes = get_vProbes( oDig_CompleteHelmet );
    matFiducials = get_matFiducials( oDig_CompleteHelmet );
    
    fid_elp = fopen(PathFileName, 'w');
    fprintf(fid_elp, '3\t2\n');
    fprintf(fid_elp, '//Probe file\n');
    fprintf(fid_elp, '//Minor revision number\n');
    fprintf(fid_elp, '1\n');
    fprintf(fid_elp, '//ProbeName\n');
    fprintf(fid_elp, '%%%%N\tName\n');
    fprintf(fid_elp, '//Probe type, number of sensors\n');
    fprintf(fid_elp, '0\t%d\n', numel(v_Probes));
    fprintf(fid_elp, '//Position of fiducials X+, Y+, Y- on the subject\n');
    fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', 0, 0, 0 );
    fprintf(fid_elp, '//\tNZ\t%.6f\t%.6f\t%.6f\n', matFiducials(1,1)*1000, matFiducials(1,2)*1000, matFiducials(1,3)*1000 );  %m to mm (*1000)
    fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', 0, 0, 0 );
    fprintf(fid_elp, '//\tLPA\t%.6f\t%.6f\t%.6f\n', matFiducials(2,1)*1000, matFiducials(2,2)*1000, matFiducials(2,3)*1000 ); %m to mm (*1000)
    fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', 0, 0, 0 );
    fprintf(fid_elp, '//\tRPA\t%.6f\t%.6f\t%.6f\n', matFiducials(3,1)*1000, matFiducials(3,2)*1000, matFiducials(3,3)*1000 ); %m to mm (*1000)
    
    Counter = 0;
    for( Pos=1:numel(v_Probes) )
        Counter = Counter+1;
        fprintf(fid_elp, '//Sensor type\n');
        fprintf(fid_elp, '%%S\t%d\n', v_Probes(Pos).Type );
        fprintf(fid_elp, '//Sensor name and data for sensor# %d\n', Counter );
        fprintf(fid_elp, '%%N\t%s\n', v_Probes(Pos).Label );
        fprintf(fid_elp, '%.6f\t%.6f\t%.6f\n', v_Probes(Pos).Coord(1), ...
                                               v_Probes(Pos).Coord(2), ...
                                               v_Probes(Pos).Coord(3) );
                                           
        fprintf(fid_elp, '//\t%s\t%.6f\t%.6f\t%.6f\n', v_Probes(Pos).Label, ...
                                                       v_Probes(Pos).Coord(1)*1000, ...%m to mm (*1000)
                                                       v_Probes(Pos).Coord(2)*1000, ...%m to mm (*1000)
                                                       v_Probes(Pos).Coord(3)*1000 );  %m to mm (*1000)
    end

    fclose(fid_elp);
