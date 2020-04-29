function vMarkers = Read_iMagic_Mrk( FullFileName )
    
    nMarkers = 0;
    vMarkers = [];

    %Ouverture du fichier .esp
    fid_MRK = fopen(FullFileName);
    if (fid_MRK < 2)
        msgbox( 'File could ne be found' );
        return;
    end
    
    %Lecture de la 1ere ligne de coordonnees
    a_line = fgets(fid_MRK);
    
    %Boucle de lecture de lignes du fichier
    while( isempty( find(a_line == -1) ) );

        stuff = strread(a_line, '%s' );
        
        if( strcmp( stuff{1}, 'Point' ) )
            %disp('Point found');
            nMarkers = nMarkers+1;
            a_line=fgets(fid_MRK);
            
            if( isempty( find(a_line == -1) ) )
                cellLabel = strread(a_line, '%s' );
                strLabel = cellLabel{1};
                %disp( sprintf( '   Name:%s', strLabel ) );
                vMarkers(nMarkers).Label = strLabel;
                a_line=fgets(fid_MRK);
            
                if( isempty( ~find(a_line == -1) ) )
                    nElems = strread(a_line);
                    %disp( sprintf( '      NumElems:%d', nElems ) );
                    a_line=fgets(fid_MRK);
                    
                    if( isempty( ~find(a_line == -1))  && (nElems > 0) )
                        Coord = strread( a_line, '%d' );
                        vMarkers(nMarkers).Coord = [ Coord(1), Coord(2), Coord(3), 1 ];
                    end
                end
            end
        end
        %Lecture de la ligne suivante
        a_line=fgets(fid_MRK);
    end
    

    fclose(fid_MRK);
end
    
    
