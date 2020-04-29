function Write_Mtg_Pos( oHelmet, PathName, FileName )

    %Attention: Conversion en m!
    
    v_Holes = get_vHoles(oHelmet);
    sMtg = get_Mtg(oHelmet);
    
    if( ~FileName ) 
        return;
    end
    
    pItems = [sMtg.v_pSrc,sMtg.v_pDet];
    nItems = numel(pItems);
    
    fid_pos = fopen([PathName,FileName], 'w');
    
    fprintf(fid_pos, '[FileInfo]\n');
    fprintf(fid_pos, 'ID=3DHeadShape\n');
    fprintf(fid_pos, 'VER=1.20\n');
    fprintf(fid_pos, 'MODE=1\n');
    fprintf(fid_pos, 'TYPE=0\n');
    fprintf(fid_pos, 'ProbeNum=%d\n', nItems);
    fprintf(fid_pos, '[User]\n');
    fprintf(fid_pos, 'Name=John Doe\n');
    fprintf(fid_pos, 'ID=0001\n');
    fprintf(fid_pos, 'Comment=None\n');
    fprintf(fid_pos, 'Sex=Unknown\n');
    fprintf(fid_pos, '[LeftEar]\n');
    fprintf(fid_pos, 'X=%.2f\n', sMtg.matFiducials(2,1)*1000 );
    fprintf(fid_pos, 'Y=%.2f\n', sMtg.matFiducials(2,2)*1000 );
    fprintf(fid_pos, 'Z=%.2f\n', sMtg.matFiducials(2,3)*1000 );
    fprintf(fid_pos, '[RightEar]\n');
    fprintf(fid_pos, 'X=%.2f\n', sMtg.matFiducials(3,1)*1000 );
    fprintf(fid_pos, 'Y=%.2f\n', sMtg.matFiducials(3,2)*1000 );
    fprintf(fid_pos, 'Z=%.2f\n', sMtg.matFiducials(3,3)*1000 );
    fprintf(fid_pos, '[Nasion]\n');
    fprintf(fid_pos, 'X=%.2f\n', sMtg.matFiducials(1,1)*1000 );
    fprintf(fid_pos, 'Y=%.2f\n', sMtg.matFiducials(1,2)*1000 );
    fprintf(fid_pos, 'Z=%.2f\n', sMtg.matFiducials(1,3)*1000 );
%     fprintf(fid_pos, '[Inion]\n');
%     fprintf(fid_pos, 'X=%.2f\n', 0 );
%     fprintf(fid_pos, 'Y=%.2f\n', 0 );
%     fprintf(fid_pos, 'Z=%.2f\n', 0 );
%     fprintf(fid_pos, '[Top]\n');
%     fprintf(fid_pos, 'X=%.2f\n', 0 );
%     fprintf(fid_pos, 'Y=%.2f\n', 0 );
%     fprintf(fid_pos, 'Z=%.2f\n', 0 );
    
    for( iPos=1:nItems )
        p = pItems(iPos);
        fprintf(fid_pos, '[Probe%s]\n', cell2mat(strread(get_HoleFiberID(oHelmet,p), '%s', 1)) );
        fprintf(fid_pos, 'X=%.2f\n', v_Holes(p).Coord.x*1000 );
        fprintf(fid_pos, 'Y=%.2f\n', v_Holes(p).Coord.y*1000 );
        fprintf(fid_pos, 'Z=%.2f\n', v_Holes(p).Coord.z*1000 );
%         fprintf(fid_pos, 'A=%.2f\n', 0 );
%         fprintf(fid_pos, 'E=%.2f\n', 0 );
%         fprintf(fid_pos, 'R=%.2f\n', 0 );
        fprintf(fid_pos, 'NX=%.5f\n', v_Holes(p).Normal.x );
        fprintf(fid_pos, 'NY=%.5f\n', v_Holes(p).Normal.y );
        fprintf(fid_pos, 'NZ=%.5f\n', v_Holes(p).Normal.z );
    end
    
%     fprintf(fid_pos, '[Angle]\n');
%     fprintf(fid_pos, 'X=%.2f\n', 0 );
%     fprintf(fid_pos, 'Y=%.2f\n', 0 );
%     fprintf(fid_pos, 'Z=%.2f\n', 0 );

    fclose(fid_pos);
    