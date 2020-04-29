function Write_Brainsight_SEQ( oHelmet, PathFileName )

    if( ~PathFileName )
        return;
    end

    vHoles = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );
    
    fid_BrSightSeq = fopen( PathFileName, 'w' );
    
    fprintf( fid_BrSightSeq, '3\t2\n' );
    fprintf( fid_BrSightSeq, '1\n' );
    fprintf( fid_BrSightSeq, 'NZ\t1\n' );
    fprintf( fid_BrSightSeq, 'LPA\t2\n' );
    fprintf( fid_BrSightSeq, 'RPA\t4\n' );
    
    pItemsRef = [ find( [vHoles.Type] == get_PtTypeNo( oHelmet, 'HelmetRef' ) ) ];
    pItemsMtg = sort([ sMtg.v_pDet, sMtg.v_pSrc, sMtg.v_pEle ]);
    
    pItems = [ pItemsRef, pItemsMtg ];

    for( iItem=1:numel(pItems) )
        
        p = pItems(iItem);
        
        IsLast = ( iItem == numel(pItems) );
        
        if( vHoles(p).Type == get_PtTypeNo( oHelmet, 'HelmetRef' ) )
            fprintf( fid_BrSightSeq, '%s\t%d\n', vHoles(p).Label, vHoles(p).Type );
            fprintf( fid_BrSightSeq, '%s\t%d\n', vHoles(p).Label, vHoles(p).Type );
            if( ~IsLast ) 
                fprintf( fid_BrSightSeq, '%s\t%d\n', vHoles(p).Label, vHoles(p).Type );
            else
                fprintf( fid_BrSightSeq, '%s\t%d', vHoles(p).Label, vHoles(p).Type );
            end
        else
            if( ~IsLast ) 
                fprintf( fid_BrSightSeq, '%s\t%d\n', vHoles(p).Label, vHoles(p).Type );
            else
                fprintf( fid_BrSightSeq, '%s\t%d', vHoles(p).Label, vHoles(p).Type );
            end
        end
    end

    fclose(fid_BrSightSeq);