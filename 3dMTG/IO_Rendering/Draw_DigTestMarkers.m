function Draw_DigTestMarkers( oDispOpt, oInterDlgComm, oHelmet )

    persistent v_hPrevDispItems;
    
    %Effacer les anciens elements graphiques avant d'afficher les nouveaux
    if( ~isempty( v_hPrevDispItems ) )
        delete( v_hPrevDispItems(find( ishandle(v_hPrevDispItems) & (v_hPrevDispItems ~= 0) )));
        v_hPrevDispItems = [];
    end
    
    if( get_DispOptChecked( oDispOpt, 'Dig_SubjectRegistrationTestMarkers' ) )
  
        vColor = get_ItemColor( oDispOpt, 'Dig_SubjectRegistrationTestMarkers' );
    
        Diameter = 0.005; %in meters
        nParallels_Meridians = 5; %5Parallels and 5 meridians (25 polygons)
        [VertexBuffer, IndexBuffer] = CreateSphereMesh( Diameter, nParallels_Meridians );
        
        %Matrice de transformation du MRI (transformation globale)
        vHoles = get_vHoles( oHelmet );
        vMarkers = vHoles( [vHoles.Type] == get_PtTypeNo(oHelmet, 'SkinFitTest' ) );
        
        for( iMrk=1:numel(vMarkers) )
            
            %Translation de la sphere a la position du marker
            tmpVertexBuffer = [ VertexBuffer(:,1)+vMarkers(iMrk).Coord.x, VertexBuffer(:,2)+vMarkers(iMrk).Coord.y, VertexBuffer(:,3)+vMarkers(iMrk).Coord.z ];     
        
            v_hPrevDispItems(iMrk) = patch( 'Vertices',  tmpVertexBuffer, ...
                                            'Faces', IndexBuffer, ...
                                            'FaceColor', vColor, ...
                                            'EdgeColor', 'none', ...
                                            'FaceLighting', 'gouraud' );
        end
    end
    