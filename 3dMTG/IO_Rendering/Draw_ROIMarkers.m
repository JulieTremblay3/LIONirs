function Draw_ROIMarkers( oDispOpt, oInterDlgComm, oMRI )

    persistent v_hPrevDispItems;
    
    %Effacer les anciens elements graphiques avant d'afficher les nouveaux
    if( ~isempty( v_hPrevDispItems ) )
        delete( v_hPrevDispItems(find( ishandle(v_hPrevDispItems) & (v_hPrevDispItems ~= 0) )));
        v_hPrevDispItems = [];
    end
    
    vMarkers = get_vMarkers( oMRI );
    if( get_DispOptChecked( oDispOpt, 'MRI_ROIMarkers' ) && ( (numel(vMarkers) > 1) || ~isempty( find( vMarkers(1).Coord ) ) ) )
  
        Diameter = 0.004; %in meters
        nParallels_Meridians = 6; %6Parallels and 6 meridians (36 polygons)
        [VertexBuffer, IndexBuffer] = CreateSphereMesh( Diameter, nParallels_Meridians );
        
        %Matrice de transformation du MRI (transformation globale)
        T = get_matTransform( oMRI );
        
        FontWeight = 'normal';
        if( get_DispOptChecked( oDispOpt, 'App_BoldLabels' ) )
            FontWeight = 'bold';
        end
        vColor = get_ItemColor( oDispOpt, 'MRI_ROIMarkers' );
    
        for( iMrk=1:numel(vMarkers) )
            
            %Matrice de transformation du MRI (transformation globale)
            Position = vMarkers(iMrk).Coord*T;
            
            %Translation de la sphere a la position du marker
            tmpVertexBuffer = [ VertexBuffer(:,1)+Position(1), VertexBuffer(:,2)+Position(2), VertexBuffer(:,3)+Position(3) ];     
        
            v_hPrevDispItems((iMrk*2)-1) = patch( 'Vertices',  tmpVertexBuffer, ...
                                                  'Faces', IndexBuffer, ...
                                                  'FaceColor', vColor, ...
                                                  'EdgeColor', 'none', ...
                                                  'FaceLighting', 'gouraud' );
            v_hPrevDispItems((iMrk*2)) = text(  Position(1), Position(2),  Position(3), ...
                                                vMarkers(iMrk).Label, 'FontSize',10, 'Color', vColor, ...
                                                'FontWeight', FontWeight, ...
                                                'VerticalAlignment', 'top');
        end
    end
    