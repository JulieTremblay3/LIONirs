function Draw_HoleSurfaceIntersection( oDispOpt, oInterDlgComm, oHelmet )

    persistent v_hPrevDispItems;
    
    %Effacer les anciens elements graphiques avant d'afficher les nouveaux
    if( ~isempty( v_hPrevDispItems ) )
        delete( v_hPrevDispItems(find( ishandle(v_hPrevDispItems) & (v_hPrevDispItems ~= 0) )));
        v_hPrevDispItems = [];
    end
    
    vColor = get_ItemColor( oDispOpt, 'MRI_SurfaceHolesIntersectionMarkers' );
    
    if( get_DispOptChecked( oDispOpt, 'MRI_SurfaceHolesIntersectionMarkers' ) )
        %[VertexBuffer, IndexBuffer] = CreateSphereMesh( 0.0054, 4 );
        [VertexBuffer, IndexBuffer] = CreateSphereMesh( 0.01, 16 );
        vHoles = get_vHoles( oHelmet );
        
        %Boucle d'affichage de point d'intersection entre le cortex et la normale
        if( get_DispOptChecked( oDispOpt, 'MRI_SurfaceCortexHiRes' ) || get_DispOptChecked( oDispOpt, 'MRI_SurfaceCortexLowRes' ) )
            for( p=1:numel(vHoles) )
                v_hPrevDispItems(p) = patch( 'Vertices',  [VertexBuffer(:,1)+(vHoles(p).Coord.x-vHoles(p).Normal.x*vHoles(p).CortexDepth),...
                                                           VertexBuffer(:,2)+(vHoles(p).Coord.y-vHoles(p).Normal.y*vHoles(p).CortexDepth), ...
                                                           VertexBuffer(:,3)+(vHoles(p).Coord.z-vHoles(p).Normal.z*vHoles(p).CortexDepth)], ...
                                             'Faces', IndexBuffer, ...
                                             'FaceColor', vColor, ...
                                             'EdgeColor', 'none', ...
                                             'FaceLighting', 'gouraud' );
            end
        %Boucle d'affichage de point d'intersection entre la peau et la normale
        elseif get_DispOptChecked( oDispOpt, 'MRI_SurfaceSkin' )
            for( p=1:numel(vHoles) )
                v_hPrevDispItems(p) = patch( 'Vertices', [VertexBuffer(:,1)+(vHoles(p).Coord.x-vHoles(p).Normal.x*vHoles(p).SkinDepth),...
                                                          VertexBuffer(:,2)+(vHoles(p).Coord.y-vHoles(p).Normal.y*vHoles(p).SkinDepth), ...
                                                          VertexBuffer(:,3)+(vHoles(p).Coord.z-vHoles(p).Normal.z*vHoles(p).SkinDepth)], ...
                                             'Faces', IndexBuffer, ...
                                             'FaceColor', vColor, ...
                                             'EdgeColor', 'none', ...
                                             'FaceLighting', 'gouraud' );
            end
        end
    end
    