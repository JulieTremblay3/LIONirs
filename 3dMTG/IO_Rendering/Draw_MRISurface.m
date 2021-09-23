function Draw_MRISurface( oDispOpt, oInterDlgComm, oMRI )

    persistent v_hPrevDispItems;
     
    %Effacer les anciens elements graphiques avant d'afficher les nouveaux
    if( ~isempty( v_hPrevDispItems ) )
        delete( v_hPrevDispItems(find( ishandle(v_hPrevDispItems) & (v_hPrevDispItems ~= 0) )));
        v_hPrevDispItems = [];
    end
 
    if( get_DispOptChecked( oDispOpt, 'MRI_SurfaceCortexHiRes' ) )
        [VertexBuffer, IndexBuffer] = get_CortexMeshHiRes( oMRI );
        %vColor = get_ItemColor( oDispOpt, 'MRI_SurfaceCortex' );
        vColor = get_CortexHiResVcolor(oMRI);
        atlas_zone = str2num(get_DispAtlasZone(oDispOpt));
        nb_color = max(vColor);
        if isempty(vColor)
            return
        end
        if nb_color == 0
            colormap([88/256 88/256 88/256]);
        else
            mapcolorcube = colorcube(nb_color*2);
            map = mapcolorcube(1:nb_color,:); 
            map(18,:)=[0,64/256,0];
            map(19,:)=[0,0,1];
            map(20,:)=[1,0,0];
            caxis([0,nb_color]);
        end
      if ~isempty(atlas_zone)
          map = lines(100) %ones(nb_color,3)*(159/256); 
          nb_enumatlas = numel(atlas_zone);
          mixcolor = lines(nb_enumatlas); 
          for ind_map = 1:numel(atlas_zone)
              map(atlas_zone(ind_map)+1,:) = mixcolor(ind_map,:);
%               if ind_map == 17
%                    map(atlas_zone(ind_map)+1,:)= 
          end
          colormap(map);      
      else
          p=mfilename('fullpath');
          [pathstr, name, ext]=fileparts(p);
          [num, txt, raw] = xlsread(fullfile(pathstr, 'Brodmann48ColorDisplay.xls'));
          %map = [0,0,0]
         map = num(:,3:5);
          map = lines(100) %ones(nb_color,3)*(159/256); 
          colormap(map); 
      end
                        
    elseif( get_DispOptChecked( oDispOpt, 'MRI_SurfaceCortexLowRes' ) )
        [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI );
        vColor = get_ItemColor( oDispOpt, 'MRI_SurfaceCortex' );   
        colormap([88/256 88/256 88/256]);
    elseif( get_DispOptChecked( oDispOpt, 'MRI_SurfaceSkin' ) )
        [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
        vColor = get_ItemColor( oDispOpt, 'MRI_SurfaceSkin' );
        colormap([255/256 255/256 179/256]);
        vColor = [255/256 255/256 179/256];

    else
        return;
    end
 
    
 if( ~isempty( VertexBuffer) && ~isempty( IndexBuffer) )
    VertexBuffer(:,4) = 1;
    T = get_matTransform(oMRI);
    Vertex_tmp = VertexBuffer*T;   
    [m n] = size(IndexBuffer);
%     Vertex_tmp2 =  Vertex_tmp*T^-1;   
%     vColor = zeros(m,1); 


%     if get_DispOptChecked( oDispOpt, 'MRI_Image' )
%          name = get_DispOptMRIImage( oDispOpt )
%          [n,m] = size(name);         
%          for nb_load = 1:n;
%                 load(name(nb_load,:),'-mat');
%                 matTransform = DOT{1,1}.matTransform;
%                 Image = DOT{1,1}.IMG.HbO;
%                 v_center = (Vertex_tmp(IndexBuffer(:,1),1:4) + Vertex_tmp(IndexBuffer(:,2),1:4)+ Vertex_tmp(IndexBuffer(:,3),1:4))/3;
%                 v_center_trp = v_center * matTransform;
%                 x_min =  DOT{1}.IMG.Medium.CompVol.X(1)/100;
%                 x_max =  DOT{1}.IMG.Medium.CompVol.X(end)/100;
%                 y_min =  DOT{1}.IMG.Medium.CompVol.Y(1)/100;
%                 y_max =  DOT{1}.IMG.Medium.CompVol.Y(end)/100;
%                 mat_SrcDet = [DOT{1}.SD.SrcPos(:,3);DOT{1}.SD.DetPos(:,3)];
%                 ind = find(mat_SrcDet(:));
%                 z_min = min(mat_SrcDet(ind))/100 - 0.01; 
%                 z_max = max(mat_SrcDet(ind))/100 + 0.01;    
%                  %   par pixel d'affichage
%                 [num_y num_x]= size(Image);
%                 for x = 1:(num_x)
%                     x_min1 = x_min + (x-1)*(x_max-x_min)/num_x;
%                     x_max1 = x_min + (x)*(x_max-x_min)/num_x;
%                     for y = 1:(num_y)
%                         y_min1 = y_min + (y-1)*(y_max-y_min)/num_y;
%                         y_max1 = y_min +(y)*(y_max-y_min)/num_y;
%                          p_in = find(v_center_trp(:,1) > x_min1 & v_center_trp(:,1) < x_max1 & v_center_trp(:,2) > y_min1 & v_center_trp(:,2) < y_max1 & v_center_trp(:,3) > z_min & v_center_trp(:,3) < z_max);
%                          for i = 1 : numel(p_in)
%                             vColor(p_in(i),1) = Image(y,x);
%                          end   
% 
%                     end
%                 end 
%          end
%     end   


        if( get_DispOptChecked( oDispOpt, 'App_Transparency' ) )
            v_hPrevDispItems = patch('Vertices',Vertex_tmp(:,1:3),...
                                     'Faces',IndexBuffer, ...
                                     'EdgeColor', 'none', ...
                                     'FaceVertexCData', vColor, ...
                                     'Facecolor','flat',...
                                     'FaceLighting', 'gouraud', ...
                                     'SpecularColorReflectance', 0.5, ...
                                     'DiffuseStrength', 0.5, ...
                                     'SpecularStrength', 0.5, ...
                                     'FaceAlpha', 0.5, ...
                                     'AmbientStrength', 0.5);
       
        else
            
            v_hPrevDispItems = patch('Vertices',Vertex_tmp(:,1:3),...
                                     'Faces',IndexBuffer, ...
                                     'EdgeColor', 'none', ...
                                     'FaceVertexCData', vColor, ...
                                     'Facecolor','flat',...
                                     'FaceLighting', 'gouraud', ...
                                     'SpecularColorReflectance', 0.5, ...
                                     'DiffuseStrength', 0.5, ...
                                     'SpecularStrength', 0.5, ...
                                     'AmbientStrength', 0.5);                                
        end
               colormap([255/256 255/256 179/256]);

end