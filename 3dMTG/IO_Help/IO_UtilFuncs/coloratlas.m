function oMRI = coloratlas(oMRI, imgatlasname)

[VertexBuffer, IndexBuffer] = get_CortexMeshHiRes( oMRI );
vColor = zeros(size(VertexBuffer,1),1);% get_CortexHiResVcolor(oMRI);
VertexBuffer(:,4) = 1;
Vertex_tmp = VertexBuffer;
v_center = (Vertex_tmp(IndexBuffer(:,1),1:4) + Vertex_tmp(IndexBuffer(:,2),1:4)+ Vertex_tmp(IndexBuffer(:,3),1:4))/3;
try
    V = spm_vol(imgatlasname);
catch
    return
end
if 1 %function prenant les voxels voisins pour définir l'atlas (si segmentation est la même que le casque ça fonction très bien
    clear VertexBuffer
    clear IndexBuffer
    Vertex_r_x = round(Vertex_tmp(:,1));
    Vertex_r_y = round(Vertex_tmp(:,3));
    Vertex_r_z = round(Vertex_tmp(:,2));
    z_ind = [];
    Vertexlist =[];
    %les slices utilisés et leurs vertexs respectifs
    for z = 1:V.dim(3)
        if find(Vertex_r_z==z)
            Vertexlist = [Vertexlist,{find(Vertex_r_z==z)}];
            z_ind = [z_ind,z];
        end
    end
    for z = 1:size(z_ind,2) %V.dim(3)
        z_line = z_ind(z);
        %     z_line = V.dim(3)-z_line
        atlasz1 = spm_slice_vol(V,spm_matrix([0 0 (z_line-1)]),V.dim(1:2),0);
        atlasz2 = spm_slice_vol(V,spm_matrix([0 0 z_line]),V.dim(1:2),0);
        atlasz3 = spm_slice_vol(V,spm_matrix([0 0 (z_line+1)]),V.dim(1:2),0);
        listxy = Vertexlist{z};
        for list_ind = 1:size(listxy,1)
            vertex_ind = listxy(list_ind);
            x = Vertex_r_x(vertex_ind);
            y = Vertex_r_y(vertex_ind);
            y = V.dim(2)-y;
            x = V.dim(1)-x;
            c1 = atlasz1(x-1:x+1,y-1:y+1);
            c2 = atlasz2(x-1:x+1,y-1:y+1);
            c3 = atlasz3(x-1:x+1,y-1:y+1);
            atlas = [ c1,c2,c3];
            id = find(atlas);
            atlas = sort(atlas(id));
            if ~isempty(id)
                clear val
                for i=1:numel(id)
                    val(i)=sum(atlas==atlas(i));
                end
                [x,j]=max(val);
                c = atlas(j);
            else
                c=0;
            end
            vColor(vertex_ind,1) = c;
        end
        clear atlasz
    end

else %fonction projetent la normal de la surface sur l'atlas (idéal pour la peau ! )
    clear VertexBuffer
    clear IndexBuffer
    Vertex_r_x = V.dim(1)-round(Vertex_tmp(:,1)); % pour etre dans le même référentiel que la résonance
    Vertex_r_y = V.dim(2)-round(Vertex_tmp(:,3));
    Vertex_r_z = round(Vertex_tmp(:,2));

    % Normal
    norm = sqrt((Vertex_r_x-V.dim(1)/2).^2+(Vertex_r_y-V.dim(2)/2).^2+(Vertex_r_z-V.dim(3)/2).^2); 
    
    Normal_x = (Vertex_r_x-V.dim(1)/2)./norm;
    Normal_y = (Vertex_r_y-V.dim(2)/2)./norm;
    Normal_z = (Vertex_r_z-V.dim(3)/2)./norm;   
    
    MaxDepth_vox = 20; % pour eviter le cas ou la normal ne touche jamais le cortex
    DepthSteps = 1:1:MaxDepth_vox;
    tic

    for pElem=1:numel(Normal_x);
        pElem
        HCoordsStepsX =Vertex_r_x(pElem)-Normal_x(pElem)*DepthSteps;
        HCoordsStepsY = Vertex_r_y(pElem)-Normal_y(pElem)*DepthSteps;
        HCoordsStepsZ = Vertex_r_z(pElem)-Normal_z(pElem)*DepthSteps;
          
        for( iDepth=1:numel(DepthSteps) )
            z_line = HCoordsStepsZ(iDepth);
            atlasz1 = spm_slice_vol(V,spm_matrix([0 0 (z_line)]),V.dim(1:2),0);
            c = 0;
            if atlasz1(round(HCoordsStepsX(iDepth)),round(HCoordsStepsY(iDepth)))> 0
                c = atlasz1(round(HCoordsStepsX(iDepth)),round(HCoordsStepsY(iDepth)));
                break
            end
        end
         vColor(pElem,1) = c;
    end
    toc
end
%figure;plot(vColor)
    oMRI = set_CortexHiResVcolor(oMRI,vColor);