function [VertexBuffer, IndexBuffer] = CreateSphereMesh( Diameter, nFaces )
    %Creates a vector of points disposed in a circle on the XY plane
    %Each line of the matrix is a 4 element vector
    
    
    [x,y,z] = sphere(nFaces);
    Vx = x';
    Vy = y';
    Vz = z';
    VertexBuffer = [Vx(:).*Diameter/2, Vy(:).*Diameter/2, Vz(:).*Diameter/2, ones(numel(Vz(:)),1)];

    IndexBuffer = [];
    nt = 0;
    for(iH=1:nFaces+1)
        for(iV=1:nFaces)
            iH1 = iH;
            iH2 = iH+1;
            if( iH2 > nFaces+1 )
                iH2 = 1;
            end
            IndexBuffer(nt+1,:) = [(iV-1)*(nFaces+1)+iH1, (iV-1)*(nFaces+1)+iH2, (iV)*(nFaces+1)+iH2, (iV)*(nFaces+1)+iH1, (iV-1)*(nFaces+1)+iH1];
            nt = nt+1;
        end
    end