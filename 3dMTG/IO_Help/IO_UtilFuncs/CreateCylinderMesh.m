function [VertexBuffer, IndexBuffer] = CreateCylinderMesh( Diameter, Depth, nPts )
    %Creates a cylinder mesh
    
    t = (0:(2*pi/nPts):2*pi)';
    
    Face1VertexBuffer = [ (Diameter/2)*cos(t), (Diameter/2)*sin(t), zeros(numel(t),1), ones(numel(t),1) ];
    Face1IndexCenter = size(Face1VertexBuffer,1)+1;
    Face1VertexBuffer( Face1IndexCenter, : ) = [ 0, 0, 0, 1 ];
    
    VertexBuffer = [ Face1VertexBuffer(:,:); ...
                     Face1VertexBuffer(:,1), Face1VertexBuffer(:,2), Face1VertexBuffer(:,3)-Depth, Face1VertexBuffer(:,4) ];

	Face1IndexBuffer = zeros( nPts, 5 );
    SidePanelsIndexBuffer = zeros( nPts, 5 );
    
    for( iPoly=1:nPts )
        Face1IndexBuffer(iPoly,:) = [ Face1IndexCenter, iPoly, iPoly+1, Face1IndexCenter, Face1IndexCenter ];
        SidePanelsIndexBuffer(iPoly,:) = [ iPoly, iPoly+1, iPoly+1+Face1IndexCenter, iPoly+Face1IndexCenter, iPoly ];
    end
    
    Face2IndexBuffer = Face1IndexBuffer + size( Face1VertexBuffer,1 );
    
    IndexBuffer = [ Face1IndexBuffer;...
                    Face2IndexBuffer; ...
                    SidePanelsIndexBuffer];
               
    
    
    
    