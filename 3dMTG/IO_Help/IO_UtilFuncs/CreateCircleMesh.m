function [VertexBuffer, IndexBuffer] = CreateCircleMesh( Diameter, nPts )
    %Creates a vector of points disposed in a circle on the XY plane
    %Each line of the matrix is a 4 element vector
    
    t = (0:(2*pi/nPts):2*pi)';
    
    VertexBuffer = [ (Diameter/2)*cos(t), (Diameter/2)*sin(t), zeros(numel(t),1), ones(numel(t),1) ];

    IndexBuffer = find( ones( size(VertexBuffer,1), 1 ) )';
    