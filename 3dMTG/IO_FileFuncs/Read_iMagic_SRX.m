function [vertex_matrix, faces_matrix, NV, NT] = Read_iMagic_SRX( filename )

    fid_Srx=fopen(filename);
    [Packet, count] =  fread(fid_Srx, 2, 'uint32');
    NV = Packet(1);
    NT = Packet(2);
    [Coordinates, count] =  fread(fid_Srx, NV*3, 'single');
    [IndexBuffer, count] =  fread(fid_Srx, NT*3, 'uint32');
    XCoords = Coordinates(1:3:NV*3);
    YCoords = Coordinates(2:3:NV*3);
    ZCoords = Coordinates(3:3:NV*3);
    vertex_matrix = [[XCoords],[YCoords],[ZCoords]];
    faces_matrix = reshape( IndexBuffer+1, 3, NT )';
    %patch('Vertices',vertex_matrix,'Faces',faces_matrix);