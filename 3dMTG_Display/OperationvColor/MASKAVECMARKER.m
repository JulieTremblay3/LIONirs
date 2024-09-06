%Creation D'un masque de la région de Broca
filevertex ='F:\Fluence\MTG_MRI\NORS_NIRS\skinlow.SRX'%ouvrir la topo pour avoir les vertex
%utiliser un ma
filevmrk = 'F:\Fluence\MTG_MRI\NORS_NIRS\L_Wernicke.MRK'; %suprimer les premières lignes

[VertexBuffer faces_matrix, NV, NT] = Read_iMagic_SRX( filevertex )



%Open VMRK
fid=fopen(filevmrk)
fgetl(fid)
fgetl(fid)
fgetl(fid)
i = 1;
while 1
    line = fgetl(fid)
    if line ==-1
       break
    end
    A(i,:)=str2num(line);
    i = i+1;
end
fclose(fid)

ivcolor = [];
dx = 3; 
dy = 3;
dz = 3;

%Trouver les vertexs proches des markeurs et mettre la le vcolor à 1 pour
%faire la régions d'intérêts sur la Topographie. 
for i=1:size(A,1)
    id = find((A(i,1)-dx)< VertexBuffer(:,1) &...
        (A(i,1)+dx)> VertexBuffer(:,1)&...
        (A(i,2)-dy)< VertexBuffer(:,2) &...
        (A(i,2)+dy)> VertexBuffer(:,2)&...
        (A(i,3)-dz)< VertexBuffer(:,3) &...
        (A(i,3)+dz)> VertexBuffer(:,3))       
        ivcolor = [ivcolor;id];
end
vColor = zeros(size( VertexBuffer,1),1);
vColor(ivcolor)=1;
save([filevmrk(1:end-3),'vcolor'],'vColor','-mat')


%Trouver les vertexs proches des markeurs et mettre la le vcolor à 1 pour
%faire la régions d'intérêts sur la Topographie symétrique 
xsym = 181;
ivcolor  = [];

for i=1:size(A,1)
    id = find((181-A(i,1)-dx)< VertexBuffer(:,1) &...
        (181-A(i,1)+dx+xsym)> VertexBuffer(:,1)&...
        (A(i,2)-dy)< VertexBuffer(:,2) &...
        (A(i,2)+dy)> VertexBuffer(:,2)&...
        (A(i,3)-dz)< VertexBuffer(:,3) &...
        (A(i,3)+dz)> VertexBuffer(:,3))       
        ivcolor = [ivcolor;id];
end
vColor = zeros(size( VertexBuffer,1),1);
vColor(ivcolor)=1;
save([filevmrk(1:end-4),'symetric','.vcolor'],'vColor','-mat')
