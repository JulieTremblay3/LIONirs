function v_Holes = OrientationVHoles(v_Holes,type)

matvHoles = ones(numel(v_Holes),4);
matvHolesNorm =  ones(numel(v_Holes),4);
for ind = 1:numel(v_Holes)
    matvHoles(ind,1)=v_Holes(ind).Coord.x;
    matvHoles(ind,2)=v_Holes(ind).Coord.y;
    matvHoles(ind,3)=v_Holes(ind).Coord.z;
    matvTransformation = v_Holes(ind).Transformation;
    matvHolesNorm(ind,1)=v_Holes(ind).Normal.x;
    matvHolesNorm(ind,2)=v_Holes(ind).Normal.y;
    matvHolesNorm(ind,3)=v_Holes(ind).Normal.z;
end
alpha = pi/6;
mRpa = makehgtform('xrotate', alpha)
matvHoles2 = matvHoles * mRpa';  
matvHolesNorm2 =  matvHolesNorm * mRpa'; 

for ind = 1:numel(v_Holes)
%    v_Holes(ind).Transformation = matvTransformation*mRpa';
    v_Holes(ind).Coord.x = matvHoles2(ind,1);
    v_Holes(ind).Coord.y = matvHoles2(ind,2);
    v_Holes(ind).Coord.z = matvHoles2(ind,3);
%     v_Holes(ind).Normal.x = matvHolesNorm2(ind,1);
%     v_Holes(ind).Normal.y = matvHolesNorm2(ind,2);
%     v_Holes(ind).Normal.z = matvHolesNorm2(ind,3);
end