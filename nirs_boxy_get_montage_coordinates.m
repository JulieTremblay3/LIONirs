function [mat_Mtg mat_Ele] = nirs_boxy_get_montage_coordinates(pSrc,...
    pDet,pEle,holes,HolesMtg)
%get real world x, y, z coordinates of sources, detectors & electrodes
p_montage = [pSrc,pDet,pEle];
max_Col= size(p_montage,2);
mat_Mtg = zeros(max_Col,1,4); %pre-allocate

for i = 1:max_Col
    %p: index that identifies the source, detector or electrode on helmet
    p = p_montage(i);
    mat_Mtg(i,1,1) = holes(p).Coord.x;
    mat_Mtg(i,1,2) = holes(p).Coord.y;
    mat_Mtg(i,1,3) = holes(p).Coord.z;
    mat_Mtg(i,1,4) = HolesMtg(p);   
end %end for i = 1:max_Col

%Separate the electrodes from the sources and detectors
mat_Ele =mat_Mtg((size([pSrc,pDet],2)+1):max_Col,1,:);
mat_Mtg =mat_Mtg(1:size([pSrc,pDet],2),1,:);
end
