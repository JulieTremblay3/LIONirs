%                 
% matTransformManual: 4x4 transformation matrix with translation components
%               at the bottom, only part of the manual displacement is save here
%               and will be use to displace mri in the skin projection calcul;
%               unit in voxel
% 
%
function oMRI = set_matTransformManual( oMRI, matTransformManual )
    oMRI.matTransformManual = matTransformManual;
 