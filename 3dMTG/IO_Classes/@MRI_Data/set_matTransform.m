%                 
% matTransform: 4x4 transformation matrix with translation components
%               at the bottom (vDisplay=v*matTransform);
%
% 
%
function oMRI = set_matTransform( oMRI, matTransform )
    oMRI.matCoRegistration = matTransform;
 