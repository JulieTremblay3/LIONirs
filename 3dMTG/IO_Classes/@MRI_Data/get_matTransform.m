%                 
% matTransform: 4x4 transformation matrix with translation components
%               at the bottom (vDisplay=v*matTransform);
%
% 
%
function matTransform = get_matTransform( oMRI )
    matTransform = oMRI.matCoRegistration;
 