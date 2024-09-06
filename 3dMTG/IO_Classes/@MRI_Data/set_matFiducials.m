%                 _                _
% matFiducials = |  NASx NASy NASz  |
%                |  LPAx LPAy LPAz  |
%                |_ RPAx RPAy RPAz _|
%
function oMRI = set_matFiducials( oMRI, matFiducials )
    oMRI.matFiducials = matFiducials;
 