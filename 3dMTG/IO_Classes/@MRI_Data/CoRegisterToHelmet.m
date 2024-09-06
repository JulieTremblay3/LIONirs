% function CoRegisterToHelmet
%
% Permet de construire la matrice de coregistration interne servant
% a alligner le MRI ainsi que les surfaces et les segmentation sur l'espace
% des coordonnees du casque (grace aux fiducies: NAS, LPA, RPA). 
%
function oMRI = CoRegisterToHelmet( oMRI, oHelm )
    
    sMtg = get_Mtg(oHelm);
    
    %References communes: [NAS;LPA;RPA]
    matSourceRefs = oMRI.matFiducials;
    matDestinationRefs = sMtg.matFiducials;
    
    %Dimensions physiques
    vSourcePhysicalDim = get_PhysicalVoxDim( oMRI );
    vDestPhysicalDim = [1,1,1];
    
%     matSourceRefs      = [154,198,95;305,220,265;16,206,270];
%     matDestinationRefs = sMtg.matFiducials;
%     vSourcePhysicalDim = [.0005,.0005,.0005];
%     vDestPhysicalDim   = [1,1,1];
    %Coregistration 6DOF "Rigid Body"
    oMRI.matCoRegistration = Create_6DOF_CoRegistration_Matrix( matSourceRefs, matDestinationRefs, vSourcePhysicalDim, vDestPhysicalDim, false, true );
    
    
 