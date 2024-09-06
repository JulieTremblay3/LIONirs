%For every Helmet Hole, the position and the normal is 
%calculated in MRI Voxel-space and returned as two vector outputs. The position of
%elements in Helm_MRI_Norms and Helm_MRI_Coords is in the same order as in
%Helmet::vHoles.
%Projection : 1 = to helmet, 2= to skin, 3 = to cortex
function [Helm_MRI_Coords, Helm_MRI_Norms ] = CoRegisterHelmetToMRI( oHelmet, oMRI_Data,projection )

    vHolesH = get_vHoles( oHelmet );
    sMtg = get_Mtg(oHelmet);

    %Transfert des coordonnees sous forme de matrices
    CoordsH = [vHolesH.Coord];
    NormalsH = [vHolesH.Normal];
    matNormalsH = [NormalsH.x;NormalsH.y;NormalsH.z]';
    if projection == 1   
        matCoordsH = [CoordsH.x;CoordsH.y;CoordsH.z]';
    elseif projection == 2       
        for p = 1:numel(vHolesH(:))
                   CoordsH(p).x = vHolesH(p).Coord.x-vHolesH(p).Normal.x*vHolesH(p).SkinDepth;
                   CoordsH(p).y = vHolesH(p).Coord.y-vHolesH(p).Normal.y*vHolesH(p).SkinDepth;
                   CoordsH(p).z = vHolesH(p).Coord.z-vHolesH(p).Normal.z*vHolesH(p).SkinDepth;
        end
        matCoordsH = [CoordsH.x;CoordsH.y;CoordsH.z]';
    elseif projection == 3   
        for p = 1:numel(vHolesH(:))
                   CoordsH(p).x = vHolesH(p).Coord.x-vHolesH(p).Normal.x*vHolesH(p).CortexDepth;
                   CoordsH(p).y = vHolesH(p).Coord.y-vHolesH(p).Normal.y*vHolesH(p).CortexDepth;
                   CoordsH(p).z = vHolesH(p).Coord.z-vHolesH(p).Normal.z*vHolesH(p).CortexDepth;
        end
        matCoordsH = [CoordsH.x;CoordsH.y;CoordsH.z]';  
    end 
    
    matHelmetRefs = sMtg.matFiducials;
    matMRIRefs = get_matFiducials( oMRI_Data );
    if matMRIRefs == 0
       msgbox('Your mri .vox don''t have fiducial')
       return 
    end
	%Construction de la matrice de coregistration
    IsHelmetRH = true;
    IsMRIRH = false;
    Physical_HelmCoordSystemSize = [1,1,1]; %Le casque provient d'un systeme orthonorme et metrique
    matCoordsH(:,4) = 1;
    matNormalsH(:,4) = 1;

    %Transformation a 6 Degres de liberte (Rigid body transform)
	%(6DOF = rotation xyz et translation xyz)
    %Application de la matrice de coregistration pour les positions
% 	matMRIRefs(:,4) = 1
%     matMRIRefs * (matCo-1)

    Helm_MRI_Coords = matCoordsH*Create_6DOF_CoRegistration_Matrix( matHelmetRefs,                matMRIRefs, ...
                                                                Physical_HelmCoordSystemSize, get_PhysicalVoxDim( oMRI_Data ), ...
                                                                 IsHelmetRH,                   IsMRIRH );
%     matMRIRefs = matMRIRefs.*0.001;
%     Helm_MRI_Coords = matCoordsH*(Create_6DOF_CoRegistration_Matrix( matHelmetRefs,                matMRIRefs, ...
%                                                                 Physical_HelmCoordSystemSize, Physical_HelmCoordSystemSize, ...
%                                                        IsHelmetRH,IsMRIRH ));
% %   matCo = get_matTransform(oMRI_Data);    
%  defaultMRITF  = Create_6DOF_CoRegistration_Matrix( matHelmetRefs,...
%          matMRIRefs, Physical_HelmCoordSystemSize, get_PhysicalVoxDim( oMRI_Data ), ...
%          IsHelmetRH, IsMRIRH )
% %  Helm_MRI_Coords=  matCoordsH*matCo^-1;
% %  defaultMRITF - matCo^-1
%       Helm_MRI_Coords = matCoordsH * defaultMRITF^-1;    
      
      
    %Transformation a 3 Degres de liberte (Affine transform)
	%(3DOF = rotation xyz) <-Bref, pas de translation pour les normales
    %Application de la matrice de coregistration pour les normales.
    %Pas de remise a l'echelle car les 2 systemes sont cubiques
    %(anisotropiques??).
    Helm_MRI_Norms = matNormalsH*Create_3DOF_CoRegistration_Matrix( matHelmetRefs,                matMRIRefs, ...
        IsHelmetRH,                   IsMRIRH );
    

    
    %Normales unitaires (en voxel-space)
    for( iNorm=1:size(Helm_MRI_Norms,1) )
        if( norm(Helm_MRI_Norms(iNorm,1:3)) )
            Helm_MRI_Norms(iNorm,1:3) = Helm_MRI_Norms(iNorm,1:3) / norm(Helm_MRI_Norms(iNorm,1:3));
        end
    end