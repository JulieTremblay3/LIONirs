% Classe de gestion des options d'affichage
function oDispOpt = IO_DisplayOptions(args)

    %Constructeur par defaut
    if( nargin == 0 )
        
        %MATLAB n'offre pas la possibilite de definir des constantes.
        %Ainsi, les champs ci-bas (indices) servent de constantes
        %afin de definir chaque element 
        oDispOpt.Subscripts.Dig_HelmetHoles                     = 01;
        oDispOpt.Subscripts.Dig_HelmetHolesNormals              = 02;
        oDispOpt.Subscripts.Dig_HelmetRegistrationPoints        = 03;
        oDispOpt.Subscripts.Dig_SubjectFiducials                = 04;
        oDispOpt.Subscripts.Dig_SubjectRegistrationTestMarkers  = 05;
        oDispOpt.Subscripts.Dig_CoordinateSystemAxes            = 06;
        oDispOpt.Subscripts.Lbl_HelmetHoleId                    = 07;
        oDispOpt.Subscripts.Lbl_OptodesHoleId                   = 08;
        oDispOpt.Subscripts.Lbl_OptodesFibersId                 = 09;
        oDispOpt.Subscripts.App_BoldLabels                      = 10;
        oDispOpt.Subscripts.App_BoldLines                       = 11;
        oDispOpt.Subscripts.App_InvertedColorScheme             = 12;
        oDispOpt.Subscripts.App_ItemsAlwaysVisible              = 13;
        oDispOpt.Subscripts.App_Reflections                     = 14;
        oDispOpt.Subscripts.App_Transparency                    = 15;
        oDispOpt.Subscripts.Opt_FittedOnSkin                    = 16;
        oDispOpt.Subscripts.Opt_CortexIntersectionMarker        = 17;
        oDispOpt.Subscripts.Opt_MuxChannel                      = 18;
        oDispOpt.Subscripts.Ele_PositionOnHelmet                = 19;
        oDispOpt.Subscripts.MRI_SurfaceCortexHiRes              = 20;
        oDispOpt.Subscripts.MRI_SurfaceCortexLowRes             = 21;
        oDispOpt.Subscripts.MRI_SurfaceSkin                     = 22;
        oDispOpt.Subscripts.MRI_SurfaceHolesIntersectionMarkers = 23;
        oDispOpt.Subscripts.MRI_ROIMarkers                      = 24;
        oDispOpt.Subscripts.MRI_SubjectFiducials                = 25;
        oDispOpt.Subscripts.MRI_Image                           = 26;
        
        oDispOpt.Subscripts.SphereHoles                         = 27;
        oDispOpt.Subscripts.line_holes                          = 28;
        oDispOpt.Subscripts.App_BackgroundColor                 = 29;

        %Nombre d'elements dans la strcuture d'indices
        oDispOpt.NbItems = length( fields(oDispOpt.Subscripts) );
        
        %Vecteur de booleens: 
        % o=item inactif (unchecked)
        % 1=item actif (checked)
        oDispOpt.v_bItemChk = zeros( 1, oDispOpt.NbItems );
        
        %Initialisation: Ces elements sont actives par defaut
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Dig_HelmetHoles                    ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Dig_SubjectFiducials               ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Dig_SubjectRegistrationTestMarkers ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.App_BoldLabels                     ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Opt_FittedOnSkin                   ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Opt_CortexIntersectionMarker       ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Ele_PositionOnHelmet               ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.MRI_SurfaceSkin                    ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.MRI_ROIMarkers                     ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Ele_PositionOnHelmet               ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Lbl_HelmetHoleId                   ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Lbl_OptodesHoleId                  ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.Lbl_OptodesFibersId                ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.App_Reflections                    ) = true; % Items 
        oDispOpt.v_bItemChk(oDispOpt.Subscripts.App_BackgroundColor                ) = false; % Items 

        %Intervalle d'affichage du string d'identification de trou
        oDispOpt.Lbl_HelmetHoleId_Interval = 5;
        
        oDispOpt.AtlasZone = '';
        oDispOpt.SphereDiameter = '0.5';
        % Modèles de couleurs                                         =    [  R ,  G ,  B  ] 
        %-------------------------------------------------------
        % COLOR SCHEMES 1 and 2
        %-------------------------------------------------------
        oDispOpt.ColorSchemes(:,:).Background =                          [ [ 1.0, 1.0, 1.0 ];... %Color Scheme 1
                                                                           [ 0.0, 0.0, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Opt_Detector =                        [ [ 0.0, 0.7, 0.0 ];... %Color Scheme 1
                                                                           [ 0.0, 0.7, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Opt_Source =                          [ [ 0.9, 0.0, 0.0 ];... %Color Scheme 1
                                                                           [ 0.9, 0.0, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Dig_Holes_Complete =                  [ [ 0.0, 0.0, 0.0 ];... %Color Scheme 1
                                                                           [ 1.0, 1.0, 1.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Dig_Holes_Partial =                   [ [ 0.8, 0.0, 0.0 ];... %Color Scheme 1
                                                                           [ 0.8, 0.0, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Opt_LabelSrcHole =                    [ [ 0.1, 0.8 ,0.7 ];... %Color Scheme 1
                                                                           [ 0.1, 0.8 ,0.7 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Opt_LabelSrcFiber =                   [ [ 0.2, 0.2, 1.0 ];... %Color Scheme 1
                                                                           [ 0.2, 0.2, 1.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Opt_LabelDetHole =                    [ [ 1.0, 0.4, 0.4 ];... %Color Scheme 1
                                                                           [ 1.0, 0.4, 0.4 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Opt_LabelDetFiber =                   [ [ 1.0, 0.3, 0.1 ];... %Color Scheme 1
                                                                           [ 1.0, 0.3, 0.1 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Selection =                           [ [ 0.0, 0.0, 1.0 ];... %Color Scheme 1
                                                                           [ 0.0, 0.0, 1.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Identification =                      [ [ 0.0, 0.0, 1.0 ];... %Color Scheme 1
                                                                           [ 0.0, 0.0, 1.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Ruler =                               [ [ 0.8, 0.8, 0.0 ];... %Color Scheme 1
                                                                           [ 0.8, 0.8, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Ele_PositionOnHelmet =                [ [ 0.8, 0.8, 0.0 ];... %Color Scheme 1
                                                                           [ 0.8, 0.8, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Lbl_HelmetHoleId =                    [ [ 0.1, 0.1, 1.0 ];... %Color Scheme 1
                                                                           [ 0.1, 0.1, 1.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Dig_SubjectRegistrationTestMarkers =  [ [ 0.0, 0.8, 0.0 ];... %Color Scheme 1
                                                                           [ 0.0, 0.8, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Dig_CoordinateSystemAxes =            [ [ 0.0, 0.0, 1.0 ];... %Color Scheme 1
                                                                           [ 0.0, 0.0, 1.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Dig_HelmetRef =                       [ [ 0.0, 0.8, 0.0 ];... %Color Scheme 1
                                                                           [ 0.0, 0.8, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).Dig_Fiducials =                       [ [ 0.0, 0.8, 0.0 ];... %Color Scheme 1
                                                                           [ 0.0, 0.8, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).MRI_Fiducials =                       [ [ 0.9, 0.0, 0.0 ];... %Color Scheme 1
                                                                           [ 0.9, 0.0, 0.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).MRI_SurfaceCortex =                   [ [ 0.5, 0.5, 0.5 ];... %Color Scheme 1
                                                                           [ 0.5, 0.5, 0.5 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).MRI_SurfaceSkin =                     [ [ 0.5, 0.5, 0.5 ];... %Color Scheme 1
                                                                           [ 0.5, 0.5, 0.5 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).MRI_SurfaceHolesIntersectionMarkers = [ [ 0.1, 0.1, 1.0 ];... %Color Scheme 1
                                                                           [ 0.1, 0.1, 1.0 ] ];  %Color Scheme 2
        oDispOpt.ColorSchemes(:,:).MRI_ROIMarkers =                      [ [ 0.9, 0.9, 0.0 ];... %Color Scheme 1
                                                                           [ 0.9, 0.9, 0.0 ] ];  %Color Scheme 2
        
        oDispOpt = class( oDispOpt, 'IO_DisplayOptions' );
        
	%Constructeur copieur
    elseif( isa(args, 'IO_DisplayOptions') )
        oDispOpt.Indices    = args.Indices;
        oDispOpt.NbItems    = args.NbItems;
        oDispOpt.v_bItemChk = args.v_bItemChk;
        oDispOpt.Lbl_HelmetHoleId_Interval = args.Lbl_HelmetHoleId_Interval;
        oDispOpt = class( oDispOpt, 'IO_DisplayOptions' );
        
    %Constructeur parametrique
    else
        disp( 'IO_DisplayOptions: Constructeur parametrique inexistant' );
    end
    
    