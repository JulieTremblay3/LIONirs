%Fonction permettant de remettre a jour les valeurs de profondeur de trous
%selon le mode de calcul: Si la profondeur est uniforme pour tous les
%trous, alors l'information du MRI est ignoree. Autrement, la profondeur
%sera calculee individuellement pour chaque trou (par recalage sur la
%segmentation de l'IRM)
function oHelmet = Update_vHolesDepth( oHelmet, oMRI )

    vHoles = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );
    
    %Calcul de profondeurs uniformes
    if( sMtg.Gen_Params.bUseUniformSkinDepth )
        for( p=1:numel( vHoles ) )
            vHoles(p).SkinDepth = sMtg.Gen_Params.UniformSkinDepth;
        end
        oHelmet = set_vHoles( oHelmet, vHoles );
    
    %On ne peut pas calculer les profondeurs individuelles: mettre zero
    elseif( isempty(get_SkinSegmentation(oMRI)) )
        for( p=1:numel( vHoles ) )
            vHoles(p).SkinDepth = 0;
        end
        oHelmet = set_vHoles( oHelmet, vHoles );
    
    %Calcul de profondeurs individuelles selon la segmentation du MRI (si
    %disponible)
    else
        oHelmet = Fit_Helmet_OnMRISegmentations( oHelmet, oMRI, true, true, false );
    end
        
