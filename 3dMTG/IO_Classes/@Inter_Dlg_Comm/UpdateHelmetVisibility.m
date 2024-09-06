%**************************************************************************
% function UpdateHelmetVisibility
%
% Fonction permettant de mettre à jour la matrice de visibilité des trous
% du casque en fonction de la position de la caméra. Les trous de casque
% dont la normale pointe dans même direction que celle de la caméra
% sont flaggés comme non-visibles si le casque n'est pas en mouvement.
%
% Particularité: Tous les trous sont toujours visibles quand le casque
%                est en cours de rotation (IsAxeRotating == true).
%
%**************************************************************************
function Obj = UpdateHelmetVisibility( Obj, oDispOpt, oHelmet )
     
    if( ~ishandle(Obj.HelmetAxeDisp.axe_handle) || ~ishandle(Obj.HelmetAxeDisp.axe_handle) )
        disp('UpdateHelmetVisibility:Error - axe handle uninitialised');
        return; 
    end
     
    CameraPosition = get( Obj.HelmetAxeDisp.axe_handle, 'CameraPosition' );
    vHoles = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );

    
    Obj.HelmetAxeDisp.v_bVisible = zeros(1,numel(vHoles));
	Obj.HelmetAxeDisp.v_bDisplayLabels = zeros(1,numel(vHoles));
        
    for( p = 1:numel(vHoles) )

        Obj.HelmetAxeDisp.v_bVisible(p) = false;
        Obj.HelmetAxeDisp.v_bDisplayLabels(p) = false;

        if( vHoles(p).Type == 400 )   
            vn = [ vHoles(p).Normal.x; vHoles(p).Normal.y; vHoles(p).Normal.z ];

            %Verifier si la camera se trouve du cote de la normale
            if( CameraPosition*vn > 0  || Obj.HelmetAxeDisp.IsAxeRotating || (vn')*vn == 0 || get_DispOptChecked( oDispOpt, 'App_ItemsAlwaysVisible' ))
                Obj.HelmetAxeDisp.v_bVisible(p) = true;
            end %if

            %Gestion de l'intervalle d'affichage des identificateurs
            strLbl = vHoles(p).Label;
            iNumBegin = find( (int8(strLbl) >= int8('0')) & (int8(strLbl) <= int8('9')) );

            if( iNumBegin )
                strNum = strLbl( iNumBegin:numel(strLbl) );
                iNumEnd = find( ( int8(strNum) < int8('0') ) | ( int8(strNum) > int8('9') ) )-1;

                if( iNumEnd )
                    strNum = strNum(1:iNumEnd);
                end

                Num = strread( strNum );

                if( ~mod( Num, get_HoleIdInterval( oDispOpt ) ) )
                    Obj.HelmetAxeDisp.v_bDisplayLabels(p) = true;
                end
            end
        end %if
    end %for
    
    