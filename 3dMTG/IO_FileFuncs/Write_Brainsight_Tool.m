function Write_Brainsight_Tool( oHelmet, oDig_CompleteHelmet, bExportSrcDetOnly, PathFileName )

    if( ~PathFileName )
        return;
    end
    
    %Trous du casque
    vHoles = get_vHoles( oHelmet );
    
    %Montage du casque (identificateurs de sources et de detecteurs, etc.)
    sMtg = get_Mtg( oHelmet );
    
    %Determination des elements a ecrire dans le fichier
    if( bExportSrcDetOnly )
        pElems = [sMtg.v_pDet, sMtg.v_pSrc];
    else
        pElems = find([vHoles.Type] == get_PtTypeNo(oHelmet, 'NormalHole'));
    end
    
    %Cellules de Labels des selections
    vCellLabels = {vHoles(pElems).Label};
    
    %Booleens de provenance
    v_bFromCompleteHelmet = [vHoles(pElems).IsFromCompleteHelmet];
    
    %Verifier que la provenance des elements est bel et bien d'un fichier
    %de casque complet (LOCATOR elp).
    if( numel(find(v_bFromCompleteHelmet)) ~= numel(v_bFromCompleteHelmet) )
        warndlg_modal( 'Only elements from Complete Helmet LOCATOR .elp file will be exported.' );
    end
    
    %Construction du casque contenant uniquement les elements en provenance
    %du fichier .elp de casque complet (.elp format LOCATOR).
    %Le constrcuteur parametreique de Helmet calcule automatiquement les
    %normales lors de la constrcution.
    oHelmetLOC = Helmet(oDig_CompleteHelmet);
    
    %Determination des elements LOCATOR a ecrire (les autres sont ignores).
    %Les Labels sont utilises afin de determiner les elements a ecrire.
    vHolesLOC = get_vHoles( oHelmetLOC );
    pElemsLOC = [];
    for( p=1:numel(vHolesLOC) )
        if( ~isempty(find(strcmp( vCellLabels, vHolesLOC(p).Label ))) )
            pElemsLOC(numel(pElemsLOC)+1) = p;
        end
    end
    
    %Ecriture du fichier
    fid_BrSight = fopen( PathFileName, 'w' );
    
    fprintf( fid_BrSight, 'OptodeMesh Start\n' );
    fprintf( fid_BrSight, 'Mesh Name\n' );
    fprintf( fid_BrSight, 'HSJ Export Test\n' );
    fprintf( fid_BrSight, 'Mesh Type\n' );
    fprintf( fid_BrSight, 'Projection Helmet\n' );
    fprintf( fid_BrSight, 'Conform method\n' );
    fprintf( fid_BrSight, 'project along normals\n' );
    
    %Boucle d'ecriture pour chaque element
    for( iElem=1:numel(pElemsLOC) )
        
        p = pElemsLOC(iElem);
        
        fprintf( fid_BrSight, 'Optode Start\n' );

        fprintf( fid_BrSight, 'Optode Property Label\n' );
        fprintf( fid_BrSight, '%s\n', vHolesLOC(p).Label );

        fprintf( fid_BrSight, 'Optode Start Location to Origin Transform\n' );

        if( vHolesLOC(p).Normal.x || vHolesLOC(p).Normal.y || vHolesLOC(p).Normal.z )

            % Vecteurs n et z
            vn = [ vHolesLOC(p).Normal.x; vHolesLOC(p).Normal.y; vHolesLOC(p).Normal.z ];
            vn = vn./(vn'*vn)^0.5;
            vz = [ 0; 0; 1; ];

            % Produit vectoriel de n et z (axe de rotation)
            vABC = cross( vn, vz );

            % Normalisation de l'axe de rotation
            vABC = vABC./(vABC'*vABC)^0.5;

            % Determination de l'angle de rotation
            tetha = 2*atan2( ((vn-vz)'*(vn-vz))^0.5, ((vn+vz)'*(vn+vz))^0.5 );

            % Creation du Quaternion de rotation: Q=(X,Y,Z,W);
            sin_a = sin( tetha / 2 );
            cos_a = cos( tetha / 2 );
            quat_rotation = [ vABC(1) * sin_a, vABC(2) * sin_a, ...
                              vABC(3) * sin_a, cos_a ];

            %Normalisation du quaternion
            quat_rotation = quat_rotation ./ (quat_rotation*quat_rotation')^0.5;

            XX = quat_rotation(1)*quat_rotation(1);
            XY = quat_rotation(1)*quat_rotation(2);
            XZ = quat_rotation(1)*quat_rotation(3);
            XW = quat_rotation(1)*quat_rotation(4);
            YY = quat_rotation(2)*quat_rotation(2);
            YZ = quat_rotation(2)*quat_rotation(3);
            YW = quat_rotation(2)*quat_rotation(4);
            ZZ = quat_rotation(3)*quat_rotation(3);
            ZW = quat_rotation(3)*quat_rotation(4);

            %Matrice de transformation du trou courant
            %Conversion des coordonnees m to mm (*1000)
            T1 = [ 1 - 2 * ( YY + ZZ ), 2 * ( XY - ZW ),     2 * ( XZ + YW ),      0; ...
                   2 * ( XY + ZW ),     1 - 2 * ( XX + ZZ ), 2 * ( YZ - XW ),      0; ...
                   2 * ( XZ - YW ),     2 * ( YZ + XW ),     1 - 2 * ( XX + YY ),  0; ...  
                   vHolesLOC(p).Coord.x*1000,   vHolesLOC(p).Coord.y*1000,   vHolesLOC(p).Coord.z*1000     1  ];


           fprintf( fid_BrSight, '%12.8f %12.8f %12.8f %12.8f\n', T1(1,1), T1(1,2), T1(1,3), T1(1,4) );
           fprintf( fid_BrSight, '%12.8f %12.8f %12.8f %12.8f\n', T1(2,1), T1(2,2), T1(2,3), T1(2,4) );
           fprintf( fid_BrSight, '%12.8f %12.8f %12.8f %12.8f\n', T1(3,1), T1(3,2), T1(3,3), T1(3,4) );
           fprintf( fid_BrSight, '%12.8f %12.8f %12.8f %12.8f\n', T1(4,1), T1(4,2), T1(4,3), T1(4,4) );
        else
           fprintf( fid_BrSight, '%12.8f %12.8f %12.8f %12.8f\n', 1, 0, 0, 0 );
           fprintf( fid_BrSight, '%12.8f %12.8f %12.8f %12.8f\n', 0, 1, 0, 0 );
           fprintf( fid_BrSight, '%12.8f %12.8f %12.8f %12.8f\n', 0, 0, 1, 0 );
           fprintf( fid_BrSight, '%12.8f %12.8f %12.8f %12.8f\n', vHolesLOC(p).Coord.x*1000, vHolesLOC(p).Coord.y*1000, vHolesLOC(p).Coord.z*1000, 1 );
        end

        fprintf( fid_BrSight, 'Optode End\n' );
    end
    
    fprintf( fid_BrSight, 'OptodeMesh End\n' );

    fclose(fid_BrSight);