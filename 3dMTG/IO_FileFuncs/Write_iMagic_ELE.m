%Exportation des positions de sources et de Detecteurs vers le format .ELE
%de iMagic. Le fichier cree peut ensuite etre ouvert dans un projet iMagic.
function Write_iMagic_ELE( oHelmet, oMRI, PathFileName,projection  )

    if( ~PathFileName )
        return;
    end

    vHoles = get_vHoles( oHelmet );
    sMtg = get_Mtg( oHelmet );
    
    %pElems = vecteur d'indices lineaires correspondants aux src/dets
    %utilises
    pElems = [ sMtg.v_pDet(find(sMtg.v_pDet)), sMtg.v_pSrc(find(sMtg.v_pSrc)) ];
    
    %Fonction de transfert des coordonnees de trous en voxel-space du MRI:
    %HCoords sont les coordonnees du casque en voxel-space (Helmet Coords)
    %(HNorms inutilise dans cette fonction)
    [HCoords, HNorms ] = CoRegisterHelmetToMRI( oHelmet, oMRI,projection );
    
    %Transfert des positions correspondants aux sources et aux detecteurs
    mat_Positions = HCoords(pElems,1:3);
    
    %Construction du vecteur de Labels de source et de detecteurs
    v_Labels(1,:) = '          ';
    for( iElem=1:numel(pElems) )
       v_Labels(iElem,:) = sprintf( '%-10s', cell2mat(strread(get_HoleFiberID(oHelmet,pElems(iElem)), '%s', 1)) );
    end
    
    fid_ele = fopen( PathFileName, 'w' );
    
    for( iLine=1:size(mat_Positions,1) )
        stuff = strread( v_Labels(iLine,:), '%s' );
        fprintf(fid_ele, '%s %4.0f %4.0f %4.0f\r\n', stuff{1}, mat_Positions(iLine,1), mat_Positions(iLine,2), mat_Positions(iLine,3) );
        %normal
        if 0
           ind = pElems(iLine);
           fprintf(fid_ele, '%s %1.4f %1.4f %1.4f\r\n', stuff{1}, vHoles(ind).Normal.x, vHoles(ind).Normal.y, vHoles(ind).Normal.z );
        end
    end
  
        

    fclose(fid_ele);