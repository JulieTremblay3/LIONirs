function oDig = Load_Digitization( FullFileName_ELP )

    %Nouveau casque
    oDig = Digitization;
    
    vProbes = get_vProbes( oDig );
    sEmptyProbe = vProbes(1);
    
    [matPolarisElems, matFids] = Read_LocatorElp( FullFileName_ELP );
    
    oDig = set_matFiducials( oDig, matFids );
    
    %Coords de fichier .ESP
    for( iPolarElem=1:size(matPolarisElems,1) )  
        
        CurProbe = sEmptyProbe;
        CurProbe.Type  = matPolarisElems{iPolarElem,1};
        CurProbe.Label = cell2mat(matPolarisElems(iPolarElem,2));
        CurProbe.Coord = [ matPolarisElems{iPolarElem,3}, matPolarisElems{iPolarElem,4}, matPolarisElems{iPolarElem,5}, 1 ];
        
        vProbes( iPolarElem ) = CurProbe;
    end
    
    oDig = set_vProbes( oDig, vProbes );
