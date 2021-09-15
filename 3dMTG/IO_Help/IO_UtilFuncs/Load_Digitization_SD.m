function oDig = Load_Digitization_SD(NIRS)

    %Nouveau casque
    oDig = Digitization;
    
    vProbes = get_vProbes( oDig ); 
    sEmptyProbe = vProbes(1);
    
    %[matPolarisElems, matFids] = Read_LocatorElp( FullFileName_ELP );
    
    matFids = [0	-0.1	0
                0.1	0	0
                -0.10 0	0];
    oDig = set_matFiducials( oDig, matFids );
    NbElems = 1;
    matPolarisElems = cell(numel(NIRS.SD.SrcPos)/3+numel(NIRS.SD.DetPos)/3,5);
    %NIRS.SD.SrcPos
    %project unit must be in meter here
    if isfield(NIRS.SD,'SpatialUnit')
        if   strcmp(NIRS.SD.SpatialUnit,'mm')
            convscale = 0.001
        elseif  strcmp(NIRS.SD.SpatialUnit,'cm')
            convscale = 0.01
        else
            convscale = 1
        end
    else
        convscale = 1
    end

    for i = 1:numel(NIRS.SD.SrcPos)/3
     matPolarisElems(NbElems,1) = {400};
     matPolarisElems(NbElems,2) = {['S',num2str(i)]};
     matPolarisElems(NbElems,3) = {NIRS.SD.SrcPos(i,1)*convscale};
     matPolarisElems(NbElems,4) = {NIRS.SD.SrcPos(i,2)*convscale};
     matPolarisElems(NbElems,5) = {NIRS.SD.SrcPos(i,3)*convscale};
     NbElems = NbElems+1;
    end
    
       %NIRS.SD.SrcPos
    for i = 1:numel(NIRS.SD.DetPos)/3
     matPolarisElems(NbElems,1) = {400};
     matPolarisElems(NbElems,2) = {['D',num2str(i)]};
     matPolarisElems(NbElems,3) = {NIRS.SD.DetPos(i,1)*convscale};
     matPolarisElems(NbElems,4) = {NIRS.SD.DetPos(i,2)*convscale};
     matPolarisElems(NbElems,5) = {NIRS.SD.DetPos(i,3)*convscale};
     NbElems = NbElems+1;
    end
    
      
    %Coords de fichier .ESP
    for( iPolarElem=1:size(matPolarisElems,1) )  
        
        CurProbe = sEmptyProbe;
        CurProbe.Type  = matPolarisElems{iPolarElem,1};
        CurProbe.Label = cell2mat(matPolarisElems(iPolarElem,2));
        CurProbe.Coord = [ matPolarisElems{iPolarElem,3}, matPolarisElems{iPolarElem,4}, matPolarisElems{iPolarElem,5}, 1 ];
        
        vProbes( iPolarElem ) = CurProbe;
    end
    
    oDig = set_vProbes( oDig, vProbes );
