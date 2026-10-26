function oDig = Load_Digitization_SD(NIRS,vararging)
% Description:  add the information of the digitalization (coordinate) include in .nirs
% data file the filed SD.SrcPos and SD.DetPos already include de coordinate
% these coordinates are usually in 2d and will be used to create the vHoles
% field of the anatomy project (.prj).
% could include the type varaging
%

    %Nouveau casque
    oDig = Digitization;
    if numel(vararging)>0
        type = vararging(1);
    else
        type =2; %polaris by default
    end
    
    vProbes = get_vProbes( oDig ); 
    sEmptyProbe = vProbes(1);
    
    %[matPolarisElems, matFids] = Read_LocatorElp( FullFileName_ELP );
    disp('Set default fiducial point')
    if type==2 %polaris 
              matFids = [0	-0.1	0 %nasion en m 10cm
                 0.1	0	  0      %LPA en m 10cm
                -0.10  0	 0 ];         
    elseif type==1 %NIRS x 
           matFids = [01	 0	0 %nasion en m 10cm
                0	0.1	  0      %LPA en m 10cm
                0  -0.1	  0];         
    end
    oDig = set_matFiducials( oDig, matFids );

    NbElems = 1;
    matPolarisElems = cell(numel(NIRS.SD.SrcPos)/3+numel(NIRS.SD.DetPos)/3,5);
    %NIRS.SD.SrcPos
    %project unit must be in meter here
    if isfield(NIRS.SD,'SpatialUnit')
        if   strcmp(NIRS.SD.SpatialUnit,'mm')
            convscale = 0.001;
        elseif  strcmp(NIRS.SD.SpatialUnit,'cm')
            convscale = 0.01;
        else
            convscale = 1;
        end
    else
        convscale = 1;
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
   
        if  type ==1 %NIRx defaults nirs file coordinate
              CurProbe.Coord = [-matPolarisElems{iPolarElem,5},  matPolarisElems{iPolarElem,3}, matPolarisElems{iPolarElem,4}, 1 ];      
        else type==2 %Default axes polaris
            CurProbe.Coord = [ matPolarisElems{iPolarElem,3}, matPolarisElems{iPolarElem,4}, matPolarisElems{iPolarElem,5}, 1 ];
        end
        vProbes( iPolarElem ) = CurProbe;
    end

    oDig = set_vProbes( oDig, vProbes );
