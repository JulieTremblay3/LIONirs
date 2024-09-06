function obj = MRI_Data(args)
%MRI_Data class constructor - Voxels, Markers, Segmentation and Surface Data
 
%Constructeur par défaut
if( nargin == 0 )
    
    obj.matVoxels = [];       %matrice de voxels
    obj.matCortexSegmentation = []; %matrice de booleens (taille=voxels)
    obj.matSkinSegmentation = []; %matrice de booleens (taille=voxels)
    obj.matFiducials = [];    % [NAS;LPA;RPA], voxel-space 
    obj.matCoRegistration = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    obj.matTransformManual = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    obj.VoxDim.XColumn = 0.001; %dimension physique 
    obj.VoxDim.YSlice  = 0.001; %d'un voxel,
    obj.VoxDim.ZRow    = 0.001; %(en metre)
    obj.VoxDepth = 8;
    obj.Vox_IMax = 256;
    obj.Vox_IMin = 0;
    obj.PathFile = '';
    obj.vMarkers(1).Label = '';   %Marker Label
    obj.vMarkers(1).Coord = [0,0,0,0];
    obj.Surfaces.CortexHiRes.VertexBuffer = [];
    obj.Surfaces.CortexHiRes.IndexBuffer = [];
    obj.Surfaces.CortexLowRes.VertexBuffer = [];
    obj.Surfaces.CortexLowRes.IndexBuffer = [];
    obj.Surfaces.Skin.VertexBuffer = [];
    obj.Surfaces.Skin.IndexBuffer = [];
    obj.Surfaces.Skin.Vcolor = [];
    obj.Surfaces.CortexHiRes.Vcolor = [];
    obj.Surfaces.CortexLowRes.Vcolor = [];
    obj.sn_file = [];
    obj = class(obj, 'MRI_Data');
    
%Constructeur copieur
elseif( isa(args, 'MRI_Data' ) )
    obj.matVoxels = args.matVoxels;
    obj.matCortexSegmentation = args.matCortexSegmentation;
    obj.matSkinSegmentation = args.matSkinSegmentation;
    obj.matFiducials = args.matFiducials;
    obj.matCoRegistration = args.matCoRegistration;
    obj.matTransformManual = args.matTransformManual
    obj.VoxDim.XColumn = args.VoxDim.XColumn;
    obj.VoxDim.YSlice  = args.VoxDim.YSlice;
    obj.VoxDim.ZRow    = args.VoxDim.ZRow;
    obj.VoxDepth = args.VoxDepth;
    obj.Vox_IMax = args.Vox_IMax;
    obj.Vox_IMin = args.Vox_IMin;
    obj.PathFile = args.PathFile;
    obj.vMarkers = args.vMarkers;
    obj.Surfaces = args.Surfaces;
    obj.sn_file = args.sn_file;
    obj = class(obj, 'MRI_Data');
    
%Constructeur paramétrique (prend un objet nii (Analyse MRI struct) en param)
else
    obj.matVoxels = [];
    obj.matCortexSegmentation = [];
    obj.matSkinSegmentation = [];
    obj.matFiducials = [];
    obj.matCoRegistration = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    obj.matTransformManual = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];

    obj.VoxDim.XColumn = 0.001; %dimension physique 
    obj.VoxDim.YSlice  = 0.001; %d'un voxel,
    obj.VoxDim.ZRow    = 0.001; %(en metre)
    obj.VoxDepth = 8;
    obj.Vox_IMax = 256;
    obj.Vox_IMin = 0;
    obj.PathFile = '';
    obj.vMarkers(1).Label = '';   %Marker Label
    obj.vMarkers(1).Coord = [0,0,0,0];
    obj.Surfaces.CortexHiRes.VertexBuffer = [];
    obj.Surfaces.CortexHiRes.IndexBuffer = [];
    obj.Surfaces.CortexHiRes.Vcolor = [];
    obj.Surfaces.CortexLowRes.VertexBuffer = [];
    obj.Surfaces.CortexLowRes.IndexBuffer = [];
    obj.Surfaces.CortexLowRes.Vcolor = [];
    obj.Surfaces.Skin.VertexBuffer = [];
    obj.Surfaces.Skin.IndexBuffer = [];
    obj.Surfaces.Skin.Vcolor = [];
    obj.sn_file = [];
    if( isfield( args, 'hdr' ) )
        if( isfield( args.hdr, 'dime' ) )
            if( isfield( args.hdr.dime, 'pixdim' ) )
                if( numel( args.hdr.dime.pixdim ) >= 4 )
                    %Copnversion mm - metre (1/1000)
                    obj.VoxDim.XColumn = args.hdr.dime.pixdim(2)/1000;
                    obj.VoxDim.YSlice  = args.hdr.dime.pixdim(3)/1000;
                    obj.VoxDim.ZRow    = args.hdr.dime.pixdim(4)/1000;
                end
            end
            if( isfield( args.hdr.dime, 'bitpix' ) )
                obj.VoxDepth = args.hdr.dime.bitpix;
            end
            if( isfield( args.hdr.dime, 'glmax' ) )
                obj.Vox_IMax = args.hdr.dime.glmax;
            end
            if( isfield( args.hdr.dime, 'glmin' ) )
                obj.Vox_IMin = args.hdr.dime.glmin;
            end
        end
    end
    
    if( isfield( args, 'fileprefix' ) )
        obj.PathFile = args.fileprefix;
    end
    
    if( isfield( args, 'img' ) )
        obj.matVoxels = args.img;
    end
    
    
    obj = class(obj, 'MRI_Data');
end