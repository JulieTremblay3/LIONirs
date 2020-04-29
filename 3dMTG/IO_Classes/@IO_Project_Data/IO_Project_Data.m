function obj = IO_Project_Data(args)
%IO_Project_Data Optical imaging project data object.

%Constructeur par défaut
if nargin == 0
    
    %Complete Helmet Raw Digitization data. Probes positions, names and
    %types only.
    obj.m_oDig_CompleteHelmet = Digitization;
    
    %Subject Fiducials Digitization data. Probes Positions, names and
    %types only.
    obj.m_oDig_SubjectFids = Digitization;
    
    %Helmet: Complete Helmet fitted on Subject Fids, montage info, normals,
    %etc. Includes CompleteHelmet digitization registered in SubjectFids
    %space.
    obj.m_Helmet = Helmet;
    
    %MRI Data: Voxels, Surfaces, Segmentation, Markers, CoRegistration
    %matrix
    obj.m_MRI_Data = MRI_Data;
    
    %Display Options: Graphical Elements checked in the Display Option
    %window
    obj.m_oDispOpt = IO_DisplayOptions;
    
%     bDig_CompleteHelmet_Loaded = false;
%     bDig_SubjectFids_Loaded = false;
%     bMRI_Data_Loaded = false;
%     bSeg_Skin_Loaded = false;
%     bSeg_Cortex_Loaded = false;
    
    obj = class(obj,'IO_Project_Data');

%Constructeur copieur
elseif isa(args,'IO_Project_Data')
    obj.m_oDig_CompleteHelmet  = args.m_oDig_CompleteHelmet;
    obj.m_oDig_SubjectFids     = args.m_oDig_SubjectFids;
    obj.m_Helmet               = args.m_Helmet;
    obj.m_MRI_Data             = args.m_MRI_Data;
    obj.m_oDispOpt             = args.m_oDispOpt;
%     bDig_CompleteHelmet_Loaded = args.bDig_CompleteHelmet_Loaded;
%     bDig_SubjectFids_Loaded    = args.bDig_SubjectFids_Loaded;
%     bMRI_Data_Loaded           = args.bMRI_Data_Loaded;
%     bSeg_Skin_Loaded           = args.bSeg_Skin_Loaded;
%     bSeg_Cortex_Loaded         = args.bSeg_Cortex_Loaded;
    obj = class(obj,'IO_Project_Data');

%Constructeur paramétrique - inutilisé pour l'instant.
else
    disp( 'IO_Project_Data Constructor:Error:Unknown parameter' );
end