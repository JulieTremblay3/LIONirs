function obj = Helmet(args)
    %HELMET class constructor - Gestion des informations et des matrices de 
    %                           coordonnées du casque.

    %Constructeur par défaut
    if nargin == 0

        %Informations générales du casque
        obj.ID = 0;   % Numéro du casque ex: 00001
        obj.NeighborhoodOk = 0; %Voisinage calcule
        obj.Type = 'Empty'; % 'Empty' / 'Labo IO HSJ' / 'Unknown'
        obj.HoleTypeList.SkinFitTest    = 100;
        obj.HoleTypeList.Correction     = 200;
        obj.HoleTypeList.HelmetRef      = 300;
        obj.HoleTypeList.NormalHole     = 400;
        obj.HoleTypeList.RegionPoint    = 500;
        obj.mat_ilHoles = [];
        obj.v_Center = [0,0,0];  % Centre du casque

        %Vecteur de structs d'information sur chaque trou du casque
        obj.v_Holes(1).Type = 0;    %Elem Type (100,200,300,400, etc)
        obj.v_Holes(1).Label = '';  %Elem Name
        obj.v_Holes(1).Coord.x = 0; % Coordonnées du trou
        obj.v_Holes(1).Coord.y = 0; % 
        obj.v_Holes(1).Coord.z = 0; % 
        obj.v_Holes(1).Normal.x = 0; % Coordonnées de la normale du trou
        obj.v_Holes(1).Normal.y = 0; % 
        obj.v_Holes(1).Normal.z = 0; % 
        obj.v_Holes(1).Transformation = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
        obj.v_Holes(1).SkinDepth = 0; % 
        obj.v_Holes(1).TalSkin.x = 0;
        obj.v_Holes(1).TalSkin.y = 0;
        obj.v_Holes(1).TalSkin.z = 0;
        obj.v_Holes(1).CortexDepth = 0; % 
        obj.v_Holes(1).TalCortex.x = 0;
        obj.v_Holes(1).TalCortex.y = 0;
        obj.v_Holes(1).TalCortex.z = 0;
        obj.v_Holes(1).Neighbors.Nb = 0; % Nombre de voisins <= 'NeighborhoodRange'
        obj.v_Holes(1).Neighbors.v_Near(1).p  = 0; % Indice 1D ds matHoles du voisin
        obj.v_Holes(1).Neighbors.v_Near(1).Dist = 0; % Distance entre ce voisin et le trou de la matrice matHoles.
        obj.v_Holes(1).CanBeDet = true;
        obj.v_Holes(1).CanBeSrc = true;
        obj.v_Holes(1).IsFromCompleteHelmet = false;   %Elem Type
        obj.v_Holes(1).RegionsID = 0; % Flag de presence ds region d'interet

        %Structure d'information sur le montage
        obj.Mtg_Data.mat_RegionNorms = zeros(1,1,3);
        obj.Mtg_Data.mat_RegionBorder = zeros(1,1,3);
        obj.Mtg_Data.matFiducials = []; % Points de fiducie (X+; Y+; Y-)
        obj.Mtg_Data.v_HolesMtg(1) = 0; % 0XX0XX (src1src2) ou XX000000 det1
        obj.Mtg_Data.v_HolesMtg_Selected = []; %JT 1 = selectionne 0 = non selectionne
        obj.Mtg_Data.v_HolesEle = [];%JT cell de label pour les électrodes
        obj.Mtg_Data.v_pSrc = []; % Liste de sources utilisées (indice 1D du tableau matHoles)
        obj.Mtg_Data.v_pDet = []; % Liste de détecteurs utilisés (indice 1D du tableau matHoles)
        obj.Mtg_Data.v_pEle = []; % Liste de electrodes utilisées (indice 1D du tableau matHoles)
        obj.Mtg_Data.Gen_Params.Min_Dist_Mux          = 1/100; % Distance minimale voulue pour un mux.
        obj.Mtg_Data.Gen_Params.Max_Dist_Mux          = 10/100; % Distance maximale voulue pour un mux.
        obj.Mtg_Data.Gen_Params.DistContamination     = 1/100; % Distance maximale de contamination. %changé pour mux 32
        obj.Mtg_Data.Gen_Params.Nb_Sources_Disp       = 128;      % Nombre de sources disponibles.
        obj.Mtg_Data.Gen_Params.Nb_Detect_Disp        = 16;       % Nombre de détecteurs disponibles.
        obj.Mtg_Data.Gen_Params.Nb_Longueur_Onde      = 2;       % Nombre de longueurs d'onde différentes.
        obj.Mtg_Data.Gen_Params.NbBanks = 8;           % Nombre de banques de sources s'allumant en même temps.
        obj.Mtg_Data.Gen_Params.NbSrcPerBank = 16;     % Nombre de sources dans chaque banque de sources.
        obj.Mtg_Data.Gen_Params.UniformSkinDepth = 0;    % Profondeur pour le calcul de la distance(en metres)
        obj.Mtg_Data.Gen_Params.bUseUniformSkinDepth = false;    % Flag: utiliser ou non une profondeur uniforme (voir "obj.v_Holes(1).SkinDepth": profondeurs individuelles)
        obj.Mtg_Data.Gen_Params.MontageAutomatic = 0; %Permet d'activer ou de désactiver le placement automatique des sources
        obj.Mtg_Data.Gen_Params.ElectrodeType = '10-10'; %Type 10-10 system, E128 system Neuroscan, personnal
        obj.Mtg_Data.Gen_Params.AcqSystem = 'ISS'; %ISS Format de donnée couplage Hopital Ste-Justine 018001
        % IMAGINC format de donnée prototype 002001, 004003
        obj = class(obj, 'Helmet');

    %Constructeur copieur
    elseif( isa(args, 'Helmet' ) )
        disp( 'Copie Helmet en cours...' );
        obj.ID                = args.ID;
        obj.NeighborhoodOk    = args.NeighborhoodOk;
        obj.Type              = args.Type;
        obj.HoleTypeList      = args.HoleTypeList;
        obj.mat_ilHoles       = args.mat_ilHoles;
        obj.v_Center          = args.v_Center;
        obj                   = set_matHoles( obj, get_matHoles(args));
        obj                   = set_Mtg( obj, get_Mtg(args));
        obj                   = class(obj, 'Helmet');
    
    
    %Constructeur parametrique
    elseif( isa(args, 'Digitization' ) )
        obj = Helmet();
        sEmptyHole = obj.v_Holes(1);
        
        vProbes = get_vProbes( args );
        for( p=1:numel( vProbes ) )
            obj.v_Holes(p) = sEmptyHole;
            obj.v_Holes(p).Coord.x = vProbes(p).Coord(1);
            obj.v_Holes(p).Coord.y = vProbes(p).Coord(2);
            obj.v_Holes(p).Coord.z = vProbes(p).Coord(3);
            obj.v_Holes(p).Type = vProbes(p).Type;
            obj.v_Holes(p).Label = vProbes(p).Label;
        end
        
        obj.Mtg_Data.matFiducials = get_matFiducials( args );
        obj = ComputeTransform( ComputeHNormals( ComputeNeighborhood( obj ) ) );
    else
        disp( 'Helmet Constructor: Unknown parameter' );
    end