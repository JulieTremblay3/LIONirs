function varargout = IO_Dlg_MainWindow(varargin)
% IO_DLG_MAINWINDOW M-menuparameters for IO_Dlg_MainWindow.fig

    % Last Modified by GUIDE v2.5 08-Sep-2021 15:50:09

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @IO_Dlg_MainWindow_OpeningFcn, ...
                       'gui_OutputFcn',  @IO_Dlg_MainWindow_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
    
    
%**************************************************************************
%  S E C T I O N   0 0 1 - L A N C E M E N T   D E   L ' I N T E R F A C E
%**************************************************************************


% --------------------------------------------------------------------
% function IO_Dlg_MainWindow_OpeningFcn
%
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IO_Dlg_MainWindow (see VARARGIN)
%
% --- Executes just before IO_Dlg_MainWindow is made visible.
function IO_Dlg_MainWindow_OpeningFcn(hObject, eventdata, handles, varargin)

    % Choose default command line output for IO_Dlg_MainWindow
    handles.output = hObject;

    %Projet vide (Appel au constructeur de la classe IO_Project_Data)
    oProjectData = IO_Project_Data();
    
    %Initialisation de l'objet de communication inter-fenetres (Controlleur)
    %(Appel au constructeur de la classe Inter_Dlg_Comm)
    oInterDlgComm = Inter_Dlg_Comm();
    
    %Initialisation des options d'affichage
    %(Appel au constructeur de la classe IO_DisplayOptions
    oDispOpt = IO_DisplayOptions();
     
    %Casque vide pour l'instant (Contient les coordonnees numerisees)
    oHelmet = get_Helmet( oProjectData );
   
    %Initialisation de l'affichage 3D (Helmet Axe Display)
    oInterDlgComm = InitHelmetAxeDisp( oInterDlgComm, handles.axes_Helmet, oHelmet );

    axes( handles.axes_Helmet );
    
    %Initialisation des options d'affichage
    oInterDlgComm = Refresh_All( oDispOpt, oInterDlgComm, oProjectData );
    
    %Opengl = affichage plus rapide
    set(gcf, 'Renderer', 'opengl');
   
    % Update handles structure
    handles.m_oDisplayOptions = oDispOpt;
    handles.m_oInterDlgComm = oInterDlgComm;
    handles.m_oProjectData = oProjectData;
 
    guidata(hObject, handles);
    %Initialisation des icones du toolbar
    Init_Toolbar(handles);
    
% --------------------------------------------------------------------
function varargout = IO_Dlg_MainWindow_OutputFcn(hObject, eventdata, handles)
    % Get default command line output from handles structure
    varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function MenuFileOpen_Callback(hObject, eventdata, handles)
    file = uigetfile('*.fig');
    if ~isequal(file, 0)
        open(file);
    end

% --------------------------------------------------------------------
function axes_Helmet_ButtonDownFcn(hObject, eventdata, handles)

% --------------------------------------------------------------------
function axes_Helmet_CreateFcn(hObject, eventdata, handles)
    
% --------------------------------------------------------------------
function figure_MainWindow_CreateFcn(hObject, eventdata, handles)

%**************************************************************************
%  S E C T I O N   0 0 2 - B O U T O N S   E T   M E N U
%**************************************************************************

%002MF - MenuFile
% --------------------------------------------------------------------
function MenuFileNewProject_Callback(hObject, eventdata, handles)
    NewProject(handles);

% --------------------------------------------------------------------
function MenuFileOpenProject_Callback(hObject, eventdata, handles)
    OpenProjectFile( handles );

% --------------------------------------------------------------------
function MenuFileSaveProject_Callback(hObject, eventdata, handles)
    SaveProjectFile(handles);
    
% --------------------------------------------------------------------
function MenuFileExport_Callback(hObject, eventdata, handles)
    IO_Dlg_Export( handles.m_oProjectData );
    
% --------------------------------------------------------------------
function MenuFilePrint_Callback(hObject, eventdata, handles)
    %printpreview(handles.figure_MainWindow);
    print -depsc2 -painters -r300 'tmpPrint.eps';
    system(sprintf( '%s\\%s', pwd, 'tmpPrint.eps' ));

% --------------------------------------------------------------------
function MenuFileExit_Callback(hObject, eventdata, handles)
    selection = questdlg(['Close ' get(handles.figure_MainWindow,'Name') '?'],...
                         ['Close ' get(handles.figure_MainWindow,'Name') '...'],...
                         'Yes','No','Yes');
    if strcmp(selection,'No')
        return;
    end

    delete(handles.figure_MainWindow)

%002ME - MenuEdit
% --------------------------------------------------------------------
function MenuEditSelectAll_Callback(hObject, eventdata, handles)
    
    oDispOpt = handles.m_oDisplayOptions;
    oInterDlgComm = handles.m_oInterDlgComm;
    oHelmet = get_Helmet( handles.m_oProjectData );
    vHoles = get_vHoles(oHelmet);
    
    oInterDlgComm = set_HelmetAxeDisp_vSelection( oInterDlgComm, find( [vHoles.Type] == get_PtTypeNo( oHelmet, 'NormalHole' ) ) );
    Refresh_Selection( oDispOpt, oInterDlgComm, oHelmet );
    
    handles.m_oInterDlgComm = oInterDlgComm;
    guidata( handles.figure_MainWindow, handles );
    
% --------------------------------------------------------------------
function MenuEditSelectNone_Callback(hObject, eventdata, handles)
    
    oDispOpt = handles.m_oDisplayOptions;
    oInterDlgComm = handles.m_oInterDlgComm;
    oHelmet = get_Helmet( handles.m_oProjectData );
    
    oInterDlgComm = ClearSelections( oInterDlgComm );
    Refresh_Selection( oDispOpt, oInterDlgComm, oHelmet );
    
    handles.m_oInterDlgComm = oInterDlgComm;
    guidata( handles.figure_MainWindow, handles );

% --------------------------------------------------------------------
function MenuEditSetAsSource_Callback(hObject, eventdata, handles)
    
    oInterDlgComm = handles.m_oInterDlgComm;
    oHelmet = get_Helmet(handles.m_oProjectData);
    
    oHelmet = SetHoleUsage( oInterDlgComm, handles.m_oDisplayOptions, oHelmet, false );
    handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
    oInterDlgComm = Refresh_All( handles.m_oDisplayOptions, oInterDlgComm, handles.m_oProjectData );
    
    handles.m_oInterDlgComm = oInterDlgComm;
    guidata( handles.figure_MainWindow, handles );

% --------------------------------------------------------------------
function MenuEditSetAsDetector_Callback(hObject, eventdata, handles)

    oInterDlgComm = handles.m_oInterDlgComm;
    oHelmet = get_Helmet(handles.m_oProjectData);
    
    oHelmet = SetHoleUsage( oInterDlgComm, handles.m_oDisplayOptions, oHelmet, true );
    handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
    oInterDlgComm = Refresh_All( handles.m_oDisplayOptions, oInterDlgComm, handles.m_oProjectData );
    
    handles.m_oInterDlgComm = oInterDlgComm;
    guidata( handles.figure_MainWindow, handles );
    
% --------------------------------------------------------------------
function oHelmet = SetHoleUsage( oInterDlgComm, oDispOpt, oHelmet, bCanBeDet )
    vHoles = get_vHoles(oHelmet);
    v_Selected = get_HelmetAxeDisp_vSelection( oInterDlgComm );
    
    if( find( v_Selected ) )
        
        v_pItems = v_Selected(find( v_Selected ));
        
        NbItemsToProcess = length(find( v_Selected ));
        
        for( iFind = 1:NbItemsToProcess )
            p = v_pItems(iFind);
            
            if( p && p <= numel(vHoles) )
                vHoles(p).CanBeDet = bCanBeDet;
                vHoles(p).CanBeSrc = true;              
            end
        end

        oHelmet = set_vHoles(oHelmet, vHoles);
    end

% --------------------------------------------------------------------
function MenuEditSendToCompleteHelmet_Callback(hObject, eventdata, handles)
    %Motte de 'get'
    oHelmet = get_Helmet( handles.m_oProjectData );
    oDigSource = get_Dig_SubjectFids( handles.m_oProjectData );
    oDigDestination = get_Dig_CompleteHelmet( handles.m_oProjectData );
    v_pItemsFromSource = get_HelmetAxeDisp_vSelection( handles.m_oInterDlgComm );
    
    %Transfert des selections de la digitalisation partielle vers la
    %digitalisation complete
    [oHelmet, oDigSource, oDigDestination] = TransferDigitizationItem( oHelmet, oDigSource, oDigDestination, v_pItemsFromSource );
    
    %Motte de 'set'
    handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
    handles.m_oProjectData = set_Dig_SubjectFids( handles.m_oProjectData, oDigSource );
    handles.m_oProjectData = set_Dig_CompleteHelmet( handles.m_oProjectData, oDigDestination );
    handles.m_oInterDlgComm = ClearSelections( handles.m_oInterDlgComm );
    guidata( handles.figure_MainWindow, handles );
    
    %Tout reafficher
    Refresh_All(handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
    

% --------------------------------------------------------------------
function MenuEditSendToPartialHelmet_Callback(hObject, eventdata, handles)
    %Motte de 'get'
    oHelmet = get_Helmet( handles.m_oProjectData );
    oDigSource = get_Dig_SubjectFids( handles.m_oProjectData );
    oDigDestination = get_Dig_CompleteHelmet( handles.m_oProjectData );
    v_pItemsFromSource = get_HelmetAxeDisp_vSelection( handles.m_oInterDlgComm );
    
    %Transfert des selections de la digitalisation partielle vers la
    %digitalisation complete
    [oHelmet, oDigSource, oDigDestination] = TransferDigitizationItem( oHelmet, oDigSource, oDigDestination, v_pItemsFromSource );
    
    %Motte de 'set'
    handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
    handles.m_oProjectData = set_Dig_SubjectFids( handles.m_oProjectData, oDigSource );
    handles.m_oProjectData = set_Dig_CompleteHelmet( handles.m_oProjectData, oDigDestination );
    handles.m_oInterDlgComm = ClearSelections( handles.m_oInterDlgComm );
    guidata( handles.figure_MainWindow, handles );
    
    %Tout reafficher
    Refresh_All(handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
    
% --------------------------------------------------------------------
function MenuEditDelete_Callback(hObject, eventdata, handles)
     
    oInterDlgComm = handles.m_oInterDlgComm;
    oHelmet = get_Helmet(handles.m_oProjectData);
    vHoles = get_vHoles(oHelmet);
    sMtg = get_Mtg(oHelmet);
    v_Selected = get_HelmetAxeDisp_vSelection( oInterDlgComm );
    
    EmptyHelmet = Helmet;
    sEmpty_vHoles = get_vHoles(EmptyHelmet);
    sEmpty_Hole = sEmpty_vHoles(1);
    
    if( find( v_Selected ) )
        
        v_pItems = v_Selected(find( v_Selected ));
        
        NbItemsToErase = length(find( v_Selected ));
        
        strQuestion = ' ';
        if( NbItemsToErase == 1 )
            strQuestion = sprintf( 'You are about to erase %s. Are you sure?', vHoles(v_pItems(1)).Label );
        else
            strQuestion = sprintf('You are about to erase %d items. Are you sure?', NbItemsToErase );
        end
        
        strbutton = questdlg( strQuestion, 'Warning', 'Ok', 'Cancel', 'Ok');
        if( strcmp( strbutton, 'Cancel' ) )
            return;
        end
        
        for( Pos = 1:length(v_pItems) )
            p = v_pItems(Pos);
            
            %Mettre a jour le voisinage...
            if( p && p <= numel(vHoles) )
                %vHoles(p) = sEmpty_Hole;
                if( find( sMtg.v_pDet == p ) )
                    ToggleDetector( p, handles );
                end
                if( find( sMtg.v_pSrc == p ) )
                    ToggleSource( p, handles );
                end                
            end
        end
        
        ValidElems = ones( size( vHoles ) );
        ValidElems( v_pItems ) = zeros(size(v_pItems));
        vHoles = vHoles( find(ValidElems) );

        %Mise a jour de la structure
        oHelmet = set_Mtg(oHelmet, sMtg);
        oHelmet = set_vHoles(oHelmet, vHoles);
        oHelmet = Build2DMap(oHelmet);
        oHelmet = ComputeNeighborhood( oHelmet );
        %oHelmet = ComputeHNormals( oHelmet );
        
        oInterDlgComm = UpdateHelmetVisibility( oInterDlgComm, handles.m_oDisplayOptions, oHelmet);
        oInterDlgComm = ClearSelections( oInterDlgComm );
        handles.m_oInterDlgComm = oInterDlgComm;
        handles.m_oProjectData = set_Helmet(handles.m_oProjectData, oHelmet);
        Refresh_All( handles.m_oDisplayOptions, oInterDlgComm, handles.m_oProjectData );
        
        guidata( handles.figure_MainWindow, handles );
    end
    
% --------------------------------------------------------------------
function oInterDlgComm = ClearSelections( oInterDlgComm )
    oInterDlgComm = set_HelmetAxeDisp_vSelection( oInterDlgComm, [] );
    
%002MD - MenuDisplay
% --------------------------------------------------------------------
function MenuDisplayDisplayOptions_Callback(hObject, eventdata, handles)
    
    %Affichage de la fenetre d'options d'affichage
    oDispOpt = IO_Dlg_DisplayOptions( handles.m_oDisplayOptions );
    
    handles.m_oDisplayOptions = oDispOpt;
    
    %sauvegarde des options
    guidata( handles.figure_MainWindow, handles );
    
    %Tout reafficher
    Refresh_All( oDispOpt, handles.m_oInterDlgComm, handles.m_oProjectData )

%002MM - MenuMontage
% --------------------------------------------------------------------
function MenuMontageClear_Callback(hObject, eventdata, handles)
    handles.m_oProjectData = set_Helmet( handles.m_oProjectData, ClearMtg( get_Helmet(handles.m_oProjectData) ) );
    handles.m_oInterDlgComm = Reset_Ruler( handles.m_oInterDlgComm );
    handles.m_oInterDlgComm = Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
    guidata( handles.figure_MainWindow, handles );
    
    
% --------------------------------------------------------------------
function MenuMontageParams_Callback(hObject, eventdata, handles)
    
    [oInterDlgComm, oHelmet] = MontageParams( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
    handles.m_oInterDlgComm = oInterDlgComm;
    handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
    guidata( handles.figure_MainWindow, handles );

%002MH - MenuHelp
% --------------------------------------------------------------------
function MenuHelpIOMtgHelp_Callback(hObject, eventdata, handles)
    [filepath, ~, ~] = fileparts(mfilename('fullpath')); % Get the path absolute path of the IOMtg main window and only keep filepath.
    system(fullfile(filepath, '..', 'IO_Help', 'Help.doc' )); % Construct the absolute path to the help document and open it.

% --------------------------------------------------------------------
function MenuHelpAboutIOMtg_Callback(hObject, eventdata, handles)
    iconData = zeros(5,5);
    iconData(1:5,3) = ones(5,1);
    
    strBmp = 'About.bmp';
    matpix = imread( strBmp );
    h = msgbox([{'Optical Imaging Montage Creation Assistant Software'};{''};{'           Version 1.65     (06/05/2007)'}],'About IOMtg...','custom',matpix);
    waitfor(h);
        

%002TB - Toolbar
% --------------------------------------------------------------------
function Init_Toolbar(handles)
    
    % Add head picture near the coordinate modification buttons
    matpix = imread( 'iMagic_SkinMesh_Icon.bmp' ); % Use imagic skin mesh icon as it is a head.
    set(handles.pushbutton76,'CData',matpix); % Add the image to non functional push button.

    hPushBtns = findobj(handles.figure_MainWindow,'Style','pushbutton');
    hToggleBtns = findobj(handles.figure_MainWindow,'Style','togglebutton');
    hAllBtns = hPushBtns;
    hAllBtns(length(hPushBtns)+1:length(hPushBtns)+length(hToggleBtns)) = hToggleBtns(1:length(hToggleBtns));
    prefix = 'toolbar';
    
    for( iBtn=1:length(hAllBtns) )
        
        BtnTag = get( hAllBtns(iBtn), 'Tag' );
        if( length(BtnTag) >= length(prefix) )
            
            %Verifier si le tag commence par "Toolbar"
            if( strncmp( prefix, BtnTag, length(prefix) ) )
                
                %Chercher l'image dont le nom est le meme que le Tag
                strBmp = sprintf( '%s.bmp',BtnTag );
                matpix = imread( strBmp );
                set( hAllBtns(iBtn), 'CData', matpix );
            end
        end
    end
    
    hBtnCamOrbit = findobj('Tag','toolbar_CamOrbit');
    if get(hBtnCamOrbit,'value' )
        set( hBtnCamOrbit, 'Value', get_MouseModeActive( handles.m_oInterDlgComm, 'CamOrbit' ) );
        
        %set( hBtnCamOrbit, 'Value', get_HelmetAxeItem_Active( handles.m_oInterDlgComm, 'CamOrbit' ) );
    end
    
% --------------------------------------------------------------------
function toolbar_NewProject_Callback(hObject, eventdata, handles)
    NewProject(handles);
    
% --------------------------------------------------------------------
function toolbar_Open_Callback(hObject, eventdata, handles)
    OpenProjectFile( handles );
    
% --------------------------------------------------------------------
function toolbar_Save_Callback(hObject, eventdata, handles)
    SaveProjectFile(handles);

% --------------------------------------------------------------------
function toolbar_Export_Callback(hObject, eventdata, handles)
    IO_Dlg_Export( handles.m_oProjectData );

% --------------------------------------------------------------------
function toolbar_Print_Callback(hObject, eventdata, handles)    
    print -depsc2 -painters -r300 'tmpPrint.eps';
    system(sprintf( '%s\\%s', pwd, 'tmpPrint.eps' ));
    
% --------------------------------------------------------------------
function toolbar_CamOrbit_Callback(hObject, eventdata, handles)
    PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    handles.m_oInterDlgComm = set_MouseMode( handles.m_oInterDlgComm, 'CamOrbit' );
    set( handles.toolbar_CamOrbit, 'Value', 1 );
    guidata( handles.figure_MainWindow, handles );
        
% --------------------------------------------------------------------
function toolbar_Zoom_Callback(hObject, eventdata, handles)
    PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    handles.m_oInterDlgComm = set_MouseMode( handles.m_oInterDlgComm, 'CamZoom' );
    set( handles.toolbar_Zoom, 'Value', 1 );
    guidata( handles.figure_MainWindow, handles );
    
% --------------------------------------------------------------------
function toolbar_Ruler_Callback(hObject, eventdata, handles)
    PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    handles.m_oInterDlgComm = set_MouseMode( handles.m_oInterDlgComm, 'Ruler' );
    set( handles.toolbar_Ruler, 'Value', 1 );
    
    %Réinitialisation des items de la règle
    sHelmetAxeDisp = get_HelmetAxeDisp(handles.m_oInterDlgComm);
    sHelmetAxeDisp.pRulerFirst = 0;
    sHelmetAxeDisp.pRulerLast = 0;
    handles.m_oInterDlgComm = set_HelmetAxeDisp( handles.m_oInterDlgComm, sHelmetAxeDisp );
    
    %Affichage des changements
    Refresh_Ruler( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet( handles.m_oProjectData ) );
    guidata( handles.figure_MainWindow, handles );
    
% --------------------------------------------------------------------
function toolbar_AddSrc_Callback(hObject, eventdata, handles)
    PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    handles.m_oInterDlgComm = set_MouseMode( handles.m_oInterDlgComm, 'AddSrc' );
    guidata( handles.figure_MainWindow, handles );
    set( handles.toolbar_AddSrc, 'Value', 1 );
        
% --------------------------------------------------------------------
function toolbar_AddDet_Callback(hObject, eventdata, handles)
    PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    handles.m_oInterDlgComm = set_MouseMode( handles.m_oInterDlgComm, 'AddDet' );
    guidata( handles.figure_MainWindow, handles );
    set( handles.toolbar_AddDet, 'Value', 1 );

% --------------------------------------------------------------------
function toolbar_AddEle_Callback(hObject, eventdata, handles)
    PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    handles.m_oInterDlgComm = set_MouseMode( handles.m_oInterDlgComm, 'AddEle' );
    guidata( handles.figure_MainWindow, handles );
    set( handles.toolbar_AddEle, 'Value', 1 );

% --------------------------------------------------------------------
function toolbar_ShowInfo_Callback(hObject, eventdata, handles)
    handles.m_oInterDlgComm = PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    handles.m_oInterDlgComm = set_MouseMode( handles.m_oInterDlgComm, 'GetInf' );
    guidata( handles.figure_MainWindow, handles );
    set( handles.toolbar_ShowInfo, 'Value', 1 );
    
% --------------------------------------------------------------------
function oInterDlgComm = PopOutAllButtons( oDispOpt, oInterDlgComm, oHelmet )
    set( findobj('Tag','toolbar_Ruler'),    'Value', 0 );
    set( findobj('Tag','toolbar_Zoom'),     'Value', 0 );
    set( findobj('Tag','toolbar_CamOrbit'), 'Value', 0 );
    set( findobj('Tag','toolbar_ShowInfo'), 'Value', 0 );
    set( findobj('Tag','toolbar_AddSrc'),   'Value', 0 );
    set( findobj('Tag','toolbar_AddDet'),   'Value', 0 );
    set( findobj('Tag','toolbar_AddEle'),   'Value', 0 );
    set( findobj('Tag','toolbar_ShowInfoFid'),   'Value', 0 );
    set( findobj('Tag','toolbar_datacursor'),   'Value', 0 );
    datacursormode off
    Refresh_Ruler( oDispOpt, oInterDlgComm, oHelmet );
    
% --------------------------------------------------------------------
function toolbar_Parameters_Callback(hObject, eventdata, handles)
    [oInterDlgComm, oHelmet] = MontageParams( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
    handles.m_oInterDlgComm = oInterDlgComm;
    handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
    guidata( handles.figure_MainWindow, handles );

% --------------------------------------------------------------------
function toolbar_DisplayOptions_Callback(hObject, eventdata, handles)

    oDispOpt = IO_Dlg_DisplayOptions( handles.m_oDisplayOptions );
    
    handles.m_oDisplayOptions = oDispOpt;
    
    %sauvegarde des options
    guidata( handles.figure_MainWindow, handles );
    
    %Tout reafficher
    Refresh_All( oDispOpt, handles.m_oInterDlgComm, handles.m_oProjectData );
    
% --------------------------------------------------------------------
function toolbar_IRMView_Callback(hObject, eventdata, handles)
    handles.m_oProjectData = IO_Dlg_MRIView( handles.m_oProjectData );
    guidata( handles.figure_MainWindow, handles );
% --------------------------------------------------------------------
function toolbar_Help_Callback(hObject, eventdata, handles)
    MenuHelpIOMtgHelp_Callback(hObject, eventdata, handles); % Call the dropdown menu help callback function since the outcome is the same (avoid code duplication).

%002IN - Inputs
% --------------------------------------------------------------------
function figure_MainWindow_KeyPressFcn(hObject, eventdata, handles)
    handles.m_oInterDlgComm;
    handles.m_oInterDlgComm = set_KeyPressed( handles.m_oInterDlgComm, lower(get( handles.figure_MainWindow, 'CurrentCharacter' )) );

% --------------------------------------------------------------------
function figure_MainWindow_WindowButtonDownFcn(hObject, eventdata, handles)

    oInterDlgComm = handles.m_oInterDlgComm;
	oHelmet = get_Helmet(handles.m_oProjectData);
    oMRI = get_MRI_Data(handles.m_oProjectData);
    oDispOpt = handles.m_oDisplayOptions;
if get( handles.toolbar_datacursor, 'Value' )
    return
end
    %Graphique 3D:
    [pointerI, in] = PositionSouris( handles.figure_MainWindow, handles.panel_HelmetManipSurface );
   % cameratoolbar
     if( in )
        %Si l'input correspond a une commande connue: (Camorbit, Zoom, Ruler, Src, Det, Ele, Inf)
        if( ...%    ~strcmp( get( handles.figure_MainWindow, 'SelectionType' ), 'open') ...   %N'est pas: Open et Camorbit et Zoom
            ...%&&  ~get_MouseModeActive( oInterDlgComm, 'CamOrbit' ) ...                            
            ...%&&  ~get_MouseModeActive( oInterDlgComm, 'CamZoom'  ) ...                    %Mais est quelque chose parmi:
            ... %&& (    strcmp( get( handles.figure_MainWindow, 'SelectionType' ), 'alt' ) ...       %-Sélection
               ( strcmp( get( handles.figure_MainWindow, 'SelectionType' ), 'alt' ) ...        %-Sélection
                 || strcmp( get( handles.figure_MainWindow, 'SelectionType' ), 'extend' ) ...     %-Identification (SHIFT)
                 || get( handles.toolbar_ShowInfo, 'Value' ) ...                                  %-Identification (btn INF)
                 || get( handles.toolbar_Ruler,    'Value' ) ...                                  %-Regle
                 || get( handles.toolbar_AddDet,   'Value' ) ...                                  %-Det
                 || get( handles.toolbar_AddEle,   'Value' ) ...                                  %-Src
                 || get( handles.toolbar_AddSrc,   'Value' )...                                  %-Ele
                 || get(handles.toolbar_ShowInfoFid,'Value') ))  
            
            %Determination de l'indice lineaire ('pointeur') de l'item
            %selectionne dans le casque.
            pClosestFound = get_pItemClicked( oInterDlgComm, oHelmet, get( handles.axes_Helmet,  'CurrentPoint' ) );
            
            %Quelque chose a ete cliqué
            if( pClosestFound )
                if( strcmp( get( handles.figure_MainWindow, 'SelectionType' ), 'alt' ) )
                    oInterDlgComm = ToggleSelection( oInterDlgComm, oHelmet, pClosestFound )
                    Refresh_Selection( oDispOpt, oInterDlgComm, oHelmet );
                    
                elseif( strcmp( get( handles.figure_MainWindow, 'SelectionType' ), 'extend' ) || get( handles.toolbar_ShowInfo, 'Value') )
                    oInterDlgComm = ToggleIdentification( oInterDlgComm, oHelmet, pClosestFound, handles );
                    Refresh_Information( oDispOpt, oInterDlgComm, oHelmet, oMRI);                
                elseif( get(handles.toolbar_AddDet, 'Value' ) )
                    [oInterDlgComm, oHelmet ] = ToggleDetector( oInterDlgComm, oHelmet, pClosestFound, handles.chk_AutoFibDet, handles.chk_AutoFibSrc );
                    Refresh_SrcDetEle( oDispOpt, oInterDlgComm, oHelmet );                    
                elseif( get(handles.toolbar_AddSrc, 'Value' ) )
                    [oInterDlgComm, oHelmet ] = ToggleSource( oInterDlgComm, oHelmet, pClosestFound, handles.chk_AutoFibDet, handles.chk_AutoFibSrc );
                    Refresh_SrcDetEle( oDispOpt, oInterDlgComm, oHelmet );                
                elseif( get(handles.toolbar_AddEle, 'Value' ) )
                    oHelmet = ToggleElectrode( oInterDlgComm, oHelmet, pClosestFound );
                    Refresh_SrcDetEle( oDispOpt, oInterDlgComm, oHelmet );                
                elseif(  get(handles.toolbar_Ruler, 'Value' ) )
                    oInterDlgComm = ToggleRuler( oInterDlgComm, pClosestFound );
                    Refresh_Ruler( oDispOpt, oInterDlgComm, oHelmet );
                elseif  get( handles.toolbar_ShowInfoFid, 'Value') 
                    oInterDlgComm = ToggleIdentification( oInterDlgComm, oHelmet, pClosestFound, handles );
                    Refresh_Information_Fid( oDispOpt, oInterDlgComm, oHelmet, oMRI);
                end
            end
            
        %Affichage complet pour le temps de la rotation
        else
            if( get_MouseModeActive( oInterDlgComm, 'CamOrbit' ) )
                oInterDlgComm = set_HelmetAxeDisp_AxeRotating(oInterDlgComm,true);
            end
            
            %Gestion de souris
            MouseData = get_MouseData(oInterDlgComm);
            MouseData.IsLBtnDown = true;
            oInterDlgComm = set_MouseData(oInterDlgComm, MouseData);
            
            %Mise a jour de la visibiliten (tout afficher durant la
            %rotation)
            oInterDlgComm = Refresh_Visibility( oDispOpt, oInterDlgComm, oHelmet );
        end
         sMtg = get_Mtg(oHelmet);
         set(handles.text_srs_number,'string',num2str(numel(sMtg.v_pSrc)));
         set(handles.text_det_number,'string',num2str(numel(sMtg.v_pDet)));
         set(handles.text_ele_number,'string',num2str(numel(sMtg.v_pEle)));

        %Sauvegarde des changements 
        handles.m_oInterDlgComm = oInterDlgComm;
        handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
        guidata( handles.figure_MainWindow, handles );
    end
   

% --------------------------------------------------------------------
function figure_MainWindow_WindowButtonMotionFcn(hObject, eventdata, handles)
    
    oInterDlgComm = handles.m_oInterDlgComm;
	MouseData = get_MouseData(oInterDlgComm); 
     
    if( MouseData.IsLBtnDown )

        [pointerI, in] = PositionSouris( handles.figure_MainWindow, handles.panel_HelmetManipSurface );

        if( in && ~strcmp( get( handles.figure_MainWindow, 'SelectionType' ), 'extend') ...
               && ~strcmp( get( handles.figure_MainWindow, 'SelectionType' ), 'alt') )
           
            persistent LastRedrawTime;
            persistent LastPoint;
            
            %Initialisation de premier appel
            if( isempty(LastRedrawTime) || isempty(LastPoint) )
                LastRedrawTime = cputime;
                LastPoint = get( handles.figure_MainWindow, 'CurrentPoint' );
            end
            
            CurPoint = get( handles.figure_MainWindow, 'CurrentPoint' );
            
            DeltaX = LastPoint(1,1) - CurPoint(1,1);
            DeltaY = LastPoint(1,2) - CurPoint(1,2);
            
            %Protection de mouvements trop brusques
            if( abs(DeltaX) > 0.05 || abs(DeltaY) > 0.05 )
                DeltaX = 0;
                DeltaY = 0;
            end

            %Zoom ou camorbit selon le mode souris et l'amplitude de mvt
            if( get_MouseModeActive( oInterDlgComm, 'CamOrbit' ) )
                camorbit(DeltaX*350,DeltaY*350);
            elseif( get_MouseModeActive( oInterDlgComm, 'CamZoom' ) )
                camzoom(1.0-DeltaY*5);
            end

            LastPoint = CurPoint;
        end
    end

% --------------------------------------------------------------------
function figure_MainWindow_WindowButtonUpFcn(hObject, eventdata, handles)

    oInterDlgComm = handles.m_oInterDlgComm;
    oHelmet = get_Helmet( handles.m_oProjectData );
    MouseData = get_MouseData(oInterDlgComm);
    MouseData.IsLBtnDown = false;
    oInterDlgComm = set_MouseData( oInterDlgComm, MouseData );
	oInterDlgComm = set_HelmetAxeDisp_AxeRotating( oInterDlgComm, false);
            
    if( get_MouseModeActive( oInterDlgComm, 'CamOrbit' )  )
        oInterDlgComm = Refresh_Visibility( handles.m_oDisplayOptions, oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    end
    
    handles.m_oInterDlgComm = oInterDlgComm;
	guidata( handles.figure_MainWindow, handles );


% --------------------------------------------------------------------
function OpenProjectFile( handles )
    
    [FileName,PathName] = uigetfile({'*.prj;*.mat','IOMtg Project File  (*.prj;*.mat)';'*.*','All Files (*.*)'},'Open...');
    if( ~FileName ) 
        return;
    end
    
    %initialisation
    LoadedStruct = [];
    LoadedStruct = load([PathName,FileName], '-mat');
    set(handles.figure_MainWindow,'name', ['IOMtg -',FileName]);
    handles.m_oProjectData = Load_PrjStruct( LoadedStruct, true );

    handles.m_oInterDlgComm = InitInterfaceData( handles.m_oInterDlgComm, handles.m_oDisplayOptions, get_Helmet(handles.m_oProjectData) );
    Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
        %Ratio de l'IRM
    oMRI = get_MRI_Data(  handles.m_oProjectData ); 
    vDestPhysicalDim = get_PhysicalVoxDim(oMRI);
    current_mat = get_matTransform(oMRI);
    current_size = [current_mat(3,1),current_mat(1,2),current_mat(2,3)];
    Current_ratio = max(current_size./vDestPhysicalDim);
    set(handles.IRM_ratio,'string',num2str(Current_ratio));    
    oHelmet=get_Helmet(    handles.m_oProjectData);
    sMtg = get_Mtg(oHelmet);
    set(handles.text_srs_number,'string',num2str(numel(sMtg.v_pSrc)));
    set(handles.text_det_number,'string',num2str(numel(sMtg.v_pDet)));
    set(handles.text_ele_number,'string',num2str(numel(sMtg.v_pEle)));
    
    guidata( handles.figure_MainWindow, handles );
    
% --------------------------------------------------------------------
function SaveProjectFile(handles)
    
%initialisation
    SaveStruct = Create_PrjStruct( handles.m_oProjectData );
    
    [FileName,PathName] = uiputfile({'*.prj','IOMtg Project File  (*.prj)';'*.mat','MatLab Project File  (*.mat)';'*.*','All Files (*.*)'},'Save As...');
    if( ~FileName )
        return;
    end
    save([PathName,FileName], 'SaveStruct');
    set(handles.figure_MainWindow,'name', ['IOMtg -' FileName])
   
% Lance la fenetre de creation de nouveau projet
% --------------------------------------------------------------------
function NewProject(handles)
    
    %Passage des handles du parent (fenetre courante)
    h = IO_Dlg_NewProject( handles.figure_MainWindow, handles );
    waitfor(h);
    set(handles.IRM_ratio,'string','1')
    handles = guidata( handles.figure_MainWindow );
    
    oInterDlgComm = handles.m_oInterDlgComm;
    oHelmet = get_Helmet( handles.m_oProjectData );

    oInterDlgComm = InitInterfaceData( oInterDlgComm, handles.m_oDisplayOptions, oHelmet );
    Refresh_All( handles.m_oDisplayOptions, oInterDlgComm, handles.m_oProjectData );
    
    handles.m_oInterDlgComm = oInterDlgComm;

    guidata( handles.figure_MainWindow, handles );

% Lance la fenetre de parametres du montage
% --------------------------------------------------------------------
function [oInterDlgComm, oHelmet] = MontageParams( oDispOpt, oInterDlgComm, oPrjData ) 
    %Fenetre de parametres du montage
    oHelmet = IO_Dlg_GenParams( get_Helmet(oPrjData) );
    
    %Mise a jour de la profondeur de chaque positions
    oHelmet = Update_vHolesDepth( oHelmet, get_MRI_Data( oPrjData ) );
    oPrjData = set_Helmet( oPrjData, oHelmet );
    
    %Rafraichissement de l'interface
    oInterDlgComm = Refresh_All( oDispOpt, oInterDlgComm, oPrjData );

% --------------------------------------------------------------------
% Fonction permettant d'ajouter ou d'enlever un 'x' de selection sur un
% trou du casque. Cette fonction n'affiche rien: il faut 
% utiliser Refresh_Selection pour mettre l'ecran a jour.
% --------------------------------------------------------------------
function oInterDlgComm = ToggleSelection( oInterDlgComm, oHelmet, pClosest )
  
    %Recherche de l'indice de stockage
    PosEmpty=0;
    ClickedTwice=false;

    sHelmetAxeDisp = get_HelmetAxeDisp(oInterDlgComm);
    
    %Liste des selections
    v_pSelection = sHelmetAxeDisp.v_pSelection;
    
    %Pour toutes les selections actuellement existantes
    for( Pos=1:length(sHelmetAxeDisp.v_pSelection) )

        %Verifier si la selection est deja presente : si oui
        %l'effacer
        if( sHelmetAxeDisp.v_pSelection(Pos) == pClosest )
            sHelmetAxeDisp.v_pSelection(Pos) = 0;
            ClickedTwice=true;
            break;   
        end

        %Si la position n'a pas encore ete faire et que l'item est vide
        if( ~PosEmpty && ~sHelmetAxeDisp.v_pSelection(Pos) )
            PosEmpty = Pos;
        end
    end

    if( ~ClickedTwice )
        %Si la liste est pleine, ajouter une nouvelle position
        if( ~PosEmpty )
            PosEmpty = length(sHelmetAxeDisp.v_pSelection)+1;
        end

        sHelmetAxeDisp.v_pSelection(PosEmpty) = pClosest;        
    end

    %Sauvegarde des changements
    oInterDlgComm = set_HelmetAxeDisp( oInterDlgComm, sHelmetAxeDisp );

% --------------------------------------------------------------------
% Fonction permettant d'activer ou de desactiver l'identification de 
% trou (le bouton 'INF'). Cette fonction n'affiche rien: il faut 
% utiliser Refresh_Information pour mettre l'ecran a jour.
% --------------------------------------------------------------------
function oInterDlgComm = ToggleIdentification( oInterDlgComm, oHelmet, p, handles )

    sHelmetAxeDisp = get_HelmetAxeDisp( oInterDlgComm );

    if( ~p || ( sHelmetAxeDisp.pIdentification == p ) )
        sHelmetAxeDisp.pIdentification = 0;
    else
        sHelmetAxeDisp.pIdentification = p;   
    end
    
    oInterDlgComm = set_HelmetAxeDisp(oInterDlgComm, sHelmetAxeDisp);

% --------------------------------------------------------------------
% Fonction permettant d'activer ou de desactiver et de determiner les
% points de mesure de la regle. Cette fonction n'affiche rien: il faut 
% utiliser Refresh_Ruler pour mettre l'ecran a jour.
% --------------------------------------------------------------------
function oInterDlgComm = ToggleRuler( oInterDlgComm, p )
    
    sHelmetAxeDisp = get_HelmetAxeDisp(oInterDlgComm);
    
    if( sHelmetAxeDisp.pRulerFirst )
        if( p == sHelmetAxeDisp.pRulerFirst )
            sHelmetAxeDisp.pRulerFirst = 0;
            sHelmetAxeDisp.pRulerLast = 0;
        else
            sHelmetAxeDisp.pRulerLast = p;
        end
    else
        sHelmetAxeDisp.pRulerFirst = p;
    end
    
    oInterDlgComm = set_HelmetAxeDisp(oInterDlgComm, sHelmetAxeDisp);
    
    
% --------------------------------------------------------------------
% Fonction permettant de mettre ou d'enlever un detecteur sur un trou 
% de casque donne (pDet=indice lineaire ds vHoles).
% De plus, cette fonction effectue les appels necessaires a la 
% verification de la contamination et a l'attribution optimise de 
% fibres optiques. Cette fonction n'affiche rien: il faut 
% utiliser Refresh_SrcDetEle pour mettre l'ecran a jour.
% --------------------------------------------------------------------
function [oInterDlgComm, oHelmet ] = ToggleDetector( oInterDlgComm, oHelmet, pDet, bAutoFibDet, bAutoFibSrc )
    
    vHoles = get_vHoles(oHelmet);
    sMtg = get_Mtg(oHelmet);
    MaxDets = sMtg.Gen_Params.Nb_Detect_Disp;
    
    if( sMtg.Gen_Params.Nb_Longueur_Onde <=0  || sMtg.Gen_Params.Nb_Detect_Disp <=0 )
        return;
    end
    
    %Verifier si pDet valide ou si src deja la
    if( ~pDet || ~isempty( find( sMtg.v_pSrc == pDet ) ) )
        return;
    end
    
    %Verifier repetition
    v_iRepeat = find( sMtg.v_pDet == pDet );
    ClickedTwice = ~isempty( v_iRepeat );
    %Retrait du detecteur de la liste de fibres optiques si dblclk.
    sMtg.v_HolesMtg(sMtg.v_pDet(v_iRepeat)) = 0;
    %Retrait du detecteur de la liste de detecteurs si dblclk.
    sMtg.v_pDet(v_iRepeat) = zeros(1, length(v_iRepeat) );
    
    %Si la selection est un detecteur valide 
    if( pDet && ~ClickedTwice && numel(sMtg.v_pDet) < MaxDets )
        if( vHoles(pDet).CanBeDet )
            %Ajout du detecteur dans le montage
            sMtg.v_pDet(numel(sMtg.v_pDet)+1) = pDet;
        else
            return;
        end
    end
    
    %Tri des detecteurs par leur indice 1D
    matTemp = sort( sMtg.v_pDet );
    sMtg.v_pDet( 1:length(find(matTemp))) = matTemp(find(matTemp));
    sMtg.v_pDet( length(find(matTemp))+1:length(sMtg.v_pDet) ) = 0;
    
    
    %Attribution de fibres de detecteur
    if( get( bAutoFibDet, 'Value' ) )
        if sMtg.Gen_Params.MontageAutomatic
        [oHelmetDetSelected, bOk] = Set_Mtg_DefaultDetFibersCfg( set_Mtg(oHelmet,sMtg) );
        else
        [oHelmetDetSelected, bOk] = SetDetectorsFiber_Withoutoptimising(set_Mtg(oHelmet,sMtg),pDet);   
        end
    else
        for( iDet=1:numel(sMtg.v_pDet) )
            pDet = sMtg.v_pDet(iDet);
            if( pDet > numel(sMtg.v_HolesMtg) )
                sMtg.v_HolesMtg(pDet) = 0;
            end
        end
        oHelmetDetSelected = set_Mtg( oHelmet, sMtg );
    end
    
    %Verification de la faisabilite du montage et 
    %attribution des sources en fonction des detecteurs
    [oHelmetSrcAndDetSelected, bOk] = SetSourceFibers( oHelmetDetSelected );
    if( ~bOk )&(sMtg.Gen_Params.MontageAutomatic)    
        h = errordlg( sprintf( 'Contamination Detected: cannot have more \n than %d sources in range of detectors.', get_Mtg_NbSrcPeriods(oHelmet) ) ,'Warning');
        uiwait(h);
        return;
    else
        if( get( bAutoFibSrc, 'Value' ) )
            oHelmet = oHelmetSrcAndDetSelected;
        else
            oHelmet = oHelmetDetSelected;
        end 
    end

% --------------------------------------------------------------------
% Fonction permettant de mettre ou d'enlever une source sur un trou 
% de casque donne (pDet=indice lineaire ds vHoles).
% De plus, cette fonction effectue les appels necessaires a la 
% verification de la contamination et a l'attribution optimise de 
% fibres optiques. Cette fonction n'affiche rien: il faut 
% utiliser Refresh_SrcDetEle pour mettre l'ecran a jour.
% --------------------------------------------------------------------
function [oInterDlgComm,oHelmet] = ToggleSource( oInterDlgComm, oHelmet, pSrc, bAutoFibDet, bAutoFibSrc )

    vHoles = get_vHoles(oHelmet);
    sMtg = get_Mtg(oHelmet);

    %Verifier si : Parametres de montage valides et si p valide ou si det deja la
    if( sMtg.Gen_Params.Nb_Longueur_Onde <=0   || sMtg.Gen_Params.Nb_Sources_Disp <=0 || ...
        sMtg.Gen_Params.NbSrcPerBank <=0 || ~pSrc || ~isempty( find( sMtg.v_pDet == pSrc ) )  )
        h = errordlg(  'Parametres de montage valides ou selection erronee', 'Erreur' );
        waitfor( h );
        return;
    end
    
    %Verifier repetition
    v_iRepeat = find( sMtg.v_pSrc == pSrc );
    ClickedTwice = ~isempty( v_iRepeat );
    sMtg.v_pSrc(v_iRepeat) = zeros(1, length(v_iRepeat) );
    
    %Si la selection est un detecteur valide 
    if( ~ClickedTwice && numel(sMtg.v_pSrc) < get_Mtg_MaxSrcHoles(oHelmet) )
        sMtg.v_pSrc(numel(sMtg.v_pSrc)+1) = pSrc;
    end
     if( get( bAutoFibSrc, 'Value' ) )
        if ~isfield(sMtg.Gen_Params,'MontageAutomatic')
            [oModifiedHelm, bOk] = SetSourceFibers( set_Mtg(oHelmet, sMtg) );
        else
            if (sMtg.Gen_Params.MontageAutomatic)    
                 [oModifiedHelm, bOk] = SetSourceFibers( set_Mtg(oHelmet, sMtg) );
            else
                 [oModifiedHelm, bOk] = SetSourcesFiber_Withoutoptimising( set_Mtg(oHelmet, sMtg),pSrc );
            end    
        end
     else
         oModifiedHelm = set_Mtg(oHelmet, sMtg);
         bOk = true;
     end
    % Si la source ne provoque pas de contamination
    if( bOk )
        oHelmet = oModifiedHelm;
        oInterDlgComm = DisplayNewSrcDist( oInterDlgComm, oHelmet, pSrc, ClickedTwice );
    else
        hDlg = warndlg_modal( 'This source creates contamination and cannot be used' );
        waitfor( hDlg );
    end
    
    %Refresh_All_SrcDet_Related( handles );
    
    
% --------------------------------------------------------------------
% Fonction permettant de mettre ou d'enlever une electrode sur un trou 
% Cette fonction n'affiche rien: il faut 
% utiliser Refresh_SrcDetEle pour mettre l'ecran a jour.
% --------------------------------------------------------------------
function oHelmet = ToggleElectrode( oInterDlgComm, oHelmet, pEle )

    vHoles = get_vHoles(oHelmet);
    sMtg = get_Mtg(oHelmet);
    
    %Verifier repetition
    v_iRepeat = find( sMtg.v_pEle == pEle );
    ClickedTwice = ~isempty( v_iRepeat );
    sMtg.v_pEle(v_iRepeat) = zeros(1, length(v_iRepeat) );
    if numel(sMtg.v_HolesEle) ~= numel(vHoles)
        sMtg.v_HolesEle = cell(1,numel(vHoles));   
    end
     
    %Si la selection est un detecteur valide 
    if( ~ClickedTwice )
        sMtg.v_pEle(numel(sMtg.v_pEle)+1) = pEle;
    end
    [oHelmetnew, bOk] = set_Electrode(set_Mtg(oHelmet, sMtg) ,pEle);
    oHelmet = oHelmetnew;
    
	
% --------------------------------------------------------------------
% Fonction permettant de mettre a jour le 'listbox' de selection de 
% detecteurs a l'ecran.
% --------------------------------------------------------------------
function Update_list_DetUsed( oHelmet )

	sMtg = get_Mtg(oHelmet);
    vHoles = get_vHoles(oHelmet);
    
    hlist_DetUsed =  findobj('Tag','list_DetUsed');
    
    %Mise a jour du listbox de fibres de sources
    strListbox(1,:) = '                 ';
    for( iDet=1:numel(sMtg.v_pDet) )
        p = sMtg.v_pDet(iDet);
        strListbox(iDet,:) = sprintf( '%-7s%-9s ', vHoles(p).Label, get_HoleFiberID( oHelmet, p ) );
    end
    set( hlist_DetUsed, 'String', strListbox );
	set( hlist_DetUsed, 'Value', 1 );
    
    if( get( hlist_DetUsed, 'Value' ) > size( strListbox, 2 ) )
        set( hlist_DetUsed, 'Value', size( strListbox, 2 ) );
    end
    
% --------------------------------------------------------------------
% Fonction permettant de mettre a jour le 'listbox' de selection de 
% sources a l'ecran.
% --------------------------------------------------------------------
function Update_list_SrcUsed( oHelmet )

	sMtg = get_Mtg(oHelmet);
    vHoles = get_vHoles(oHelmet);
    
    hlist_SrcUsed =  findobj('Tag','list_SrcUsed');
    
    %Mise a jour du listbox de fibres de sources
    strListbox(1,:) = '                 ';
    for( iSrc=1:numel(sMtg.v_pSrc) )
        p = sMtg.v_pSrc(iSrc);
        strListbox(iSrc,:) = sprintf( '%-7s%-9s ', vHoles(p).Label, get_HoleFiberID( oHelmet, p ) );
    end
    set( hlist_SrcUsed, 'String', strListbox );
	set( hlist_SrcUsed, 'Value', 1 );
    
    if( get( hlist_SrcUsed, 'Value' ) > size( strListbox, 2 ) )
        set( hlist_SrcUsed, 'Value', size( strListbox, 2 ) );
    end
    % --------------------------------------------------------------------
% Fonction permettant de mettre a jour le 'listbox' de selection de 
% electrode a l'ecran.
% --------------------------------------------------------------------
function Update_list_EleUsed( oHelmet )

	sMtg = get_Mtg(oHelmet);
    vHoles = get_vHoles(oHelmet);
    
    hlist_EleUsed =  findobj('Tag','list_EleUsed');
    
    %Mise a jour du listbox de fibres de sources
    strListbox(1,:) = '                 ';
    for( iEle=1:numel(sMtg.v_pEle) )
        p = sMtg.v_pEle(iEle);
        if isempty(sMtg.v_HolesEle{p}) 
           Label_ele =' '; 
        else
           Label_ele = sMtg.v_HolesEle{p}; 
        end
        strListbox(iEle,:) = sprintf( '%-7s%-9s ', vHoles(p).Label, Label_ele );
    end
    set( hlist_EleUsed, 'String', strListbox );
	set( hlist_EleUsed, 'Value', 1 );
    
    if( get( hlist_EleUsed, 'Value' ) > size( strListbox, 2 ) )
        set( hlist_EleUsed, 'Value', size( strListbox, 2 ) );
    end
    
% --------------------------------------------------------------------
% Fonction intermediaire dans l'algorithme d'attribution de source:
% cette fonction tente une distribution de source classee (default)
% et optimisee si la premiere echoue.
% --------------------------------------------------------------------
function [Helm, bOk] = SetSourceFibers( Helm )

    bOk = true;
    %if( get( handles.chk_AutoFibSrc, 'Value' ) )
        
    %disp( 'Trying DefaultSrcFibersCfg' );
    [Helm, bOk] = Set_Mtg_DefaultSrcFibersCfg( Helm );
    if( ~bOk || get_MtgSrc_Contaminated( Helm ) )
	%if( get_MtgSrc_Contaminated( Helm ) )

        %disp( 'Using Optimised Fibers CFG' );
        [Helm, bOk] = Set_Mtg_OptimisedSrcFibersCfg( Helm );

        %Validation supplementaire
        if( bOk )
            bOk = ~get_MtgSrc_Contaminated( Helm );
        end
    end
    
    
% --------------------------------------------------------------------
% Permet de faire pointer la regle sur une source nouvellement ajoutee.
% --------------------------------------------------------------------
function oInterDlgComm = DisplayNewSrcDist( oInterDlgComm, oHelmet, pNewSrc, bClickedTwice )
    
    sHAD = get_HelmetAxeDisp(oInterDlgComm);
    sMtg = get_Mtg( oHelmet );
    
    %Affichage de la regle:
    if( ~bClickedTwice )
        pClosestDet = 0;
        dClosestDet = sMtg.Gen_Params.DistContamination;
             
        for( j=1:numel(sMtg.v_pDet ) )
            pDet = sMtg.v_pDet(j);

            if( pDet && pNewSrc )
                dist = get_DistHoles( oHelmet, pNewSrc, pDet );

                if( dist < dClosestDet )
                    dClosestDet = dist;
                    pClosestDet = pDet;
                end
            end
        end

        if( pClosestDet )
            %Distance de la derniere source placee du detecteur
            sHAD.pRulerFirst = pClosestDet;
            sHAD.pRulerLast = pNewSrc;            
        end
    else
        %Réinitialisation des items de la règle
        sHAD.pRulerFirst = 0;
        sHAD.pRulerLast  = 0;
    end
    
    oInterDlgComm = set_HelmetAxeDisp( oInterDlgComm, sHAD );
    
% --------------------------------------------------------------------
% Fonction appellee quand un element du listbox de detecteurs est
% selectionne. 
% --------------------------------------------------------------------
function list_DetUsed_Callback(hObject, eventdata, handles)
    
    %Ne pas permettre si l'option n'est pas cochee
    if( get( handles.chk_AutoFibDet, 'Value' ) )
        h = warndlg( 'Automatic Fibers Attribution must be unchecked to manually change Detector fiber' );
        waitfor(h);
        return;
    end
    
    oHelmet = get_Helmet( handles.m_oProjectData );
    vHoles = get_vHoles( oHelmet );
    
    strListbox = get(hObject,'String');
    strItem = strListbox( get(hObject,'Value'),: );
    pItem = Label2p( oHelmet, strItem );
    
    if( pItem )
        oHelmet = IO_Dlg_MtgDetAlpha(oHelmet, pItem);
        handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
        guidata( handles.figure_MainWindow, handles );
        Refresh_SrcDetEle( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    end
    

% --------------------------------------------------------------------
function list_DetUsed_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
% Fonction appellee quand un element du listbox de sources est
% selectionne. 
% --------------------------------------------------------------------
function list_SrcUsed_Callback(hObject, eventdata, handles)
    
    %Ne pas permettre si l'option n'est pas cochee
    if( get( handles.chk_AutoFibSrc, 'Value' ) )
        h = warndlg( 'Automatic Fibers Attribution must be unchecked to manually change Source fiber' );
        waitfor(h);
        return;
    end
    
    oHelmet = get_Helmet( handles.m_oProjectData );
    vHoles = get_vHoles( oHelmet );
    
    strListbox = get(hObject,'String');
    strItem = strListbox( get(hObject,'Value'),: );
    pItem = Label2p( oHelmet, strItem );
    
    if( pItem )
        oHelmet = IO_Dlg_MtgSrcBankNo(oHelmet,pItem);
        handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
        guidata( handles.figure_MainWindow, handles );
        Refresh_SrcDetEle( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    end
    
    
%strListbox = get(hObject,'String');
%disp( sprintf( 'Item Selectionne: %s',strListbox( get(hObject,'Value'),: ) ) );
% --------------------------------------------------------------------
function list_SrcUsed_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
% Remise a zero du montage (retrait de toutes les src, det, ele).
% --------------------------------------------------------------------
function oHelmet = ClearMtg( oHelmet )
    sMtg = get_Mtg(oHelmet); 
    vHoles = get_vHoles( oHelmet );
    
    sMtg.v_pDet = [];
    sMtg.v_pSrc = [];
    sMtg.v_pEle = [];
    sMtg.v_HolesMtg = zeros( size(vHoles) );
    sMtg.v_HolesEle = [];
    oHelmet = set_Mtg(oHelmet,sMtg);
    
% --------------------------------------------------------------------
% Donne la position du pointeur de la souris 
% in = flag
% pointerI(0) = x
% pointerI(1) = y
function [pointerI, in] = PositionSouris(figureH, itemH)    
    in=0;
    resolution= get(0, 'ScreenSize');
    pointerS= get(0, 'PointerLocation')./resolution(3:4);
    posFig= get(figureH, 'Position');
    pointerF= (pointerS - posFig(1:2))./posFig(3:4);
    posItem= get(itemH, 'Position');
    pointerA= (pointerF - posItem(1:2))./posItem(3:4);
    if (sum(pointerA < 0) == 0) & (sum(pointerA > 1) == 0)
        in= 1;
    end
    pointerI = pointerA;


% --- Executes on key press over list_SrcUsed with no controls selected.
function list_SrcUsed_KeyPressFcn(hObject, eventdata, handles)
    
%Initialisation des elements de l'interface: Camera, Controlleur
%d'interface(oInterDlgComm)
function oInterDlgComm = InitInterfaceData( oInterDlgComm, oDispOpt, oHelmet )
  
    InitCamera( oInterDlgComm, oHelmet );
    
    %Initialisation des items du graphique.
    oInterDlgComm = InitHelmetAxeDisp( oInterDlgComm, get_HelmetAxeDisp_hAxe( oInterDlgComm ), oHelmet );
   
    %Mise a jour des elements visibles selon l'angle de la camera.
    oInterDlgComm = UpdateHelmetVisibility( oInterDlgComm, oDispOpt, oHelmet );
    
   
function InitCamera( oInterDlgComm, oHelmet )
    set( get_HelmetAxeDisp_hAxe( oInterDlgComm ), 'CameraTarget', get_Center( oHelmet ) );
    set( get_HelmetAxeDisp_hAxe( oInterDlgComm ), 'CameraPositionMode', 'auto' );
    
    
% Retire tous les detecteurs du montage
function btn_ClearDets_Callback(hObject, eventdata, handles)
    button = questdlg('Are you really want to clear all detectors ?','Delete detectors'); 
    if strcmp(button,'Yes')
        oHelmet = get_Helmet(handles.m_oProjectData);
        sMtg = get_Mtg(oHelmet);
        sMtg.v_pDet = []; %Vider le vecteur de detecteurs

        % Sauvegarde des changements et affichage
        oHelmet = set_Mtg(oHelmet, sMtg);
        handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
        Refresh_SrcDetEle( handles.m_oDisplayOptions, handles.m_oInterDlgComm, oHelmet );
        guidata( handles.figure_MainWindow, handles );
    end

% Retire toutes les sources du montage
function btn_ClearSrcs_Callback(hObject, eventdata, handles)
    
    button = questdlg('Are you really want to clear all sources ?','Delete sources'); 
    if strcmp(button,'Yes')
        oHelmet = get_Helmet(handles.m_oProjectData);
        sMtg = get_Mtg(oHelmet);
        sMtg.v_pSrc = []; %Vider le vecteur de sources

        % Sauvegarde des changements et affichage
        oHelmet = set_Mtg(oHelmet, sMtg);
        handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
        Refresh_SrcDetEle( handles.m_oDisplayOptions, handles.m_oInterDlgComm, oHelmet );
        guidata( handles.figure_MainWindow, handles );
    end
    
% --- Executes when figure_MainWindow is resized.
function figure_MainWindow_ResizeFcn(hObject, eventdata, handles)
    

% --------------------------------------------------------------------
function chk_AutoFibDet_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function chk_AutoFibSrc_Callback(hObject, eventdata, handles)

%Tout redessiner
function oInterDlgComm = Refresh_All( oDispOpt, oInterDlgComm, oPrjData )
    cla; %Clear Axes
    oHelmet = get_Helmet( oPrjData );
    oMRI = get_MRI_Data( oPrjData );
    
    Update_list_DetUsed( oHelmet );
    Update_list_SrcUsed( oHelmet );
    Update_list_EleUsed( oHelmet );
    Update_Background( oDispOpt, oInterDlgComm ); 
    Update_ItemInformationPanel( oInterDlgComm, oHelmet );
    oInterDlgComm = UpdateHelmetVisibility( oInterDlgComm, oDispOpt, oHelmet );
    
    Draw_CoordSysAxes(            oDispOpt, oInterDlgComm );
    Draw_DigFiducials(            oDispOpt, oInterDlgComm, oHelmet );
    Draw_HelmetRef(               oDispOpt, oInterDlgComm, oHelmet  );
    Draw_Lights(                  oDispOpt );
    Draw_MRIFiducials(            oDispOpt, oInterDlgComm, oMRI );
    Draw_MRISurface(              oDispOpt, oInterDlgComm, oMRI );
    Draw_HelmetHoles(             oDispOpt, oInterDlgComm, oHelmet );
    Draw_HoleSurfaceIntersection( oDispOpt, oInterDlgComm, oHelmet );
    Draw_HoleSpherical(           oDispOpt, oInterDlgComm, oHelmet );
    Draw_Optodes(                 oDispOpt, oInterDlgComm, oHelmet );
    Draw_Selection(               oDispOpt, oInterDlgComm, oHelmet );
    Draw_ItemInformation(         oDispOpt, oInterDlgComm, oHelmet );
    Draw_Electrodes(              oDispOpt, oInterDlgComm, oHelmet );
    Draw_Ruler(                   oDispOpt, oInterDlgComm, oHelmet );
    Draw_MuxChannels(             oDispOpt, oInterDlgComm, oHelmet );
    Draw_ROIMarkers(              oDispOpt, oInterDlgComm, oMRI );
    Draw_DigTestMarkers(          oDispOpt, oInterDlgComm, oHelmet );
 
    
%Redessiner seulement les elements qui ont un rapport avec la position de
%la camera
function oInterDlgComm = Refresh_Visibility( oDispOpt, oInterDlgComm, oHelmet )
    oInterDlgComm = UpdateHelmetVisibility( oInterDlgComm, oDispOpt, oHelmet );
    Draw_HelmetHoles( oDispOpt, oInterDlgComm, oHelmet );
    Draw_Optodes(     oDispOpt, oInterDlgComm, oHelmet );
    Draw_Electrodes(  oDispOpt, oInterDlgComm, oHelmet );
    Draw_MuxChannels( oDispOpt, oInterDlgComm, oHelmet );
    Draw_Selection(   oDispOpt, oInterDlgComm, oHelmet );
    Draw_Ruler(       oDispOpt, oInterDlgComm, oHelmet );
    
%Redessiner seulement les elements qui ont un rapport le choix de src det
%ele
function Refresh_SrcDetEle( oDispOpt, oInterDlgComm, oHelmet )
    Draw_Optodes(     oDispOpt, oInterDlgComm, oHelmet );
    Draw_Electrodes(  oDispOpt, oInterDlgComm, oHelmet );
    Draw_MuxChannels( oDispOpt, oInterDlgComm, oHelmet );
    Draw_Ruler(       oDispOpt, oInterDlgComm, oHelmet );
    Update_list_DetUsed( oHelmet );
    Update_list_SrcUsed( oHelmet );
    Update_list_EleUsed( oHelmet);
    
function Refresh_Ruler( oDispOpt, oInterDlgComm, oHelmet )
    Draw_Ruler( oDispOpt, oInterDlgComm, oHelmet );
    
function Refresh_Selection( oDispOpt, oInterDlgComm, oHelmet );
    Draw_Selection( oDispOpt, oInterDlgComm, oHelmet );
    
function Refresh_Information( oDispOpt, oInterDlgComm, oHelmet, oMRI )
    %Affichage du panel d'information
    Update_ItemInformationPanel( oInterDlgComm, oHelmet,oMRI );
    
    %Affichage de la petite fleche d'information
    Draw_ItemInformation( oDispOpt, oInterDlgComm, oHelmet );
    
function Refresh_Information_Fid( oDispOpt, oInterDlgComm, oHelmet, oMRI)
    %Affichage du panel d'information
    Update_ItemInformationPanel_Fid( oInterDlgComm, oHelmet, oMRI) ; 
    %Affichage de la petite fleche d'information
    Draw_ItemInformation( oDispOpt, oInterDlgComm, oHelmet );
    
%Redessiner la case d'information sur les trous
function Update_ItemInformationPanel( oInterDlgComm, oHelmet,oMRI)
%display affiche l'information du trou ou 
	persistent IsAlreadyDisplayed;
    persistent hAssistedModeCtrls;
    
    sHelmetAxeDisp = get_HelmetAxeDisp( oInterDlgComm );

    if( isempty( IsAlreadyDisplayed ) )
        IsAlreadyDisplayed = false;
    end
    
    %Si l'outil d'identification pointe sur un trou
    if( sHelmetAxeDisp.pIdentification )
        
        %Et que le 'Panel' d'identification n'est pas deja affiche
        if( ~IsAlreadyDisplayed )
            %Afficher le panel d'identification  
            %(elements communs: fond blanc, en-tete et zones de textes vides)
           % hAssistedModeCtrls(1) = uicontrol( gcf, 'Tag', 'frame_ItemInfo', 'Style','frame', 'Units', 'normalized', 'Position', [0.819-0.017, 0.720-0.66, 0.176, 0.230], 'BackgroundColor', [1 1 1] );
            hAssistedModeCtrls(2) = uicontrol( gcf, 'Tag', 'static1_ItemInfo', 'Style','listbox', 'Units', 'normalized', 'Position', [0.78, 0.4-0.35, 0.20, 0.350], 'BackgroundColor', [1 1 1], 'FontSize', 14,'HorizontalAlignment','left');
            hAssistedModeCtrls(4) = uicontrol( gcf, 'Tag', 'static_Identification', 'Style','text', 'Units', 'normalized', 'Position', [0.78,0.4, 0.176, 0.030], 'FontSize', 14, 'String', 'Identification','HorizontalAlignment','left');         
            IsAlreadyDisplayed = true;            
        end
        %Afficher l'information specifique (position, type, et autre de trou)
        vHoles = get_vHoles(oHelmet);
        sMtg = get_Mtg(oHelmet);        
        p = sHelmetAxeDisp.pIdentification;
        if( vHoles(p).CanBeDet )
            strItemType = 'Detector ';
        elseif( vHoles(p).CanBeSrc )
            strItemType = 'Source   ';
        else
            strItemType = 'Unknown  ';
        end
         strIDH = sprintf( '%-9s',vHoles(p).Label );     
        SkinFit.x = vHoles(p).Coord.x-vHoles(p).Normal.x*vHoles(p).SkinDepth ;
        SkinFit.y = vHoles(p).Coord.y-vHoles(p).Normal.y*vHoles(p).SkinDepth ;
        SkinFit.z = vHoles(p).Coord.z-vHoles(p).Normal.z*vHoles(p).SkinDepth ;
        CortexFit.x = vHoles(p).Coord.x-vHoles(p).Normal.x*vHoles(p).CortexDepth ;
        CortexFit.y = vHoles(p).Coord.y-vHoles(p).Normal.y*vHoles(p).CortexDepth ;
        CortexFit.z = vHoles(p).Coord.z-vHoles(p).Normal.z*vHoles(p).CortexDepth ;
        if 1
             TalSkin.x = vHoles(p).TalSkin.x;
             TalSkin.y = vHoles(p).TalSkin.y;
             TalSkin.z = vHoles(p).TalSkin.z;
             TalCortex.x = vHoles(p).TalCortex.x;
             TalCortex.y = vHoles(p).TalCortex.y;
             TalCortex.z = vHoles(p).TalCortex.z;
        end
        
        h = findobj('Tag','static1_ItemInfo');       
        set(  h, 'String', [  {sprintf( '%s  %s x,y,z', 'Item Name:',strIDH )};...   
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'Coord: ',vHoles(p).Coord.x,vHoles(p).Coord.y,vHoles(p).Coord.z)};...
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'Norm : ',vHoles(p).Normal.x,vHoles(p).Normal.y,vHoles(p).Normal.z)};... 
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'Skin : ', SkinFit.x, SkinFit.y, SkinFit.z)};... 
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'Tai :', TalSkin.x,TalSkin.y,TalSkin.z)};...
             {sprintf( '%s%0.3f', 'Depth :', vHoles(p).SkinDepth )};...            
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'Cortex : ',CortexFit.x,CortexFit.y,CortexFit.z)};...
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'Tai :', TalCortex.x,TalCortex.y,TalCortex.z)};...
             {sprintf( '%s%0.3f', 'Depth :', vHoles(p).CortexDepth )}]);%;...
  
%              {sprintf( '%s%0.3f,%0.3f,%0.3f', 'IRM : ',matFidsIRM(3,1),matFidsIRM(3,2),matFidsIRM(3,3))};...
%              {sprintf( '%s%0.3f', [vHoles(p).Label ' Distance NAS :'] ,dist)};...
%              {sprintf( '%s%0.3f', [vHoles(p).Label ' Distance IRM :'] ,distMRI)}]);
  


    elseif( ~sHelmetAxeDisp.pIdentification && IsAlreadyDisplayed )
        delete(hAssistedModeCtrls( find( hAssistedModeCtrls .* ishandle(hAssistedModeCtrls)) ));
        IsAlreadyDisplayed = false;
    end
function Update_ItemInformationPanel_Fid( oInterDlgComm, oHelmet,oMRI)
%display affiche l'information du trou ou 
	persistent IsAlreadyDisplayed;
    persistent hAssistedModeCtrls;
    
    sHelmetAxeDisp = get_HelmetAxeDisp( oInterDlgComm );

    if( isempty( IsAlreadyDisplayed ) )
        IsAlreadyDisplayed = false;
    end 
    
    %Si l'outil d'identification pointe sur un trou
    if( sHelmetAxeDisp.pIdentification )
        
        %Et que le 'Panel' d'identification n'est pas deja affiche
        if( ~IsAlreadyDisplayed )
            %Afficher le panel d'identification  
            %(elements communs: fond blanc, en-tete et zones de textes vides)[0.819-0.017
            hAssistedModeCtrls(2) = uicontrol( gcf, 'Tag', 'static1_ItemInfo', 'Style','listbox', 'Units', 'normalized', 'Position', [0.78 0.4-0.35, 0.20, 0.35], 'BackgroundColor', [1 1 1],'FontSize', 14,'HorizontalAlignment','left');
            hAssistedModeCtrls(4) = uicontrol( gcf, 'Tag', 'static_Identification', 'Style','text', 'Units', 'normalized', 'Position', [0.78,0.4, 0.176, 0.030], 'FontSize', 14, 'String', 'Identification','HorizontalAlignment','left' );            
            IsAlreadyDisplayed = true;            
        end
        
        %Afficher l'information specifique (position, type, et autre de trou)
        vHoles = get_vHoles(oHelmet);
        sMtg = get_Mtg(oHelmet);

        p = sHelmetAxeDisp.pIdentification;

        if( vHoles(p).CanBeDet )
            strItemType = 'Detector ';
        elseif( vHoles(p).CanBeSrc )
            strItemType = 'Source   ';
        else
            strItemType = 'Unknown  ';
        end
        matFids = sMtg.matFiducials;
        
        matFidsIRM = get_matFiducials(oMRI);
        %Ajout de coord 4D et Application de la matrice de registration
        matFidsIRM = [matFidsIRM, ones(size(matFidsIRM,1),1)]*get_matTransform(oMRI);
        
        Nas = matFids(1,:);
        dist =  sqrt((vHoles(p).Coord.x  - Nas(1))^2+(vHoles(p).Coord.y  -Nas(2))^2+ (vHoles(p).Coord.z  -Nas(3))^2);
        Nas = matFidsIRM(1,:);
        distMRI =  sqrt((vHoles(p).Coord.x  - Nas(1))^2+(vHoles(p).Coord.y  -Nas(2))^2+ (vHoles(p).Coord.z  -Nas(3))^2);
        h = findobj('Tag','static1_ItemInfo');
        set(h, 'HorizontalAlignment','left')      
         set( h, 'String', [ {sprintf( '%s:  %s ,  %s,  %s', 'FID ','x','y','z')};...
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'NAS : ',matFids(1,1),matFids(1,2),matFids(1,3))};...
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'IRM : ',matFidsIRM(1,1),matFidsIRM(1,2),matFidsIRM(1,3))};... 
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'LPA : ',matFids(2,1),matFids(2,2),matFids(2,3))};...
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'IRM : ',matFidsIRM(2,1),matFidsIRM(2,2),matFidsIRM(2,3))};...
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'RPA : ',matFids(3,1),matFids(3,2),matFids(3,3))};...
             {sprintf( '%s%0.3f,%0.3f,%0.3f', 'IRM : ',matFidsIRM(3,1),matFidsIRM(3,2),matFidsIRM(3,3))};...
             {sprintf( '%s%0.3f', [vHoles(p).Label ' Distance NAS :'] ,dist)};...
             {sprintf( '%s%0.3f', [vHoles(p).Label ' Distance IRM :'] ,distMRI)}]);

    elseif( ~sHelmetAxeDisp.pIdentification && IsAlreadyDisplayed )
        delete(hAssistedModeCtrls( find( hAssistedModeCtrls .* ishandle(hAssistedModeCtrls)) ));
        IsAlreadyDisplayed = false;
    end    
% --------------------------------------------------------------------
% Rafraishissement de la couleur du Background de dessin 
% --------------------------------------------------------------------
function Update_Background( oDispOpt, oInterDlgComm )
    set( get_HelmetAxeDisp_hAxe(oInterDlgComm), 'Color',  get_ItemColor( oDispOpt, 'Background' ) );
    set( get_HelmetAxeDisp_hAxe(oInterDlgComm), 'XColor', get_ItemColor( oDispOpt, 'Background' ) );
    set( get_HelmetAxeDisp_hAxe(oInterDlgComm), 'YColor', get_ItemColor( oDispOpt, 'Background' ) );
    set( get_HelmetAxeDisp_hAxe(oInterDlgComm), 'ZColor', get_ItemColor( oDispOpt, 'Background' ) );
    set( gcf, 'Color', get_ItemColor( oDispOpt, 'Background' ) );
    
% --------------------------------------------------------------------
% Desactivation de la regle (pointe sur rien)
% --------------------------------------------------------------------
function oInterDlgComm = Reset_Ruler( oInterDlgComm )
        
    %Réinitialisation des items de la règle
    sHelmetAxeDisp = get_HelmetAxeDisp(oInterDlgComm);
    sHelmetAxeDisp.pRulerFirst = 0;
    sHelmetAxeDisp.pRulerLast = 0;
    oInterDlgComm = set_HelmetAxeDisp( oInterDlgComm, sHelmetAxeDisp );
    


% --------------------------------------------------------------------
function menu_displaymux16_Callback(hObject, eventdata, handles)
% hObject    handle to menu_displaymux16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prjname = get(handles.figure_MainWindow,'name')
display_mux16(handles.m_oProjectData,prjname)

% --------------------------------------------------------------------
function menu_displaymux32_Callback(hObject, eventdata, handles)
% hObject    handle to menu_displaymux32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prjname = get(handles.figure_MainWindow,'name')
display_mux32(handles.m_oProjectData,prjname)





% --------------------------------------------------------------------
function menu_displaymuxall32_Callback(hObject, eventdata, handles)
% hObject    handle to menu_displaymuxall32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prjname = get(handles.figure_MainWindow,'name')
display_mux32all(handles.m_oProjectData,prjname)


% --- Executes on selection change in list_EleUsed.
function list_EleUsed_Callback(hObject, eventdata, handles)
% hObject    handle to list_EleUsed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_EleUsed contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_EleUsed
    %Ne pas permettre si l'option n'est pas cochee
    
    oHelmet = get_Helmet( handles.m_oProjectData );
    vHoles = get_vHoles( oHelmet );
    
    strListbox = get(hObject,'String');
    strItem = strListbox( get(hObject,'Value'),: );
    pItem = Label2p( oHelmet, strItem );
    
    if( pItem )
        oHelmet = IO_Dlg_MtgElectrode(oHelmet, pItem);
        handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
        guidata( handles.figure_MainWindow, handles );
        Refresh_SrcDetEle( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    end

% --- Executes during object creation, after setting all properties.
function list_EleUsed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_EleUsed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn__clear_ele.
function btn_clear_ele_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_clear_ele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
    button = questdlg('Are you really want to clear all electrodes ?','Delete electrodes'); 
    if strcmp(button,'Yes')
        oHelmet = get_Helmet(handles.m_oProjectData);
        sMtg = get_Mtg(oHelmet);
        vHoles = get_vHoles(oHelmet);
        sMtg.v_pEle = []; %Vider le vecteur de sources
        sMtg.v_HolesEle = cell(1,numel(vHoles));   
        % Sauvegarde des changements et affichage
        oHelmet = set_Mtg(oHelmet, sMtg);
        handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
        Refresh_SrcDetEle( handles.m_oDisplayOptions, handles.m_oInterDlgComm, oHelmet );
        guidata( handles.figure_MainWindow, handles );
    end





% --- Executes on button press in toolbar_ShowInfoFID.
function toolbar_ShowInfoFid_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_ShowInfoFID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toolbar_ShowInfoFID
    handles.m_oInterDlgComm = PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    handles.m_oInterDlgComm = set_MouseMode( handles.m_oInterDlgComm, 'GetInf' );
    guidata( handles.figure_MainWindow, handles );
    set( handles.toolbar_ShowInfoFid, 'Value', 1 );


% --------------------------------------------------------------------
function export_boxy_Callback(hObject, eventdata, handles)
% hObject    handle to export_boxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

export_boxygraphmux32(handles.m_oProjectData)


% --- Executes on button press in toolbar_ShowInfoTai.
function toolbar_ShowInfoTai_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_ShowInfoTai (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toolbar_ShowInfoTai




% --- Executes on button press in btn_Autoattribute_srs_ace_viewfront.
function btn_viewfront_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_viewfront (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Helmet)
oHelmet = get_Helmet( handles.m_oProjectData );
% center = get_Center( ans ) 
set(handles.axes_Helmet,'CameraTarget',[0,0,0]);
set(handles.axes_Helmet,'CameraPosition',[3,0,0]);
set(handles.axes_Helmet,'CameraUpVector',[0,0,1]);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
% set(handles.axes_Helmet,'CameraTarget',center );
% set(handles.axes_Helmet,'CameraPosition',[center(1)+3,center(2),center(3)]);
% set(handles.axes_Helmet,'CameraUpVector',[center(1)+3,center(2),center(3)+1]);

%     oInterDlgComm = handles.m_oInterDlgComm;
%     oHelmet = get_Helmet( handles.m_oProjectData );
%     calc_Center( oHelmet)
%     vHoles = get_vHoles( oHelmet );
%     sMtg = get_Mtg( oHelmet );
%     oInterDlgComm = InitInterfaceData( oInterDlgComm, handles.m_oDisplayOptions, oHelmet );
%     Refresh_All( handles.m_oDisplayOptions, oInterDlgComm, handles.m_oProjectData );

% --- Executes on button press in btn_right_view.
function btn_right_view_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_right_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Helmet)
set(handles.axes_Helmet,'CameraTarget',[0,0,0]);
set(handles.axes_Helmet,'CameraPosition',[0,-3,0]);
set(handles.axes_Helmet,'CameraUpVector',[0,0,1]);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );

% --- Executes on button press in btn_Autoattribute_srs_ace_left_view.
function btn_left_view_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_left_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Helmet)
set(handles.axes_Helmet,'CameraTarget',[0,0,0]);
set(handles.axes_Helmet,'CameraPosition',[0,3,0]);
set(handles.axes_Helmet,'CameraUpVector',[0,0,1]);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );


% --- Executes on button press in btn_Autoattribute_srs_ace_back_view.
function btn_back_view_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_back_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Helmet)
set(handles.axes_Helmet,'CameraTarget',[0,0,0]);
set(handles.axes_Helmet,'CameraPosition',[-3,0,0]);
set(handles.axes_Helmet,'CameraUpVector',[0,0,1]);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );

% --- Executes on button press in btn_Autoattribute_srs_ace_top_view.
function btn_top_view_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_top_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_Helmet)
set(handles.axes_Helmet,'CameraTarget',[0,0,0]);
set(handles.axes_Helmet,'CameraPosition',[0,0,3]);
set(handles.axes_Helmet,'CameraUpVector',[-1,0,0]);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );

% --- Executes on button press in toolbar_Flecherotgauche.
function btn_rotatation_Callback(hObject, eventdata, handles,type)
% hObject    handle to toolbar_Flecherotgauche (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
alpha = str2num(get(handles.edit_alpha,'string'))*pi/180;
delta = str2num(get(handles.edit_delta,'string'))/1000;
oMRI = get_MRI_Data(  handles.m_oProjectData );
dimvoxel = 0.001;
movehole = get(handles.radio_moveHole,'value')
Helm = get_Helmet(handles.m_oProjectData );
sHelmetAxeDisp = get_HelmetAxeDisp( handles.m_oInterDlgComm );
p = sHelmetAxeDisp.pIdentification;
vHoles = get_vHoles(Helm);
T = get_matTransform(oMRI);
if type == 1;
   mRpa = makehgtform('zrotate', -alpha);
   mMan = makehgtform('xrotate', +alpha);
elseif type ==2
   mRpa = makehgtform('xrotate', -alpha);  
   mMan = makehgtform('yrotate', -alpha);  

elseif type == 3;
   mRpa = makehgtform('zrotate', +alpha);
   mMan = makehgtform('xrotate', -alpha); 
elseif type ==4
   mRpa = makehgtform('xrotate', +alpha);
   mMan = makehgtform('yrotate', +alpha); 
elseif type ==5 %-Y
    if movehole        
        vHoles(p).Coord.y = vHoles(p).Coord.y - delta;  
        vHoles(p).Transformation = vHoles(p).Transformation * [1 ,0 ,0 ,0; 0,1,0,0; 0,0,1,0;0,-delta,0,1];
    else
         mRpa = makehgtform('translate',[0 -delta 0]);
         mMan  = makehgtform('translate',[-delta*1/dimvoxel 0 0 ]);

    end    
elseif type ==6 %+z haut
    if movehole
        vHoles(p).Coord.z = vHoles(p).Coord.z + delta;  
         vHoles(p).Transformation = vHoles(p).Transformation * [1 ,0 ,0 ,0; 0,1,0,0; 0,0,1,0;0,0,delta,1];
    else
        mRpa = makehgtform('translate',[0 0 delta]);
        mMan  = makehgtform('translate',[0  delta*1/dimvoxel 0]);
    end
elseif type ==7 %+y
    if movehole
        vHoles(p).Coord.y = vHoles(p).Coord.y + delta;
        vHoles(p).Transformation = vHoles(p).Transformation * [1 ,0 ,0 ,0; 0,1,0,0; 0,0,1,0;0,delta,0,1];
    else
        mRpa = makehgtform('translate',[0 delta 0]);
        mMan  = makehgtform('translate',[delta*1/dimvoxel 0 0 ]);
    end
elseif type ==8 %-z %bas
    if movehole
       vHoles(p).Coord.z = vHoles(p).Coord.z - delta;  
       vHoles(p).Transformation = vHoles(p).Transformation * [1 ,0 ,0 ,0; 0,1,0,0; 0,0,1,0;0,0,-delta,1];
    else
         mRpa = makehgtform('translate',[0 0 -delta]);
         mMan  = makehgtform('translate',[0  -delta*1/dimvoxel 0]);
    end
elseif type ==9  %x %Avant
    if movehole
        vHoles(p).Coord.x = vHoles(p).Coord.x + delta;  
        vHoles(p).Transformation = vHoles(p).Transformation * [1 ,0 ,0 ,0; 0,1,0,0; 0,0,1,0;delta,0,0,1];
    else        
        mRpa = makehgtform('translate',[delta 0 0]);  
        mMan  = makehgtform('translate',[0 0 -delta*1/dimvoxel]);
    end
elseif type == 10 %-x Arrière
    if movehole
        vHoles(p).Coord.x = vHoles(p).Coord.x - delta;  
         vHoles(p).Transformation = vHoles(p).Transformation * [1 ,0 ,0 ,0; 0,1,0,0; 0,0,1,0;-delta,0,0,1];
    else 
        mRpa = makehgtform('translate',[-delta 0 0]);
        mMan  = makehgtform('translate',[0 0 +delta*1/dimvoxel]);
    end
elseif type == 11 
    mRpa = makehgtform('yrotate', -alpha);
    mMan = makehgtform('zrotate', -alpha);
elseif type == 12
    mRpa = makehgtform('yrotate', alpha);
    mMan = makehgtform('zrotate', alpha);
elseif type == 13
    VoxDims = get_PhysicalVoxDim( oMRI );
    ratio = str2num(get(handles.IRM_ratio,'string'));

    vDestPhysicalDim = get_PhysicalVoxDim(oMRI);
    current_mat = get_matTransform(oMRI);
    current_size = [current_mat(3,1),current_mat(1,2),current_mat(2,3)];
    Current_ratio = max(current_size./vDestPhysicalDim);
    ratio = str2num(get(handles.IRM_ratio,'string'));
    
    ratio_apply = 1/Current_ratio*ratio;
    mRpa = [ ratio_apply, 0, 0, 0;...
               0, ratio_apply, 0, 0;...
               0, 0, ratio_apply, 0;...
               0, 0, 0, 1 ]';
    mMan = [ ratio_apply, 0, 0, 0;...
               0, ratio_apply, 0, 0;...
               0, 0, ratio_apply, 0;...
               0, 0, 0, 1 ]';      
end


%Helm = Fit_Helmet_OnMRISegmentations( Helm, oMRI, 1, 1, false );
    if movehole
        Helm = set_vHoles(Helm,vHoles)
    else
        mat_transform = T * mRpa';  %ini tested
        oMRI = set_matTransform(oMRI,mat_transform);
         if type==1|type==2|type==3|type==4|type==5|type==6|type==7|type==8|type==9|type==10|type==11|type==12|type==13
            matTransformManual = get_matTransformManual(oMRI);
            matTransformManual = matTransformManual*mMan;
            oMRI = set_matTransformManual(oMRI, matTransformManual )
         end
%         Helm = Fit_Helmet_OnMRISegmentations_TF( Helm, oMRI, 1, 1, false,mat_transform );
%         %Apply same transform on fiducial also (as the fit on skin export
%         %will be coherent for the next construction !)
%                 handles.NewPrj = set_Helmet( handles.NewPrj, Helm );
%         if size(sMtg.matFiducials,1)==6
%             NAS  =  [sMtg.matFiducials(4,1), sMtg.matFiducials(4,2), sMtg.matFiducials(4,3),1];
%             RPA   =  [sMtg.matFiducials(5,1), sMtg.matFiducials(5,2), sMtg.matFiducials(5,3),1];
%             LPA   =  [sMtg.matFiducials(6,1), sMtg.matFiducials(6,2), sMtg.matFiducials(6,3),1];
%         else
%             NAS  =  [sMtg.matFiducials(1,1), sMtg.matFiducials(1,2), sMtg.matFiducials(1,3),1];
%             RPA   =  [sMtg.matFiducials(2,1), sMtg.matFiducials(2,2), sMtg.matFiducials(2,3),1];
%             LPA   =  [sMtg.matFiducials(3,1), sMtg.matFiducials(3,2), sMtg.matFiducials(3,3),1];
%         end
%         if type==1|type==2|type==3|type==4|type==5|type==6|type==7|type==8|type==9|type==10|type==11|type==12
% %         NAS2 = NAS*mRpa'; 
% %         LPA2 = LPA*mRpa'; 
% %         RPA2 = RPA*mRpa'; 
% %         sMtg.matFiducials(4,1) = NAS2(1);
% %         sMtg.matFiducials(4,2) = NAS2(2); 
% %         sMtg.matFiducials(4,3) = NAS2(3); 
% %         sMtg.matFiducials(5,1) = RPA2(1);
% %         sMtg.matFiducials(5,2) = RPA2(2);
% %         sMtg.matFiducials(5,3) = RPA2(3);
% %         sMtg.matFiducials(6,1) = LPA2(1);
% %         sMtg.matFiducials(6,2) = LPA2(2); 
% %         sMtg.matFiducials(6,3) = LPA2(3); 
%         
%         %ICI recalculé la distance pour la surface de la peau et du cortex
%         %car un déplacement à été effectué. 
%          %check if they have skin 
%         if 1; bSkinSeg = 1; else; bSkinSeg = 0; end
%         if 1; bCortexSeg = 1;else;bCortexSeg=0; end
%           Helm = Fit_Helmet_OnMRISegmentations_manual( Helm, oMRI, 1, 1, false );
%           vHoles = get_vHoles(Helm);
% %            Helm = set_Mtg(Helm, sMtg);
%            Helm = set_vHoles(Helm,vHoles);
%         end
        
 
       
    end

handles.m_oProjectData = set_Helmet( handles.m_oProjectData, Helm );
handles.m_oProjectData = set_MRI_Data(handles.m_oProjectData ,oMRI);
set(hObject,'value',0);

guidata( handles.figure_MainWindow, handles );
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );


% --- Executes on button press in btn_Autoattribute_srs_ace_rotateback.
function btn_Autoattribute_srs_ace_rotateback_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_rotateback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in toolbar_Flechegauche.
function toolbar_Flechegauche_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_Flechegauche (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in toolbar_Flechedroite.
function toolbar_Flechedroite_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_Flechedroite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in toolbar_Flechehaut.
function toolbar_Flechehaut_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_Flechehaut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton54.
function pushbutton54_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in toolbar_Flechebas.
function toolbar_Flechebas_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_Flechebas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in toolbar_Flecherothaut.
function toolbar_Flecherothaut_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_Flecherothaut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in toolbar_Flecherotbas.
function toolbar_Flecherotbas_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_Flecherotbas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in toolbar_Flecheavant.
function toolbar_Flecheavant_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_Flecheavant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in toolbar_Flechearriere.
function toolbar_Flechearriere_Callback(hObject, eventdata, handles)
% hObject    handle to toolbar_Flechearriere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha as a double


% --- Executes during object creation, after setting all properties.
function edit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_delta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_delta as text
%        str2double(get(hObject,'String')) returns contents of edit_delta as a double


% --- Executes during object creation, after setting all properties.
function edit_delta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in data_cursor.
function toolbar_datacursor_Callback(hObject, eventdata, handles)
% hObject    handle to data_cursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    oDispOpt = handles.m_oDisplayOptions;
    if get_DispOptChecked( oDispOpt, 'MRI_SurfaceCortexHiRes' ) 
    handles.m_oInterDlgComm = PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    dcm_obj = datacursormode(handles.figure_MainWindow);
    set(dcm_obj,'UpdateFcn',@findatlaslabelfct)
    set(dcm_obj,'SnapToDataVertex','on') 
    set(dcm_obj,'DisplayStyle','window') 
    set(handles.toolbar_datacursor,'value',1) 
    datacursormode on
    elseif  get_DispOptChecked( oDispOpt, 'MRI_SurfaceSkin' ) 
    handles.m_oInterDlgComm = PopOutAllButtons( handles.m_oDisplayOptions, handles.m_oInterDlgComm, get_Helmet(handles.m_oProjectData) );
    dcm_obj = datacursormode(handles.figure_MainWindow);
    set(dcm_obj,'UpdateFcn',@addaholefct)
    set(dcm_obj,'SnapToDataVertex','on') 
    set(dcm_obj,'DisplayStyle','window') 
    set(handles.toolbar_datacursor,'value',1) 
    datacursormode on  
    end

    

function txt = findatlaslabelfct(empt,event_obj)
pos = get(event_obj,'Position');
cd = get(event_obj.Target,'Vertices');
ind_cursor = find(cd(:,1)==pos(1)& cd(:,2)==pos(2) & cd(:,3)==pos(3)); 
cd_atlas = get(event_obj.Target,'FaceVertexCData');
nb_atlas= max(cd_atlas);
c = cd_atlas(ind_cursor);
if nb_atlas == 251 %69 region
    txt = {['Atlas number : ',num2str(c),' ', findAtlas69(c)]};
elseif nb_atlas >50 & nb_atlas <= 116 
    txt = {['Atlas number : ',num2str(c),' ', findAtlas116(c)]};
% elseif nb_atlas == 90
%     txt={['Atlas number : ',num2str(c),' ', findAtlas90(c)]};
elseif nb_atlas < 50
    txt = {['Atlas number : ',num2str(c),' ', findAtlasBrodmann48(c)]};
end


 function txt = addaholefct(empt,event_obj)
% 
% handles = guihandles(gcf) 
%  oHelmet = get_Helmet( handles.m_oProjectData );
%             vHoles = get_vHoles(oHelmet);
%             oMRI = get_MRI_Data(handles.m_oProjectData);
pos = get(event_obj,'Position');
cd = get(event_obj.Target,'Vertices');
ind_cursor = find(cd(:,1)==pos(1)& cd(:,2)==pos(2) & cd(:,3)==pos(3)); 
txt = {['Coordinate : ',num2str(pos(1)),num2str(pos(2)),num2str(pos(3)) ]};

% --- Executes on button press in MRI_size.
function MRI_size_Callback(hObject, eventdata, handles)
% hObject    handle to MRI_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function IRM_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to IRM_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IRM_ratio as text
%        str2double(get(hObject,'String')) returns contents of IRM_ratio as a double


% --- Executes during object creation, after setting all properties.
function IRM_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IRM_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function Menu_calcul_talairach_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_calcul_talairach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
            oHelmet = get_Helmet( handles.m_oProjectData );
            vHoles = get_vHoles(oHelmet);
            oMRI = get_MRI_Data(handles.m_oProjectData);
            sn_file = get_sn_file( oMRI);
            [name ,path] = uigetfile('.mat');
            sn_file = load([path,name]);
            [HCoordsSkin, HNormsSkin ] = CoRegisterHelmetToMRI( oHelmet, oMRI,2);
            [HCoordsCortex, HNormsCortex ] = CoRegisterHelmetToMRI( oHelmet, oMRI,3);            
            points = [HCoordsSkin(:,1:3);HCoordsCortex(:,1:3)];
            nb_coords = size(HCoordsSkin,1);
            h = waitbar(0.5);
            trans = transform_points( points ,sn_file,[],0,1);
            close(h);
            for p = 1:nb_coords
                vHoles(p).TalSkin.x = trans(p,1);
                vHoles(p).TalSkin.y = trans(p,2) ;
                vHoles(p).TalSkin.z = trans(p,3);
                vHoles(p).TalCortex.x = trans(p+nb_coords,1);
                vHoles(p).TalCortex.y = trans(p+nb_coords,2) ;
                vHoles(p).TalCortex.z = trans(p+nb_coords,3);
            end

oHelmet = set_vHoles(oHelmet, vHoles);
handles.m_oProjectData = set_Helmet(handles.m_oProjectData,oHelmet);
guidata( handles.figure_MainWindow, handles );

function Menu_calcul_projection_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_calcul_talairach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


oHelmet = get_Helmet( handles.m_oProjectData );
vHoles = get_vHoles(oHelmet);
oMRI = get_MRI_Data(handles.m_oProjectData);
oHelmet = Fit_Helmet_OnMRISegmentations( oHelmet, oMRI, 1, 1 ,0  )     
oHelmet = set_vHoles(oHelmet, vHoles);
handles.m_oProjectData = set_Helmet(handles.m_oProjectData,oHelmet);
guidata( handles.figure_MainWindow, handles );


% --------------------------------------------------------------------
function Menu_export_aceg_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_export_aceg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prjname = get(handles.figure_MainWindow,'name')
display_mux32aceg(handles.m_oProjectData,prjname)




% --------------------------------------------------------------------
function menu_review_list_Callback(hObject, eventdata, handles)
% hObject    handle to menu_review_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

oHelmet = get_Helmet( handles.m_oProjectData );
Write_reviewlist(oHelmet)


% --- Executes on button press in btn_Autoattribute_srs_ace_combinaison.
function btn_combinaison_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_combinaison (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oHelmet = get_Helmet( handles.m_oProjectData );
calcul_combinaison_srs(oHelmet)


% --- Executes on button press in btn_Autoattribute_srs_ace_place_det.
function btn_Autoattribute_srs_ace_place_det_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Autoattribute_srs_ace_place_det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oHelmet = get_Helmet( handles.m_oProjectData );
oHelmet = auto_attribution_det(oHelmet)
figure(handles.figure_MainWindow);
handles.m_oProjectData = set_Helmet(handles.m_oProjectData,oHelmet);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
 guidata( handles.figure_MainWindow, handles );



% --- Executes on button press in btn_autoattributionsrs.
function btn_autoattributionsrs_Callback(hObject, eventdata, handles)
% hObject    handle to btn_autoattributionsrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oHelmet = get_Helmet( handles.m_oProjectData );
[not_all,msg, oHelmet] = auto_attribution_srs4level(oHelmet);
figure(handles.figure_MainWindow);
handles.m_oProjectData = set_Helmet(handles.m_oProjectData,oHelmet);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
guidata( handles.figure_MainWindow, handles );
       if not_all == 0
           msgbox(['All sources placed without contamination',msg])
       else
           msgbox(['We are not able to place all srs',msg])
       end









% --- Executes on button press in btn_change_mux.
function btn_change_mux_Callback(hObject, eventdata, handles)
% hObject    handle to btn_change_mux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oHelmet = get_Helmet( handles.m_oProjectData );
oHelmet = inverse_mux(oHelmet);
figure(handles.figure_MainWindow);
handles.m_oProjectData = set_Helmet(handles.m_oProjectData,oHelmet);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
guidata( handles.figure_MainWindow, handles );



% --- Executes on button press in btn_changesrs.
function btn_changesrs_Callback(hObject, eventdata, handles)
% hObject    handle to btn_changesrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oHelmet = get_Helmet( handles.m_oProjectData );
oHelmet = inverse_srs(oHelmet);
figure(handles.figure_MainWindow);
handles.m_oProjectData = set_Helmet(handles.m_oProjectData,oHelmet);
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );
guidata( handles.figure_MainWindow, handles );



% --- Executes on button press in radio_moveHole.
function radio_moveHole_Callback(hObject, eventdata, handles)
% hObject    handle to radio_moveHole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_moveHole
1
if get(handles.radio_moveHole,'value')
    set(handles.radiobutton5,'value',0)
end


% --------------------------------------------------------------------
function menu_CFG_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_CFG_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


SaveStruct = Create_PrjStruct( handles.m_oProjectData );
[FileName,PathName] = uiputfile({'*.cfg','Imaginc Coupling File (*.cfg)'},'Save As...')
if( ~FileName )
return;
end

PathSave = strcat(PathName,FileName);
writeCouplingFile( PathSave, SaveStruct);
set(handles.figure_MainWindow,'name', ['IOMtg -' FileName])


% --- Executes on button press in btn_fit.
function btn_fit_Callback(hObject, eventdata, handles)
% hObject    handle to btn_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Helm = get_Helmet(handles.m_oProjectData );
oMRI = get_MRI_Data(handles.m_oProjectData );
vHoles = get_vHoles(Helm);
matTransform = get_matTransform(oMRI);
matTransformManual = get_matTransformManual(oMRI)

if 1; bSkinSeg = 1; else; bSkinSeg = 0; end
if 1; bCortexSeg = 1;else;bCortexSeg=0; end
Helm = Fit_Helmet_OnMRISegmentations_manual( Helm, oMRI, 1, 1, false );
vHoles = get_vHoles(Helm);
Helm = set_vHoles(Helm,vHoles);
handles.m_oProjectData = set_Helmet( handles.m_oProjectData, Helm );
handles.m_oProjectData = set_MRI_Data(handles.m_oProjectData ,oMRI);
set(hObject,'value',0);

guidata( handles.figure_MainWindow, handles );
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );


% --------------------------------------------------------------------
function menu_CreateTopoLayout_Callback(hObject, eventdata, handles)
% hObject    handle to menu_CreateTopoLayout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Helm = get_Helmet(handles.m_oProjectData );
vHoles = get_vHoles(Helm);
 sMtg = get_Mtg(Helm);
   
listsrs = [018001,020003,022005,024007,026009,028011,030013,032015,017002,019004,021006,023008,025010,027012,029014,031016,...
    050033,052035, 054037, 056039, 058041, 060043, 062045, 064047, 049034, 051036, 053038, 055040, 057042, 059044,061046, 063048,...
    82065,84067,86069,88071,90073,92075,94077,96079,81066,83068,85070,87072,89074,91076,93078,95080,114097,116099,...
    118101,120103,122105,124107,126109,128111,113098,115100,117102,119104,121106,123108,125110,127112]';
listdet = [1:32]*1000000; %maximal num of detector
%dim 1 posx
%dim 2 pozy
%dim 3 posz
%dim 4 distance
id  = 1;
for isrs=1:numel(listsrs)
    psrs = find(listsrs(isrs)==sMtg.v_HolesMtg);
    for idet = 1:numel(listdet)
        pdet = find(listdet(idet)==sMtg.v_HolesMtg);
        if numel(pdet)==1 & numel(psrs)==1
            DetPos= [vHoles(pdet).Coord.x,vHoles(pdet).Coord.y,vHoles(pdet).Coord.z];
            SrcPos= [vHoles(psrs).Coord.x,vHoles(psrs).Coord.y,vHoles(psrs).Coord.z];
             pos(id,1) = (DetPos(1)-SrcPos(1))/2+SrcPos(1);
             pos(id,2) = (DetPos(2)-SrcPos(2))/2+SrcPos(2);
             pos(id,3) = (DetPos(3)-SrcPos(3))/2+SrcPos(3);
             pos(id,4) = sqrt( (DetPos(1)-SrcPos(1))^2 + (DetPos(2)-SrcPos(2))^2+ (DetPos(3)-SrcPos(3))^2);           %'distance entre les canaux
             pos(id,5) = isrs;
             pos(id,6) = idet;                        
             id = id + 1;
        else 
            
        end
       
        
    end
end


idbad =  find( pos(:,4)>0.04);
pos(idbad,:)=[];

X = pos(:,1);
Y = pos(:,2);
Z = pos(:,3);
[THETA, PHI, R] = cart2sph(X,Y,Z);

THETA = THETA *180/pi; %azimuth longitude Axe xy
PHI = PHI *180/pi; %elevation angle latitude
%the trick is to used azimuth projection to get coordinate in 2d
%transfert cartesian to spheric and get the azimuth projection

[THETA, PHI, R] = cart2sph(X,Y,Z);
p = pi/2-PHI;
delta = THETA;
newx = p.*sin(delta);
newy = p.*cos(delta);
% figure;hold on
% for id=1:numel(newx)
%     plot(newx(id) ,newy(id),'x','displayname', ['E',num2str(pos(id,5)),' D',num2str( pos(id,6))])
% end


%Look Left and Right side and go down row by row finding local minima
matrowcol = zeros(numel(newx),2);

leftid = find((newx>0));                 % look left side
[val , leftorder] =  sort(newx(leftid)); % sort center to out
minIndexesL = imregionalmin( newy(leftid(leftorder))); %Local minimal in y (will define effective columne for each minima
idminL = find( minIndexesL );

% %Plot the left in order
% figure; plot(newx(leftid(leftorder)), newy(leftid(leftorder)));hold on
% plot(newx(leftid(leftorder(idminL))), newy(leftid(leftorder(idminL))),'+')

%write the new col and row id in the topo layout
col = 1;
row = 1;
for  id= 1:numel(leftid)
    if minIndexesL(id)==1        
        row = 1;
        col = col+1;
    end
    matrowcol(leftid(leftorder(id)),1) = row;
    matrowcol(leftid(leftorder(id)),2) = col;
    row = row+1;
end



rightid = find((newx<=0));
[val , rightorder] =  sort(newx(rightid),'descend'); % sort center to out
minIndexesR = imregionalmin( newy(rightid(rightorder)));
idminr = find( minIndexesR );

col = numel(idminr)+1;
row = 1;
for  id= 1:numel(rightid)
    if minIndexesR(id)==1        
        row = 1;
        col = col-1;
    end
    matrowcol(rightid(rightorder(id)),1) = row;
    matrowcol(rightid(rightorder(id)),2) = col;
    row = row+1;
end

%View left order
% matrowcol(leftid(leftorder),:);
% matrowcol(rightid(rightorder),:);
%Plot the right in order
% figure; plot(newx(rightid(rightorder)), newy(rightid(rightorder)));hold on;
% plot(newx(rightid(rightorder(idminr))), newy(rightid(rightorder(idminr))),'+');
% for id=1:numel(rightorder)
% plot(newx(rightid(rightorder(id))), newy(rightid(rightorder(id))),'+');
% end

%shif left to end side
matrowcol(leftid(leftorder),2) = matrowcol(leftid(leftorder),2)+ max(matrowcol(rightid(rightorder),2));
matrowcol(leftid(leftorder),:);
matrowcol(rightid(rightorder),:);   
topomat = zeros(max(matrowcol(:,1)),max(matrowcol(:,2)));

for id=1:size(pos,1)
   topomat (matrowcol(id,1),matrowcol(id,2)) = pos(id,5)*100+pos(id,6);
end 
[file,path]=uiputfile( 'new.tpl');
fid = fopen([path,file],'w');
 for i=1:size(topomat,1);        
    for j=1:size(topomat,2)
        fprintf(fid,'%04.3f\t',topomat(i,j));
    end
     fprintf(fid,'\n');
 end

fclose(fid);

% --- Executes on button press in pushbutton76.
function pushbutton76_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
if get(handles.radiobutton5,'value')
    set(handles.radio_moveHole,'value',0)
end



function editinfoholes_Callback(hObject, eventdata, handles)
% hObject    handle to editinfoholes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editinfoholes as text
%        str2double(get(hObject,'String')) returns contents of editinfoholes as a double


% --- Executes during object creation, after setting all properties.
function editinfoholes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editinfoholes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_CreateHoles.
function btn_CreateHoles_Callback(hObject, eventdata, handles)
% hObject    handle to btn_CreateHoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj=datacursormode(handles.figure_MainWindow);
c_info = getCursorInfo(dcm_obj);
if isempty(c_info)
    msgbox('First select a location using Data Cursor tool')
    return
end

oHelmet = get_Helmet( handles.m_oProjectData );
vHoles = get_vHoles(oHelmet);
nbholes = numel(vHoles);

vHoles(nbholes+1).Type = 400;
name =  Zonename;
vHoles(nbholes+1).Label = name;
vHoles(nbholes+1).Coord.x = c_info.Position(1);
vHoles(nbholes+1).Coord.y = c_info.Position(2);
vHoles(nbholes+1).Coord.z = c_info.Position(3);
NORM = c_info.Position./sqrt(c_info.Position(1).^2+c_info.Position(2).^2+c_info.Position(3).^2);
Coord =  c_info.Position+NORM*0.005;
vHoles(nbholes+1).Coord.x = Coord(1);
vHoles(nbholes+1).Coord.y = Coord(2);
vHoles(nbholes+1).Coord.z = Coord(3);

vHoles(nbholes+1).Normal.x = NORM(1);
vHoles(nbholes+1).Normal.y = NORM(2);
vHoles(nbholes+1).Normal.z = NORM(3);
vn = [NORM(1); NORM(2); NORM(3) ];
vz = [ 0; 0; 1; ];
T = AlignementVecteurs( vz, vn );
T(4,1) =  Coord(1);
T(4,2) =  Coord(2);
T(4,3) =  Coord(3);

vHoles(nbholes+1).Transformation = T;
vHoles(nbholes+1).SkinDepth = 0;
vHoles(nbholes+1).TalSkin.x = 0;
vHoles(nbholes+1).TalSkin.y = 0;
vHoles(nbholes+1).TalSkin.z = 0;
vHoles(nbholes+1).CortexDepth = 0;
vHoles(nbholes+1).TalCortex.x = 0;
vHoles(nbholes+1).TalCortex.y = 0;
vHoles(nbholes+1).TalCortex.z = 0;
vHoles(nbholes+1).Neighbors.Nb = vHoles(end-1).Neighbors.Nb;
vHoles(nbholes+1).Neighbors.v_Near = vHoles(end-1).Neighbors.v_Near;
vHoles(nbholes+1).CanBeDet = 1;
vHoles(nbholes+1).CanBeSrs = 1;
vHoles(nbholes+1).IsFromCompleteHelmet = 1;
vHoles(nbholes+1).RegionID = 0; 

sHelmetAxeDisp = get_HelmetAxeDisp(handles.m_oInterDlgComm);
sHelmetAxeDisp.v_bVisible(nbholes+1) = 1;
sHelmetAxeDisp.v_bDisplayLabels(nbholes+1) = 1;
handles.m_oInterDlgComm =  set_HelmetAxeDisp(handles.m_oInterDlgComm, sHelmetAxeDisp);
oHelmet = set_vHoles(oHelmet,vHoles);
handles.m_oProjectData = set_Helmet( handles.m_oProjectData, oHelmet );
guidata( handles.figure_MainWindow, handles );
Refresh_All( handles.m_oDisplayOptions, handles.m_oInterDlgComm, handles.m_oProjectData );


% --------------------------------------------------------------------
function MenuFileCreateproject_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileCreateproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] =uigetfile({'*.nirs','*.nirs NIRS hoMER data file with 2d coordinate ';...
    '*.snirf','*.snirf data file with 2d coordinate ';},'Select file to create a fast 2d project');
%Construct new project structure
newPrj=Createproject_fromsnirf(fullfile(PathName,FileName));
oHelmet = get_Helmet( newPrj);

%oInterDlgComm = handles.m_oInterDlgComm;
handles.m_oProjectData = newPrj;
oInterDlgComm = InitInterfaceData( handles.m_oInterDlgComm, handles.m_oDisplayOptions, oHelmet);

Refresh_All( handles.m_oDisplayOptions,handles.m_oInterDlgComm, handles.m_oProjectData );

guidata( handles.figure_MainWindow, handles );

