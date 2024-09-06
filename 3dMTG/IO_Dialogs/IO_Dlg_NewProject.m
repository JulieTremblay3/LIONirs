function varargout = IO_Dlg_NewProject(varargin)
    % IO_DLG_NEWPROJECT M-file for IO_Dlg_NewProject.fig
    %      IO_DLG_NEWPROJECT, by itself, creates a new IO_DLG_NEWPROJECT or raises the existing
    %      singleton*.
    %
    %      H = IO_DLG_NEWPROJECT returns the handle to a new IO_DLG_NEWPROJECT or the handle to
    %      the existing singleton*.
    %
    %      IO_DLG_NEWPROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in IO_DLG_NEWPROJECT.M with the given input arguments.
    %
    %      IO_DLG_NEWPROJECT('Property','Value',...) creates a new
    %      IO_DLG_NEWPROJECT or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before IO_Dlg_NewProject_OpeningFunction gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to IO_Dlg_NewProject_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help IO_Dlg_NewProject

    % Last Modified by GUIDE v2.5 19-Apr-2017 11:47:21

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @IO_Dlg_NewProject_OpeningFcn, ...
                       'gui_OutputFcn',  @IO_Dlg_NewProject_OutputFcn, ...
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


%--------------------------------------------------------------------------
% --- Executes just before IO_Dlg_NewProject is made visible.
function IO_Dlg_NewProject_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to IO_Dlg_NewProject (see VARARGIN)

%     % Choose default command line output for IO_Dlg_NewProject
   handles.output = hObject;
% 
    %Initialisation du dernier repertoire accede
    handles.LastDirectory = pwd;
    handles.load = 0;
    %Handles de la fenetre principale
    handles.hParentFigure = varargin{1};
    handles.hParentHandles = varargin{2};
    
    % Update handles structure
   guidata(hObject, handles);
    
    axes(handles.axe_ELP_Helmet_Icon);
    matpix = imread( 'ELP_Helmet_Icon.bmp' );
    image( matpix  );
    axis image;
	axis off;
    
%     axes(handles.axe_ELP_Subject_Icon);
%     matpix = imread( 'ELP_Subject_Icon.bmp' );
%     image( matpix  ); 
%     axis image;
% 	axis off;
    
    axes(handles.axe_Mtg_Icon);
    matpix = imread( 'Export_MTG_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_iMagicVOX_Icon);
    matpix = imread( 'iMagic_VOX_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_iMagic_CortexMeshHi_Icon);
    matpix = imread( 'iMagic_CortexMeshHi_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_iMagic_CortexMeshLow_Icon);
    matpix = imread( 'iMagic_CortexMeshLow_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_iMagic_SkinMesh_Icon);
    matpix = imread( 'iMagic_SkinMesh_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_iMagic_CortexSEG_Icon);
    matpix = imread( 'iMagic_CortexSEG_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_iMagic_SkinSEG_Icon);
    matpix = imread( 'iMagic_SkinSEG_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_iMagic_ROIMarkers_Icon);
    matpix = imread( 'iMagic_ROIMarkers_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axes_AtlasIRM_Icon);
    matpix = imread( 'AtlasIRM_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    %uiwait(handles.btn_Ok);
    % UIWAIT makes IO_Dlg_NewProject wait for user response (see UIRESUME)
    % uiwait(handles.btn_Ok);


%--------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = IO_Dlg_NewProject_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


% --- Executes on button press in btn_Cancel.
function btn_Cancel_Callback(hObject, eventdata, handles)
    close;

%--------------------------------------------------------------------------
function btn_Ok_Callback(hObject, eventdata, handles)
    
    %Appel du constructeur de IO_Project_Data (Projet vide)
    newPrj = IO_Project_Data;

    %Ajout de la digitalisation de casque complet dans le projet
    if( isfield( handles, 'oDig_CompleteHelmet' ) )
        newPrj = set_Dig_CompleteHelmet( newPrj, handles.oDig_CompleteHelmet );
    end

    %Ajout de la digitalisation partielle (fiducies) dans le projet
    if( isfield( handles, 'oDig_SubjectFiducials' ) )
        newPrj = set_Dig_SubjectFids( newPrj, handles.oDig_SubjectFiducials );
    end
    
    %Coregistration des deux digitalisation vers un seul objet 'Helmet' et
    %enregistrement dans le projet.
    newPrj = set_Helmet( newPrj, FitHelmetOnSubject( get_Dig_SubjectFids(newPrj), get_Dig_CompleteHelmet(newPrj) ) );

        
        Helm = get_Helmet( newPrj );
     
    %Memorisation du centre du casque (servira eventuellement a faire
    %pointer la camera vers le centre du casque)
    calc_Center( get_Helmet(newPrj));
    newPrj = set_Helmet( newPrj, calc_Center( get_Helmet(newPrj) ) );

    %Utilisation d'un montage existant pour créer le nouveau projet
    if isfield(handles, 'oDig_MTGimport')
        helmet1 = get_Helmet(newPrj); 
        helmet2 = handles.oDig_MTGimport.oHelmet;
        helmet3 = combinehelmet(helmet1,helmet2);
        newPrj = set_Helmet( newPrj,  helmet3 );
    end
    
    handles.NewPrj = newPrj;
    
    %Necessaire a toutes les autres importations iMagic
    if( isfield( handles, 'bpp' )  && isfield( handles, 'XVoxDim' ) && isfield( handles, 'YVoxDim' ) && ...
        isfield( handles, 'VOXFids' ) && isfield( handles, 'ZVoxDim' ) && isfield( handles, 'matVoxels' ) )
        Helm = get_Helmet( handles.NewPrj );
        
        %Les dimensions ne sont pas tres intuitives (X,Z,Y)...
        [nVoxX, nVoxZ, nVoxY ] = size(handles.matVoxels);
        
        %Initialisation
        if( handles.bpp == 1 )
            sMRI.img = zeros( nVoxX, nVoxY, nVoxZ, 'uint8' );
        else
            sMRI.img = zeros( nVoxX, nVoxY, nVoxZ, 'uint16' );
        end
        
        %Formattage de la matrice de voxels ->(X,Y,Z): C'est plus intuitif
        disp('Reshaping Voxels');
        for( iY=1:nVoxY )
            sMRI.img(:,iY,:) = reshape( handles.matVoxels(:,:,iY), nVoxX, 1, nVoxZ );
        end
        disp('Done');
        
        %Sauvegarde des informations                  
        sMRI.hdr.dime.bitpix = handles.bpp;
        sMRI.hdr.dime.pixdim = [ 0, handles.XVoxDim, handles.YVoxDim, handles.ZVoxDim ];
%         sMRI.hdr.dime.pixdim = [0,0.5,0.5,0.5];
        sMRI.hdr.dime.glmax = max(max(max(handles.matVoxels)));
        
        %Constructeur objet MRI_Data (prend la struct en parametre) 
        oMRI = MRI_Data( sMRI );
        
        %Sauvegarde des fiducies: handles.VOXFids = [NAS; LPA; RPA]
        oMRI = set_matFiducials( oMRI, handles.VOXFids(1:3,:) );
        
        %Liberation de memoire: 2 mat voxels en memoire = trop!
        handles.matVoxels = [];
        guidata(handles.figure_NewProject,handles);
        
        bSkinSeg = false;
        bCortexSeg = false;
        
        %Sauvegarde de la segmentation
        if( isfield( handles, 'matCortexSegmentation' ) )
            
            [nVoxCortexSegX, nVoxCortexSegZ, nVoxCortexSegY ] = size(handles.matCortexSegmentation);
            
            %Validation de la taille de la segmentation du cortex
            if( nVoxCortexSegX == nVoxX && nVoxCortexSegZ == nVoxZ && nVoxCortexSegY <= nVoxY )
            
                %Uniformisation des tailles
                if( nVoxCortexSegY < nVoxY )
                    handles.matCortexSegmentation(:,:,nVoxCortexSegY+1:nVoxY) = zeros(nVoxX,nVoxZ,nVoxY-nVoxCortexSegY, 'uint8');                    
                end
                
                %Formattage de la matrice de segmentation:
                %->(X,Y,Z): C'est plus intuitif
                tmpSeg = zeros( nVoxX, nVoxY, nVoxZ, 'uint8' );
                disp('Reshaping Segmentation');
                for( iY=1:nVoxY )
                    tmpSeg(:,iY,:) = reshape( handles.matCortexSegmentation(:,:,iY), nVoxX, 1, nVoxZ );
                end
                disp('Done');
                
                oMRI = set_CortexSegmentation( oMRI, tmpSeg );
                bCortexSeg = true;
                
            else
                warndlg( 'La segmentation n''à pu etre chargée' );
                warndlg( 'Les dimensions XZ du fichier de segmentation de cortex ne correspondent pas aux dimensions du fichier de données' );
            end
        end
        
        %Sauvegarde de la segmentation de la peau
        if( isfield( handles, 'matSkinSegmentation' ) )
            
            [nVoxSkinSegX, nVoxSkinSegZ, nVoxSkinSegY ] = size(handles.matSkinSegmentation);
            
            %Validation de la taille de la segmentation
            if( nVoxSkinSegX == nVoxX && nVoxSkinSegZ == nVoxZ && nVoxSkinSegY <= nVoxY )
            
                %Uniformisation des tailles
                if( nVoxSkinSegY < nVoxY )
                    handles.matSkinSegmentation(:,:,nVoxSkinSegY+1:nVoxY) = zeros(nVoxX,nVoxZ,nVoxY-nVoxSkinSegY, 'uint8');                    
                end
                
                %Formattage de la matrice de segmentation:
                %->(X,Y,Z): C'est plus intuitif
                tmpSeg = zeros( nVoxX, nVoxY, nVoxZ, 'uint8' );
                disp('Reshaping Segmentation');
                for( iY=1:nVoxY )
                    tmpSeg(:,iY,:) = reshape( handles.matSkinSegmentation(:,:,iY), nVoxX, 1, nVoxZ );
                end
                disp('Done');
                
                oMRI = set_SkinSegmentation( oMRI, tmpSeg );
                bSkinSeg = true;
                
            else
                h = warndlg( 'La segmentation n''à pu etre chargée' );
                waitfor(h);
                h = warndlg( 'Les dimensions XZ du fichier de segmentation de peau ne correspondent pas aux dimensions du fichier de données' );
                waitfor(h);
            end
        end

        %Etablissement de la matrice de coregistration
        oMRI = CoRegisterToHelmet( oMRI, Helm ); 
         
        %Recalage des trous sur les segmentation (suivant les normales)
        Helm = Fit_Helmet_OnMRISegmentations( Helm, oMRI, bSkinSeg, bCortexSeg, false );
        handles.NewPrj = set_Helmet( handles.NewPrj, Helm ); 
            
        %Sauvegarde des maillages
        if(    ( isfield( handles, 'CortexMeshHiVertex' )  && isfield( handles, 'CortexMeshHiIndex' ) ) ...
            || ( isfield( handles, 'CortexMeshLowVertex' ) && isfield( handles, 'CortexMeshLowIndex' ) ) ...
            || ( isfield( handles, 'SkinMeshVertex' )      && isfield( handles, 'SkinMeshIndex' ) ) )
         
          
            
            %Sauvegarde du maillage haute resolution (++polygones) du cortex
            if( isfield( handles, 'CortexMeshHiVertex' ) && isfield( handles, 'CortexMeshHiIndex' ) )
                oMRI = set_CortexMeshHiRes( oMRI, handles.CortexMeshHiVertex, handles.CortexMeshHiIndex  );
                oMRI = set_CortexHiResVcolor( oMRI, handles.CortexMeshHiVcolor);
            end

            %Sauvegarde du maillage basse resolution (--polygones) du cortex
            if( isfield( handles, 'CortexMeshLowVertex' ) && isfield( handles, 'CortexMeshLowIndex' ) )
                oMRI = set_CortexMeshLowRes( oMRI, handles.CortexMeshLowVertex, handles.CortexMeshLowIndex  );
                oMRI = set_CortexLowResVcolor( oMRI, handles.CortexMeshLowVcolor);              
            end

            %Sauvegarde du maillage de la peau
            if( isfield( handles, 'SkinMeshVertex' ) && isfield( handles, 'SkinMeshIndex' ) )
                oMRI = set_SkinMesh( oMRI, handles.SkinMeshVertex, handles.SkinMeshIndex  );
                oMRI = set_SkinVcolor( oMRI, handles.SkinVcolor);
            end
        end
        %Couleur de l'Atlas                 
        if ~isempty( get(handles.text_atlas,'string'))
            if ~strcmp(get(handles.text_atlas,'string'),' ')
            oMRI  = coloratlas(oMRI,get(handles.text_atlas,'string'));
            end
        end
        
        if( ( isfield( handles, 'vMarkers' ) ) )
            if ~isempty(get(handles.vMarkers,'string')) 
                oMRI = set_vMarkers( oMRI, handles.vMarkers );
            end
        end
        
        if (isfield(handles,'sn_file'))
            if ~isempty(get(handles.sn_file,'string'))
            SN = load(handles.sn_file);
            oMRI = set_sn_file( oMRI, SN);
            [HCoords, HNorms ] = CoRegisterHelmetToMRI( Helm, oMRI,3);
            points = HCoords(:,1:3);            
            trans = transform_points(points,handles.sn_file,[],0,1);
            vHoles = get_vHoles( Helm);
            for p = 1:numel(points(:,1))
                vHoles(p).TalCortex.x = trans(p,1);
                vHoles(p).TalCortex.y = trans(p,2);
                vHoles(p).TalCortex.z = trans(p,3);
            end
             Helm = set_vHoles(Helm,vHoles);
             handles.NewPrj  = set_Helmet( newPrj,  Helm );
            end
        end
        
        if (isfield(handles,'edit_transformationfit'))  
            if ~isempty(get(handles.edit_transformationfit,'string')) %( ~isempty(get(handles.edit_transformationfit,'string' )) )
             load(get(handles.edit_transformationfit,'string'))             
             oMRI = set_matTransformManual(oMRI, IRMAdjustment.matTransformManual );
             oMRI = set_matTransform(oMRI,IRMAdjustment.matTransform);
            end
        end 
        %Sauvegarde du MRI_Data dans le projet en cours
        handles.NewPrj = set_MRI_Data( handles.NewPrj, oMRI );
        
      
    else
        h = warndlg( 'Aucun fichier .VOX specifié. Les surfaces et les segmentations ne peuvent etre chargés' );
        waitfor(h);
    end
     btn_save_path_Callback(hObject, eventdata, handles,0)
    
        
    %Sauvegarde des donnes du projet dans les handles de la fenetre
    %principale
    handles.hParentHandles.m_oProjectData = handles.NewPrj;

    %Saucegarde des donnees dui graphic user interface
    guidata( handles.hParentFigure, handles.hParentHandles );
    
    close(handles.figure_NewProject);
    
    
%**************************************************************************
% SECTION III - GESTION DES BOITES DE TEXTE
%**************************************************************************


%--------------------------------------------------------------------------
function text_iMagicVOXfile_Callback(hObject, eventdata, handles)


%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function text_iMagicVOXfile_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to text_iMagicVOXfile (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%--------------------------------------------------------------------------
function text_SEG_Skin_Callback(hObject, eventdata, handles)


%--------------------------------------------------------------------------
function text_SEG_Skin_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to text_SEG_Skin (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%--------------------------------------------------------------------------
function text_iMagicCortexMeshHi_Callback(hObject, eventdata, handles)
    % Hints: get(hObject,'String') returns contents of text_iMagicCortexMeshHi as text
    %        str2double(get(hObject,'String')) returns contents of
    %        text_iMagicCortexMeshHi as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function text_iMagicCortexMeshHi_CreateFcn(hObject, eventdata, handles)
    % called
    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%--------------------------------------------------------------------------
function text_MRK_Callback(hObject, eventdata, handles)
    % Hints: get(hObject,'String') returns contents of text_MRK as text
    %        str2double(get(hObject,'String')) returns contents of text_MRK as
    %        a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function text_MRK_CreateFcn(hObject, eventdata, handles)
    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%--------------------------------------------------------------------------
function text_iMagicCortexMeshLow_Callback(hObject, eventdata, handles)
    % Hints: get(hObject,'String') returns contents of text_iMagicCortexMeshLow as text
    %        str2double(get(hObject,'String')) returns contents of
    %        text_iMagicCortexMeshLow as a double
%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function text_iMagicCortexMeshLow_CreateFcn(hObject, eventdata, handles)
    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%--------------------------------------------------------------------------
function text_iMagicSkinMesh_Callback(hObject, eventdata, handles)
    % Hints: get(hObject,'String') returns contents of text_iMagicSkinMesh as text
    %        str2double(get(hObject,'String')) returns contents of text_iMagicSkinMesh as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function text_iMagicSkinMesh_CreateFcn(hObject, eventdata, handles)
    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end



function text_ELP_Helmet_Callback(hObject, eventdata, handles)
% hObject    handle to text_ELP_Helmet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_ELP_Helmet as text
%        str2double(get(hObject,'String')) returns contents of text_ELP_Helmet as a double


% --- Executes during object creation, after setting all properties.
function text_ELP_Helmet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_ELP_Helmet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function text_ELP_Subject_Callback(hObject, eventdata, handles)
% hObject    handle to text_ELP_Subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_ELP_Subject as text
%        str2double(get(hObject,'String')) returns contents of text_ELP_Subject as a double


% --- Executes during object creation, after setting all properties.
function text_ELP_Subject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_ELP_Subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function text_SEG_Cortex_Callback(hObject, eventdata, handles)
% hObject    handle to text_SEG_Cortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_SEG_Cortex as text
%        str2double(get(hObject,'String')) returns contents of text_SEG_Cortex as a double


% --- Executes during object creation, after setting all properties.
function text_SEG_Cortex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_SEG_Cortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%--------------------------------------------------------------------------
function handles = btn_Browse_ELP_Helmet_Callback(hObject, eventdata, handles)
    
    %Changement de repertoire temporaire (le temps que le fichier soit
    %choisi dans la fenetre de selection de fichiers)
    if handles.load & ~get(handles.btn_Browse_ELP_Helmet,'value')
       Fullname = get(handles.text_ELP_Helmet,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName = fliplr(token);
       PathName = fliplr(remain);
    else
        CurDir = pwd;
        cd(handles.LastDirectory);
        [FileName,PathName] = uigetfile({'*.elp','Standard Polaris File (*.elp)';'*.*','All Files (*.*)'},'Open...'); 
        if FileName == 0 %si pas de fichier on met rien
             set( handles.text_ELP_Helmet, 'String',' ' );
             if isfield(handles,'oDig_CompleteHelmet')
               handles = rmfield(handles, 'oDig_CompleteHelmet');
             end
               guidata(handles.figure_NewProject,handles);
               return
        end
        handles.LastDirectory = PathName;
        cd(CurDir);
    end
    if( isempty(FileName) )
        return;
    elseif( numel(FileName) < 5 )
        disp( 'Invalid File' );
        return;
    end
    
    if(    strcmp( FileName(numel(FileName)-3:numel(FileName)), '.elp') ...
        || strcmp( FileName(numel(FileName)-3:numel(FileName)), '.esp') )
        handles.oDig_CompleteHelmet = Load_Digitization(fullfile(PathName,FileName) );
        set( handles.text_ELP_Helmet, 'String',fullfile(PathName,FileName) );
        guidata(handles.figure_NewProject,handles);
    end
    
    
%--------------------------------------------------------------------------
function handles = btn_Browse_ELP_Subject_Callback(hObject, eventdata, handles)
    
    %Changement de repertoire temporaire (le temps que le fichier soit
    %choisi dans la fenetre de selection de fichiers)
    if handles.load & ~get(handles.btn_Browse_ELP_Subject,'value')
       Fullname = get(handles.text_ELP_Subject,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName = fliplr(token);
       PathName = fliplr(remain);
    else
    CurDir = pwd;
    cd(handles.LastDirectory);
    [FileName,PathName] = uigetfile({'*.elp','Standard Polaris File (*.elp)';'*.*','All Files (*.*)'},'Open...');
    if FileName == 0 %si pas de fichier on met rien
             set( handles.text_ELP_Subject, 'String',' ' );
             if isfield(handles,'oDig_SubjectFiducials')
               handles = rmfield(handles, 'oDig_SubjectFiducials');
             end
             guidata(handles.figure_NewProject,handles);
             return
    end
    
    handles.LastDirectory = PathName;
    cd(CurDir);
    end
    
    if( isempty(FileName) )
        return;
    elseif( numel(FileName) < 5 )
        disp( 'Invalid File' );
        return;
    end
    
    if(    strcmp( FileName(numel(FileName)-3:numel(FileName)), '.elp') ...
        || strcmp( FileName(numel(FileName)-3:numel(FileName)), '.esp') )
        handles.oDig_SubjectFiducials = Load_Digitization( fullfile(PathName,FileName) );
        set( handles.text_ELP_Subject, 'String', fullfile(PathName,FileName) );
        guidata(handles.figure_NewProject,handles);
    end
    
%--------------------------------------------------------------------------
function handles = btn_BrowseVOX_Callback(hObject, eventdata, handles)

    %Init_VOXFileInfo( handles ); 
     if handles.load & ~get(handles.btn_BrowseVOX,'value')
       Fullname = get(handles.text_iMagicVOXfile,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName_Vox = fliplr(token);
       PathName = fliplr(remain);
    else
    CurDir = pwd;
    cd(handles.LastDirectory);
    [FileName_Vox,PathName] = uigetfile({'*.vox','iMagic Volume Data File (*.vox)'},'Open...');
    if FileName_Vox == 0 %si pas de fichier on met rien
             set( handles.text_iMagicVOXfile, 'String',' ' );
             if isfield(handles,'matVoxels')
               handles = rmfield(handles, 'matVoxels');
             end
             if isfield(handles,'bpp')
               handles = rmfield(handles, 'bpp');
             end
             if isfield(handles,'XVoxDim')
               handles = rmfield(handles, 'XVoxDim');
             end
             if isfield(handles,'YVoxDim')
               handles = rmfield(handles, 'YVoxDim');
             end
             if isfield(handles,'ZVoxDim')
               handles = rmfield(handles, 'ZVoxDim');
             end
             if isfield(handles,'VOXFids')
               handles = rmfield(handles, 'VOXFids');
             end
             guidata(handles.figure_NewProject,handles);
             return
    end
    handles.LastDirectory = PathName;
    cd(CurDir);
    end
    if( isempty(FileName_Vox) )
        return;
    end
    
    iExtension = strfind(lower(FileName_Vox), '.vox');
    if( numel(iExtension) == 1 )
        FileName_Hdr = [ FileName_Vox(1:iExtension-1), '.hdr' ];
    else
        warndlg( 'Error loading Volume Data File: Not a .vox file (ERR111)' );
        matVoxels = [];
        return;
    end
    
    FullFileName_Vox = fullfile(PathName,FileName_Vox);
    FullFileName_Hdr = fullfile(PathName,FileName_Hdr);
        
    if( ~exist( FullFileName_Vox, 'file') || ~exist( FullFileName_Hdr, 'file') )
        warndlg( 'Error loading Volume Data File: .hdr or .vox file not found (ERR112)' );
        return;
    end       
    
    %----------------------------------------------------------------------
    %--------------- SECTION HDR (ASCII) ----------------------------------
    %----------------------------------------------------------------------
    fid_Hdr = fopen(FullFileName_Hdr);
    if( fid_Hdr < 2 )
        warndlg( 'Error loading Volume Data File (ERR113)' );
    end

    %XRES
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR114)');
        return;
    end
    XRes=strread(a_line, '%d');
    %set( handles.text_XRes, 'String', sprintf( '%d', XRes ) );
    handles.XRes = XRes;
                
    %YRES
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR115)');
        return;
    end
    YRes=strread(a_line, '%d');
    %set( handles.text_YRes, 'String', sprintf( '%d', YRes ) );
    handles.YRes = YRes;
                
    %ZRES
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR116)');
        return;
    end
    ZRes=strread(a_line, '%d');
    %set( handles.text_ZRes, 'String', sprintf( '%d', ZRes ) );
    handles.ZRes = ZRes;
                
    %Bytes per pixels (jumped)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR117)');
        return;
    end
    BytesPerPix=strread(a_line, '%d');
    handles.bpp = BytesPerPix;
                
    %XZDim (vox dim inside a slice), in mm
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR118)');
        return;
    end
    XZVoxDim = strread(a_line);
    if( numel(XZVoxDim) == 2 )
        handles.XVoxDim = XZVoxDim(1);
        handles.ZVoxDim = XZVoxDim(2);
    else
        handles.XVoxDim = XZVoxDim; %mm
        handles.ZVoxDim = XZVoxDim; %mm
    end
    %set( handles.text_XDim, 'String', sprintf( '%1.4f', handles.XVoxDim ) );
    %set( handles.text_ZDim, 'String', sprintf( '%1.4f', handles.ZVoxDim ) );
    
    %YDim (slice thickness) in mm
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR119)');
        return;
    end
    YDim=strread(a_line);
    handles.YVoxDim = YDim; %mm
    %set( handles.text_YDim, 'String', sprintf( '%1.4f', YDim ) );

    %Plane
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR120)');
        return;
    end
    Plane=strread(a_line, '%d');
    %     if( Plane == 0 )
    %         set( handles.text_Plane, 'String', 'Axial' );
    %     elseif( Plane == 1 )
    %         set( handles.text_Plane, 'String', 'Sagittal' );
    %     elseif( Plane == 2 )
    %         set( handles.text_Plane, 'String', 'Coronal' );
    %     end
    handles.Plane = Plane;

    %Jump this line - Skin treshold
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR121)');
        return;
    end

    %Jump this line - Cortex treshold
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR122)');
        return;
    end
    
    %LPA
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR123)');
        return;
    end
    Fids(2,:)=strread(a_line);
    %if( numel( Fids(1,:) ) == 3 )
    %    set( handles.text_LPA, 'String', sprintf('[%4d %4d %4d]', Fids(2,1), Fids(2,2), Fids(2,3) ) );
    %end
                
    %RPA
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR124)');
        return;
    end
    Fids(3,:)=strread(a_line);
    %if( numel( Fids(2,:) ) == 3 )
    %    set( handles.text_RPA, 'String', sprintf('[%4d %4d %4d]', Fids(3,1), Fids(3,2), Fids(3,3) ) );
    %end
	
    %NAS
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR125)');
        return;
    end
    Fids(1,:)=strread(a_line);
    %if( numel( Fids(3,:) ) == 3 )
    %    set( handles.text_NAS, 'String', sprintf('[%4d %4d %4d]', Fids(1,1), Fids(1,2), Fids(1,3) ) );
    %end                
	
    %Jump this line (VERTEX)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR126)');
        return;
    end
	
    %Jump this line (ANTCOMM? - TAILARACH SYSTEM ORIGIN)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR127)');
        return;
    end
	
    %Jump this line - ATLAS filename associated
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR128)');
        return;
    end
	
    %Jump this line - MIP filename associated
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR129)');
        return;
    end
	
    %Jump this line (INI)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR130)');
        return;
    end
                
    %NOS
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg('Error loading Volume Data File (ERR131)');
        return;
    end
    Fids(4,:)=strread(a_line);
    %if( numel( Fids(4,:) ) == 3 )
    %    set( handles.text_NOS, 'String', sprintf('[%4d %4d %4d]', Fids(4,1), Fids(4,2), Fids(4,3) ) );
    %end
    if isempty(find(Fids))     
        warndlg('No fiducial point marked (ERR132)');
        return;
    end
    handles.VOXFids = Fids;
    
    %----------------------------------------------------------------------
    %--------------- SECTION VOX (BINARY) ---------------------------------
    %----------------------------------------------------------------------
    fid_Vox=fopen(FullFileName_Vox);
    
    if( BytesPerPix == 1 )
        sDataType = 'uint8';
    elseif( BytesPerPix == 2 )
        sDataType = 'uint16';
    else
        warndlg('Error loading Volume Data File: unsupported bbp (ERR132)');
        return;
    end

    %Initialisation
    countTotal = 0;
    MaxPacketSize = 2500000;
    VoxelsToRead = XRes*YRes*ZRes;
    Data1D = zeros(VoxelsToRead, 1, sDataType);
        
    hWaitbar = waitbar(0,'Loading data...');
    while( countTotal < VoxelsToRead )
        
        %Lit en format 'double': convertir en uint8/16
        [Packet, count] =  fread(fid_Vox, MaxPacketSize, sDataType);
        
        if( BytesPerPix == 1 )% 8Bits
            Data1D(countTotal+1:countTotal+count) = uint8(Packet);
        else% 16Bits
            Data1D(countTotal+1:countTotal+count) = uint16(Packet);
        end
        
        countTotal = countTotal+count;
        waitbar(countTotal/VoxelsToRead);
    end
    close(hWaitbar);
    
    if( BytesPerPix == 1 )
        handles.matVoxels = uint8(reshape( Data1D, XRes, ZRes, YRes ));
    else
        handles.matVoxels = uint16(reshape( Data1D, XRes, ZRes, YRes ));
    end    
    set( handles.text_iMagicVOXfile, 'String', FullFileName_Vox );
    guidata(handles.figure_NewProject,handles);

    fclose(fid_Hdr);
    fclose(fid_Vox);
    
%--------------------------------------------------------------------------
function handles = btn_BrowseCortexMeshHi_Callback(hObject, eventdata, handles)
    
    if handles.load & ~get(handles.btn_BrowseCortexMeshHi,'value')
       Fullname = get(handles.text_iMagicCortexMeshHi,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName = fliplr(token);
       PathName = fliplr(remain);
    else
        CurDir = pwd;
        cd(handles.LastDirectory);
        if get(handles.popupmenu_segmentationtype,'value')==1
            [FileName,PathName] = uigetfile({'*.srx','iMagic Surface File (*.srx)';},'Open...');
        elseif get(handles.popupmenu_segmentationtype,'value')==2
            [FileName,PathName] = uigetfile({'*.brain.dfs','BrainSuite Surface File (*.brain.dfs)';},'Open...');
        end
        if FileName == 0 %si pas de fichier on met rien
             set( handles.text_iMagicCortexMeshHi, 'String',' ' );
             if isfield( handles, 'CortexMeshHiVertex' )
               handles = rmfield(handles, 'CortexMeshHiVertex');
             end
             if isfield( handles, 'CortexMeshHiIndex'  )
               handles = rmfield(handles, 'CortexMeshHiIndex' );
             end
             guidata(handles.figure_NewProject,handles);
             return
        end
        handles.LastDirectory = PathName;
        cd(CurDir);
    end
    
    if( isempty(FileName) )
        return;
    end

    file = fullfile(PathName,FileName);
    
    set( handles.text_iMagicCortexMeshHi, 'String', file );   
    
 
    if get(handles.popupmenu_segmentationtype,'value')==1 %Imagic    
        [vertex_matrix, faces_matrix, NV, NT] = Read_iMagic_SRX( file );
        handles.CortexMeshHiVertex = vertex_matrix;
        handles.CortexMeshHiIndex = faces_matrix;
        handles.CortexMeshHiVcolor = zeros(size(handles.CortexMeshHiVertex,1),1);
    elseif get(handles.popupmenu_segmentationtype,'value')==2 %Brainsuite           
        SRX = readdfs(file)
        if isfield(handles,'matVoxels')
            handles.CortexMeshHiVertex = [size(handles.matVoxels,1)-SRX.vertices(:,1),SRX.vertices(:,3),size(handles.matVoxels,2)-SRX.vertices(:,2)] ;
            handles.CortexMeshHiIndex = SRX.faces;
            handles.CortexMeshHiVcolor = zeros(size(handles.CortexMeshHiVertex,1),1);
            nb_x = size(handles.matVoxels,1);
            nb_y = size(handles.matVoxels,2);
            nb_z = size(handles.matVoxels,3);    
            srx2seg(file,nb_x,nb_y,nb_z); 
        else
             set( handles.text_iMagicCortexMeshHi, 'String', '' );   
             msgbox('Please open .VOX file first');
        end
    end
    guidata(handles.figure_NewProject,handles);
    
%--------------------------------------------------------------------------
function handles = btn_BrowseCortexMeshLow_Callback(hObject, eventdata, handles)  
  
     if handles.load & ~get(handles.btn_BrowseCortexMeshLow,'value')
       Fullname = get(handles.text_iMagicCortexMeshLow,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName = fliplr(token);
       PathName = fliplr(remain);
    else
    	CurDir = pwd;
        cd(handles.LastDirectory);
        if get(handles.popupmenu_segmentationtype,'value')==1
            [FileName,PathName] = uigetfile({'*.srx','iMagic Surface File (*.srx)';},'Open...');
        elseif get(handles.popupmenu_segmentationtype,'value')==2
            [FileName,PathName] = uigetfile({'*.brain.dfs','BrainSuite Surface File (*.brain.dfs)';},'Open...');
        end
        if FileName == 0 %si pas de fichier on met rien
             set( handles.text_iMagicCortexMeshLow, 'String',' ' );
             if isfield( handles, 'CortexMeshLowVertex' )
               handles = rmfield(handles, 'CortexMeshLowVertex');
             end
             if isfield( handles, 'CortexMeshLowIndex'  )
               handles = rmfield(handles, 'CortexMeshLowIndex' );
             end
             guidata(handles.figure_NewProject,handles);
             return
        end
        handles.LastDirectory = PathName;
        cd(CurDir);   
    end
    if( isempty(FileName) )
        return;
    end

    file = fullfile(PathName,FileName);
    
    set( handles.text_iMagicCortexMeshLow, 'String', file );    
    
    if get(handles.popupmenu_segmentationtype,'value')==1 %Imagic  
       [vertex_matrix, faces_matrix, NV, NT] = Read_iMagic_SRX( file ); 
        handles.CortexMeshLowVertex = vertex_matrix;
        handles.CortexMeshLowIndex = faces_matrix;
        handles.CortexMeshLowVcolor = zeros(size(handles.CortexMeshLowVertex,1),1);        
    elseif get(handles.popupmenu_segmentationtype,'value')==2 %Brainsuite   
          SRX = readdfs(file)
        if isfield(handles,'matVoxels')
            handles.CortexMeshLowVertex = [size(handles.matVoxels,1)-SRX.vertices(:,1),SRX.vertices(:,3),size(handles.matVoxels,2)-SRX.vertices(:,2)] ;
            handles.CortexMeshLowIndex = SRX.faces;
            handles.CortexMeshLowVcolor = zeros(size(handles.CortexMeshLowVertex,1),1);
            nb_x = size(handles.matVoxels,1);
            nb_y = size(handles.matVoxels,2);
            nb_z = size(handles.matVoxels,3);    
            srx2seg(file,nb_x,nb_y,nb_z); 
            
            handles.load = 1;
            [pathstr, name, ext]=fileparts(file);
            fileseg = fullfile(pathstr,[name,'.seg']);
            set(handles.text_SEG_Cortex,'string',fileseg)
            handles = btn_Browse_SEG_Cortex_Callback(hObject, eventdata, handles);
            handles.load = 0;
        else
            set( handles.text_iMagicCortexMeshHi, 'String', '' );   
            msgbox('Please open .VOX file first');
        end
    end
    guidata(handles.figure_NewProject,handles); 

%--------------------------------------------------------------------------
function handles = btn_BrowseSkinMesh_Callback(hObject, eventdata, handles)
    
    if handles.load & ~get(handles.btn_BrowseSkinMesh,'value')
       Fullname = get(handles.text_iMagicSkinMesh,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName = fliplr(token);
       PathName = fliplr(remain);
    else
        CurDir = pwd;
        cd(handles.LastDirectory);
        if get(handles.popupmenu_segmentationtype,'value')==1
            [FileName,PathName] = uigetfile({'*.srx','iMagic Surface File (*.srx)';},'Open...');
        elseif get(handles.popupmenu_segmentationtype,'value')==2
            [FileName,PathName] = uigetfile({'*.scalp.dfs','BrainSuite Surface File (*.scalp.dfs)';},'Open...');
        end
        if FileName == 0 %si pas de fichier on met rien
             set( handles.text_iMagicSkinMesh, 'String',' ' );
             if isfield( handles, 'SkinMeshVertex' )
               handles = rmfield(handles, 'SkinMeshVertex');
             end
             if isfield( handles, 'SkinMeshIndex'  )
               handles = rmfield(handles, 'SkinMeshIndex' );
             end
             guidata(handles.figure_NewProject,handles);
             return
        end
        handles.LastDirectory = PathName;
        cd(CurDir);
    end
    if( isempty(FileName) )
        return;
    end

    file = fullfile(PathName,FileName);
    
    set( handles.text_iMagicSkinMesh, 'String', file );
    
    if get(handles.popupmenu_segmentationtype,'value')==1 %Imagic 
        [vertex_matrix, faces_matrix, NV, NT] = Read_iMagic_SRX( file );
        handles.SkinMeshVertex = vertex_matrix;
        handles.SkinMeshIndex = faces_matrix;
        handles.SkinVcolor = zeros(size(handles.SkinMeshVertex,1),1);
        
    elseif get(handles.popupmenu_segmentationtype,'value')==2 %Braisuite
        SRX = readdfs(file)
        if isfield(handles,'matVoxels')
            handles.SkinMeshVertex = [size(handles.matVoxels,1)-SRX.vertices(:,1),SRX.vertices(:,3),size(handles.matVoxels,3)-SRX.vertices(:,2)] ;
            handles.SkinMeshIndex = SRX.faces;
            handles.SkinVcolor = zeros(size(handles.SkinMeshVertex,1),1);     
            nb_x = size(handles.matVoxels,1);
            nb_y = size(handles.matVoxels,2);
            nb_z = size(handles.matVoxels,3);    
            srx2seg(file,nb_x,nb_y,nb_z);                  
            handles.load = 1;
            [pathstr, name, ext]=fileparts(file);
            fileseg = fullfile(pathstr,[name,'.seg']);
            
            set(handles.text_SEG_Skin,'string', fileseg)
            handles = btn_Browse_SEG_Skin_Callback(hObject, eventdata, handles);
            handles.load = 0;
        else
            set( handles.text_iMagicCortexMeshHi, 'String', '' );   
            msgbox('Please open .VOX file first');
        end   
     end
    guidata(handles.figure_NewProject,handles);
    
%--------------------------------------------------------------------------
function handles = btn_Browse_SEG_Cortex_Callback(hObject, eventdata, handles)
      
     if handles.load & ~get(handles.btn_Browse_SEG_Cortex,'value')
       Fullname = get(handles.text_SEG_Cortex,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName_Seg = fliplr(token);
       PathName = fliplr(remain);
     else
        CurDir = pwd;
        cd(handles.LastDirectory);
        [FileName_Seg,PathName] = uigetfile({'*.seg','iMagic Segmentation File (*.seg)'}, 'Open...');
         if FileName_Seg == 0 %si pas de fichier on met rien
             set( handles.text_SEG_Cortex, 'String',' ' );
             if isfield( handles, 'matCortexSegmentation')
               handles = rmfield(handles, 'matCortexSegmentation');
             end
             guidata(handles.figure_NewProject,handles);
             return
         end        
        handles.LastDirectory = PathName;
        cd(CurDir);
     end
    if( isempty(FileName_Seg) )
        return;
    end
    
    iExtension = strfind(lower(FileName_Seg), '.seg');
    if( numel(iExtension) == 1 )
        FileName_Hdr = [ FileName_Seg(1:iExtension-1), '.hdr' ];
    else
        warndlg( 'Error loading Segmentation File: Not a .seg file (ERR101)' );
        matSegmentation = [];
        return;
    end
    
    FullFileName_Seg = fullfile(PathName,FileName_Seg);
    FullFileName_Hdr = fullfile(PathName,FileName_Hdr);
        
    if( ~exist( FullFileName_Seg, 'file') || ~exist( FullFileName_Hdr, 'file') )
        warndlg( 'Error loading Segmentation File: .hdr or .seg file not found (ERR102)' );
        return;
    end
    
    set( handles.text_SEG_Cortex, 'String', FullFileName_Seg );
    
    
    
    %----------------------------------------------------------------------
    %--------------- SECTION HDR (ASCII) ----------------------------------
    %----------------------------------------------------------------------
    %    
    %                    ____________________________________
    %                   /.                                  /|
    %                  / .                                 / |
    %                 /  .                                /  |
    %                /   .                               /   |
    %               /    .                              /    |
    %              /     .                             /     |
    %             /      .                            /      |
    %            /__________________________________ /       |
    %           |        .                          |        |
    %           |        .                          |        |
    %           |        .                          |        |
    %           |        .                          |        |
    %           |origin->0 . . . . . . . . . .. . . | . . . .|
    %           |       .                           |       /
    %           |      .                            |      /
    %           |     .                             |     /
    %           |    .                              |    /
    %           |   .                               |   /
    %           |  .                                |  /
    %           | .                                 | /
    %           |___________________________________|/  ->X+        
    %      
    %           Y+ (Slice increment)
    %           |    
    %           |  
    %           | 
    % origin->  0------ X+ (Horizontal increment)
    %          /
    %         /
    %        Z+ (Vertical increment)
    %
    %
    fid_Hdr = fopen(FullFileName_Hdr);
    if( fid_Hdr < 2 )
        warndlg( 'Error loading Segmentation File (ERR103)' );
    end
        
    %Taille X de la segmentation (Horizontal)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR104)' );
        matSegmentation = [];
        return;
    end
    XSegRes=strread(a_line, '%d' );
    handles.XCortexSegRes = XSegRes;
    
    
    %Taille Y de la segmentation (Slice)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR105)' );
        matSegmentation = [];
        return;
    end
    YSegRes=strread(a_line, '%d' );
    handles.YCortexSegRes = YSegRes;
    
                
    %Taille Z de la segmentation (Vertical)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR106)' );
        matSegmentation = [];
        return;
    end
    ZSegRes=strread(a_line, '%d' );
    handles.ZCortexSegRes = ZSegRes;
    
    
    %Index des images
    vIndex = [];
    a_line=fgets(fid_Hdr);
    while( a_line ~= -1 )
        %Lecture de l'index et correction de 0-indexed a 1-indexed
        vIndex(numel(vIndex)+1) = strread(a_line, '%d' )+1;
        a_line=fgets(fid_Hdr);
    end
    MaxSliceNo = max( vIndex );
    handles.vIndexCortexSeg = vIndex;

    %----------------------------------------------------------------------
    %--------------- SECTION SEG (BINARY) ---------------------------------
    %----------------------------------------------------------------------
    fid_Seg=fopen(FullFileName_Seg);
    
    %Initialisation
    countTotal = 0;
    MaxPacketSize = 2500000;
    VoxelsToRead = XSegRes*YSegRes*ZSegRes;
    Data1D = zeros(VoxelsToRead, 1, 'uint8');
        
    hWaitbar = waitbar(0,'Loading data...');
    while( countTotal < VoxelsToRead )
        
        %Lit en format 'double': convertir en uint8/16
        [Packet, count] =  fread(fid_Seg, MaxPacketSize, 'uint8');
        
        Data1D(countTotal+1:countTotal+count) = uint8(Packet);
      
        countTotal = countTotal+count;
        
        waitbar(countTotal/VoxelsToRead);
    end
    close(hWaitbar);
    
    mat_tmp = uint8(reshape( Data1D, XSegRes, ZSegRes, YSegRes ));
    
    handles.matCortexSegmentation = zeros( XSegRes, ZSegRes, MaxSliceNo, 'uint8' );
    
    handles.matCortexSegmentation( :, :, vIndex ) = mat_tmp;
    
    guidata(handles.figure_NewProject,handles);

    fclose(fid_Hdr);
    fclose(fid_Seg);




%-------------------------------------------------------------------------
function handles = btn_Browse_SEG_Skin_Callback(hObject, eventdata, handles)
    if handles.load & ~get(handles.btn_Browse_SEG_Skin,'value')
       Fullname = get(handles.text_SEG_Skin,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName_Seg = fliplr(token);
       PathName = fliplr(remain);
    else
        CurDir = pwd;
        cd(handles.LastDirectory);
        [FileName_Seg,PathName] = uigetfile({'*.seg','iMagic Segmentation File (*.seg)'}, 'Open...');
         if FileName_Seg == 0 %si pas de fichier on met rien
             set( handles.text_SEG_Skin, 'String',' ' );
             if isfield( handles, 'matSkinSegmentation')
               handles = rmfield(handles, 'matSkinSegmentation');
             end
             guidata(handles.figure_NewProject,handles);
             return
         end  
        
        handles.LastDirectory = PathName;
        cd(CurDir);
    end
    
    if( isempty(FileName_Seg) )
        return;
    end
    
    iExtension = strfind(lower(FileName_Seg), '.seg');
    if( numel(iExtension) == 1 )
        FileName_Hdr = [ FileName_Seg(1:iExtension-1), '.hdr' ];
    else
        warndlg( 'Error loading Segmentation File: Not a .seg file (ERR101E)' );
        matSegmentation = [];
        return;
    end
    
    FullFileName_Seg = fullfile(PathName,FileName_Seg);
    FullFileName_Hdr = fullfile(PathName,FileName_Hdr);
        
    if( ~exist( FullFileName_Seg, 'file') || ~exist( FullFileName_Hdr, 'file') )
        warndlg( 'Error loading Segmentation File: .hdr or .seg file not found (ERR102)' );
        return;
    end
    
    set( handles.text_SEG_Skin, 'String', FullFileName_Seg );
    
    
    %----------------------------------------------------------------------
    %--------------- SECTION HDR (ASCII) ----------------------------------
    %----------------------------------------------------------------------
    %    
    %                    ____________________________________
    %                   /.                                  /|
    %                  / .                                 / |
    %                 /  .                                /  |
    %                /   .                               /   |
    %               /    .                              /    |
    %              /     .                             /     |
    %             /      .                            /      |
    %            /__________________________________ /       |
    %           |        .                          |        |
    %           |        .                          |        |
    %           |        .                          |        |
    %           |        .                          |        |
    %           |origin->0 . . . . . . . . . .. . . | . . . .|
    %           |       .                           |       /
    %           |      .                            |      /
    %           |     .                             |     /
    %           |    .                              |    /
    %           |   .                               |   /
    %           |  .                                |  /
    %           | .                                 | /
    %           |___________________________________|/  ->X+        
    %      
    %           Y+ (Slice increment)
    %           |    
    %           |  
    %           | 
    % origin->  0------ X+ (Horizontal increment)
    %          /
    %         /
    %        Z+ (Vertical increment)
    %
    %
    fid_Hdr = fopen(FullFileName_Hdr);
    if( fid_Hdr < 2 )
        warndlg( 'Error loading Segmentation File (ERR103)' );
    end
	
    %Taille X de la segmentation (Horizontal)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR104)' );
        matSegmentation = [];
        return;
    end
    XSegRes=strread(a_line, '%d' );
    %set( handles.text_XSegRes, 'String', sprintf( '%d', XSegRes ) );
	handles.XSkinSegRes = XSegRes;
    
    
    %Taille Y de la segmentation (Slice)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR105)' );
        matSegmentation = [];
        return;
    end
    YSegRes=strread(a_line, '%d' );
    %set( handles.text_YSegRes, 'String', sprintf( '%d', YSegRes ) );
	handles.YSegRes = YSegRes;
    
                
    %Taille Z de la segmentation (Vertical)
    a_line=fgets(fid_Hdr);
    if( a_line == -1 )
        warndlg( 'Error loading Segmentation File (ERR106)' );
        matSegmentation = [];
        return;
    end
    ZSegRes=strread(a_line, '%d' );
    %set( handles.text_ZSegRes, 'String', sprintf( '%d', ZSegRes ) );
	handles.ZSegRes = ZSegRes;
    
    
    %Index des images
    vIndex = [];
    a_line=fgets(fid_Hdr);
    while( a_line ~= -1 )
        %Lecture de l'index et correction de 0-indexed a 1-indexed
        vIndex(numel(vIndex)+1) = strread(a_line, '%d' )+1;
        a_line=fgets(fid_Hdr);
    end
    MaxSliceNo = max( vIndex );
    handles.vIndexSkinSeg = vIndex;

    %----------------------------------------------------------------------
    %--------------- SECTION SEG (BINARY) ---------------------------------
    %----------------------------------------------------------------------
    fid_Seg=fopen(FullFileName_Seg);
    
    %Initialisation
    countTotal = 0;
    MaxPacketSize = 2500000;
    VoxelsToRead = XSegRes*YSegRes*ZSegRes;
    Data1D = zeros(VoxelsToRead, 1, 'uint8');
        
    hWaitbar = waitbar(0,'Loading data...');
    while( countTotal < VoxelsToRead )
        
        %Lit en format 'double': convertir en uint8/16
        [Packet, count] =  fread(fid_Seg, MaxPacketSize, 'uint8');
        
        Data1D(countTotal+1:countTotal+count) = uint8(Packet);
      
        countTotal = countTotal+count;
        
        waitbar(countTotal/VoxelsToRead);
    end
    close(hWaitbar);
    
    mat_tmp = uint8(reshape( Data1D, XSegRes, ZSegRes, YSegRes ));
    
    handles.matSkinSegmentation = zeros( XSegRes, ZSegRes, MaxSliceNo, 'uint8' );
    
    handles.matSkinSegmentation( :, :, vIndex ) = mat_tmp;
    
    guidata(handles.figure_NewProject,handles);

    fclose(fid_Hdr);
    fclose(fid_Seg);
    
%--------------------------------------------------------------------------
function handles = btn_BrowseMRK_Callback(hObject, eventdata, handles)
      
     if handles.load & ~get(handles.btn_BrowseMRK,'value')
       Fullname = get(handles.text_MRK,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName = fliplr(token);
       PathName = fliplr(remain);
     else
        CurDir = pwd;
        cd(handles.LastDirectory);
        [FileName,PathName] = uigetfile({'*.mrk','iMagic Marker File (*.mrk)'}, 'Open...');
        
        handles.LastDirectory = PathName;
        cd(CurDir);
     end
    
    if( ~FileName )
        handles.vMarkers = [];
        set( handles.text_MRK, 'String', '' );
        return;
    else
        set( handles.text_MRK, 'String', fullfile( PathName,FileName ) );
        handles.vMarkers = Read_iMagic_Mrk( fullfile( PathName,FileName ) );
    end
    
    guidata( handles.figure_NewProject, handles );
    



function text_MTG_Helmet_Callback(hObject, eventdata, handles)
% hObject    handle to text_MTG_Helmet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_MTG_Helmet as text
%        str2double(get(hObject,'String')) returns contents of text_MTG_Helmet as a double


% --- Executes during object creation, after setting all properties.
function text_MTG_Helmet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_MTG_Helmet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_MTg_Callback(hObject, eventdata, handles)
% hObject    handle to text_MTg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_MTg as text
%        str2double(get(hObject,'String')) returns contents of text_MTg as a double


% --- Executes during object creation, after setting all properties.
function text_MTg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_MTg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in btn_Browse_MTG_Subject.
function handles = btn_Browse_MTG_Subject_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Browse_MTG_Subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if handles.load & ~get(handles.btn_Browse_MTG_Subject,'value')
       Fullname = get(handles.text_MTG,'string');
       [token, remain] = strtok( fliplr(Fullname),'\');
       FileName = fliplr(token);
       PathName = fliplr(remain);
    else
        CurDir = pwd;
        cd(handles.LastDirectory);
        [FileName,PathName] = uigetfile({'*.mtg','iomtg montage created';}, 'Open...');
        if FileName == 0 %si pas de fichier on met rien
             set( handles.text_MTG, 'String',' ' );    
         end   
        handles.LastDirectory = PathName;
        cd(CurDir);
    end
    
    if( ~FileName )
        set( handles.text_MTG, 'String', '' );
        return;
    else
        set( handles.text_MTG, 'String', fullfile( PathName,FileName ) );
    end
    handles.oDig_MTGimport = load([PathName,FileName],'-mat');
    
    guidata( handles.figure_NewProject, handles );




function text_MTG_Callback(hObject, eventdata, handles)
% hObject    handle to text_MTG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_MTG as text
%        str2double(get(hObject,'String')) returns contents of text_MTG as a double


% --- Executes during object creation, after setting all properties.
function text_MTG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_MTG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function EDIT_SPMnor_Callback(hObject, eventdata, handles)
% hObject    handle to text_spmnor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_spmnor as text
%        str2double(get(hObject,'String')) returns contents of text_spmnor as a double


% --- Executes during object creation, after setting all properties.
function EDIT_SPMnor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_spmnor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_SPMnor.
function handles = btn_SPMnor_Callback(hObject, eventdata, handles)
% hObject    handle to btn_SPMnor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if handles.load & ~get(handles.btn_SPMnor,'value')
       Fullname = get(handles.text_SPMnor,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName = fliplr(token);
       PathName = fliplr(remain);
    else
        CurDir = pwd;
        cd(handles.LastDirectory);
        [FileName,PathName] = uigetfile({'*.mat','SPM normalization transformation';}, 'Open...');
        handles.LastDirectory = PathName;
        cd(CurDir);
    end
    
    if( ~FileName )
        set( handles.text_SPMnor, 'String', '' );
        return;
    else
        set( handles.text_SPMnor, 'String', fullfile( PathName,FileName ) );
        %handles.sn_file= load( fullfile( PathName,FileName ) );
        handles.sn_file = [PathName,FileName];
    end
          
    guidata( handles.figure_NewProject, handles );



function text_SPMnor_Callback(hObject, eventdata, handles)
% hObject    handle to text_SPMnor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_SPMnor as text
%        str2double(get(hObject,'String')) returns contents of text_SPMnor as a double


% --- Executes during object creation, after setting all properties.
function text_SPMnor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_SPMnor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function btn_load_path_Callback(hObject, eventdata, handles)
% hObject    handle to btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('.path');
if filename == 0
    return
end

fid = fopen([pathname, filename], 'r');
handles.load = 1;

name = fgetl(fid);
if ~isempty(name)
    set(handles.text_ELP_Helmet,'string',name);
    %handles = btn_Browse_ELP_Helmet_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if  ~isempty(name)
    set(handles.text_ELP_Subject,'string',name);
    %handles = btn_Browse_ELP_Subject_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_MTG,'string',name);
    %handles = btn_Browse_MTG_Subject_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_iMagicVOXfile,'string', name);
    %handles = btn_BrowseVOX_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_iMagicCortexMeshHi,'string',name);
    %handles = btn_BrowseCortexMeshHi_Callback(hObject, eventdata, handles);
end 
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_iMagicCortexMeshLow,'string',name);
    %handles = btn_BrowseCortexMeshLow_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_iMagicSkinMesh,'string',name);
    %handles = btn_BrowseSkinMesh_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_SEG_Cortex,'string',name);
    %handles = btn_Browse_SEG_Cortex_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_SEG_Skin,'string',name);
    %handles = btn_Browse_SEG_Skin_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_MRK,'string',name);
    %handles = btn_BrowseMRK_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)
    set(handles.text_SPMnor,'string',name);
    %handles = btn_SPMnor_Callback(hObject, eventdata, handles);
end
name = fgetl(fid);
if ~isempty(name)& name~=-1
    set(handles.text_atlas,'string',name);
    %handles = btn_browse_atlas_Callback(hObject, eventdata, handles)
end
name = fgetl(fid);
if ~isempty(name)& name~=-1
    set(handles.popupmenu_segmentationtype,'value',str2num(name));
else
    set(handles.popupmenu_segmentationtype,'value',1);    
end

 popupmenu_segmentationtype_Callback(hObject, eventdata, handles); 

if isfield(handles,'text_ELP_Helmet')
    handles = btn_Browse_ELP_Helmet_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_ELP_Subject')
    handles = btn_Browse_ELP_Subject_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_MTG')
    handles = btn_Browse_MTG_Subject_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_iMagicVOXfile')
    handles = btn_BrowseVOX_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_iMagicCortexMeshHi')
    handles = btn_BrowseCortexMeshHi_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_iMagicCortexMeshLow')
    handles = btn_BrowseCortexMeshLow_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_iMagicSkinMesh')
    handles = btn_BrowseSkinMesh_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_SEG_Cortex')
    handles = btn_Browse_SEG_Cortex_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_SEG_Skin')
    handles = btn_Browse_SEG_Skin_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_MRK')
     handles = btn_BrowseMRK_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_SPMnor')
    handles = btn_SPMnor_Callback(hObject, eventdata, handles);
end
if isfield(handles,'text_atlas')
    handles = btn_browse_atlas_Callback(hObject, eventdata, handles);
end


handles.load = 0;
fclose(fid);
guidata( handles.figure_NewProject, handles );

% --- Executes on button press in btn_save_path.
function btn_save_path_Callback(hObject, eventdata, handles,default)
% hObject    handle to btn_load_path_save_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if default == 0
    [filename, pathname] = uiputfile('.path');
end
if filename == 0
    return
end

fid = fopen([pathname,filename],'w');
fprintf(fid,'%s\n',get(handles.text_ELP_Helmet,'string'));
fprintf(fid,'%s\n',get(handles.text_ELP_Subject,'string'));
fprintf(fid,'%s\n',get(handles.text_MTG,'string'));
fprintf(fid,'%s\n',get(handles.text_iMagicVOXfile,'string'));
fprintf(fid,'%s\n',get(handles.text_iMagicCortexMeshHi,'string'));
fprintf(fid,'%s\n',get(handles.text_iMagicCortexMeshLow,'string'));
fprintf(fid,'%s\n',get(handles.text_iMagicSkinMesh,'string'));
fprintf(fid,'%s\n',get(handles.text_SEG_Cortex,'string'));
fprintf(fid,'%s\n',get(handles.text_SEG_Skin,'string'));
fprintf(fid,'%s\n',get(handles.text_MRK,'string'));
fprintf(fid,'%s\n',get(handles.text_SPMnor,'string'));
fprintf(fid,'%s\n',get(handles.text_atlas,'string'));
fprintf(fid,'%s\n',num2str(get(handles.popupmenu_segmentationtype,'value')));
fclose(fid);
guidata( handles.figure_NewProject, handles );



function text_atlas_Callback(hObject, eventdata, handles)
% hObject    handle to text_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_atlas as text
%        str2double(get(hObject,'String')) returns contents of text_atlas as a double


% --- Executes during object creation, after setting all properties.
function text_atlas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_browse_atlas.
function handles = btn_browse_atlas_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browse_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.load & ~get(handles.btn_browse_atlas,'value')
       Fullname = get(handles.text_atlas,'string');
       [token, remain] = strtok( fliplr(Fullname),{'\','/'});
       FileName = fliplr(token);
       PathName = fliplr(remain);
    else
        CurDir = pwd;
        cd(handles.LastDirectory);
        [FileName,PathName] = uigetfile({'*.img','atlas';}, 'Open...');
        if FileName == 0
             set( handles.text_atlas, 'String',' ' );
             if isfield( handles, 'ImageAtlas')
               handles = rmfield(handles, 'ImageAtlas'); 
             end
             guidata(handles.figure_NewProject,handles);
             return
        end
        handles.LastDirectory = PathName;
        cd(CurDir);
        
    end
    handles.ImageAtlas =  [PathName,FileName];
    set(handles.text_atlas,'string',handles.ImageAtlas)   
    guidata( handles.figure_NewProject, handles );


% --- Executes on selection change in popupmenu_segmentationtype.
function popupmenu_segmentationtype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_segmentationtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_segmentationtype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_segmentationtype
if get(handles.popupmenu_segmentationtype,'value')==1 % 'imagic'
   set(handles.uipanel7,'title','iMagic Surface Files (.SRX files)');
   set(handles.text91,'string', 'Cortex Surface Mesh for atlas display ');
   set(handles.text93,'string', 'Cortex Surface Mesh (Low-Res)');
   set(handles.text92,'string', 'Skin Surface Mesh (Hi-Res)');
 elseif get(handles.popupmenu_segmentationtype,'value')==2 %BrainSuite
    set(handles.uipanel7,'title','BrainSuite Surface Files (.dfs files)')
    set(handles.text91,'string', 'Cortex Surface Mesh for atlas display *.brain.dfs');
    set(handles.text93,'string', 'Cortex Surface Mesh (Low-Res) *.brain.dfs');
    set(handles.text92,'string', 'Skin Surface Mesh (Hi-Res) *.scalp.dfs');  
end

% --- Executes during object creation, after setting all properties.
function popupmenu_segmentationtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_segmentationtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_transformationfit_Callback(hObject, eventdata, handles)
% hObject    handle to edit_transformationfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_transformationfit as text
%        str2double(get(hObject,'String')) returns contents of edit_transformationfit as a double


% --- Executes during object creation, after setting all properties.
function edit_transformationfit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_transformationfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_browseTransfromationMatFit.
function btn_browseTransfromationMatFit_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browseTransfromationMatFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  

        [FileName,PathName] = uigetfile({'*.mat','Transformation matrice manuel use in the fit'}, 'Open...');
        set( handles.edit_transformationfit  , 'String', fullfile( PathName,FileName ) );
        guidata( handles.figure_NewProject, handles );
    
