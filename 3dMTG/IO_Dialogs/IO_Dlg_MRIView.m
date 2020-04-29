function varargout = IO_Dlg_MRIView(varargin)
% IO_DLG_MRIVIEW M-file for IO_Dlg_MRIView.fig
%      IO_DLG_MRIVIEW, by itself, creates a new IO_DLG_MRIVIEW or raises the existing
%      singleton*.
%
%      H = IO_DLG_MRIVIEW returns the handle to a new IO_DLG_MRIVIEW or the handle to
%      the existing singleton*.
%
%      IO_DLG_MRIVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IO_DLG_MRIVIEW.M with the given input arguments.
%
%      IO_DLG_MRIVIEW('Property','Value',...) creates a new IO_DLG_MRIVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IO_Dlg_MRIView_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IO_Dlg_MRIView_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IO_Dlg_MRIView

% Last Modified by GUIDE v2.5 17-Apr-2007 17:24:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IO_Dlg_MRIView_OpeningFcn, ...
                   'gui_OutputFcn',  @IO_Dlg_MRIView_OutputFcn, ...
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


% --- Executes just before IO_Dlg_MRIView is made visible.
function IO_Dlg_MRIView_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to IO_Dlg_MRIView (see VARARGIN)

    % UIWAIT makes IO_Dlg_MRIView wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

    if( ~numel( varargin ) || ~isa( varargin{1}, 'IO_Project_Data' ) )
        hDlg = errordlg( 'Error in IO_Dlg_MRIView_OpeningFcn: varargin not a ''IO_Project_Data''. Loading Defaults' );
        uiwait(hDlg);
        oPrjData = IO_Project_Data;
    else
        oPrjData = varargin{1};
    end
    
    oMRI = get_MRI_Data( oPrjData );
    
    %Fit_SrcDets_OnSkinSegmentation( handles, get_Helmet(oPrjData), oMRI)
    oHelmet = Fit_Helmet_OnMRISegmentations( get_Helmet(oPrjData), oMRI, true, true, true );
    %vHoles_SkinDepth = Fit_SrcDets_OnSkinSegmentation( get_Helmet(oPrjData), oMRI);
    
    oPrjData = set_Helmet( oPrjData, oHelmet );
    
    cameratoolbar( handles.figure1 );
    
    %Sauvegarde de l'objet
    handles.m_oPrjData = oPrjData;
    
    % Choose default command line output for IO_Dlg_DisplayOptions
    handles.output = oPrjData;
    
    
    % Update handles structure
    guidata(hObject, handles);
    
    
    
    %[Helm_MRI_Coords, Helm_MRI_Norms ] = CoRegisterHelmetToMRI( get_Helmet(oPrjData), oMRI );
    
    %handles.Helm_MRI_Coords = Helm_MRI_Coords(:,1:3);
    %handles.Helm_MRI_Norms = Helm_MRI_Norms(:,1:3);
    
    matSkinSegmentation = get_SkinSegmentation( oMRI );
    [nX,nY,nZ] = size( matSkinSegmentation );
    
    
    
    oMRI = set_SkinSegmentation( oMRI, matSkinSegmentation );
    handles.m_oPrjData = set_MRI_Data( oPrjData, oMRI );
    IMax = get_IMaxVox( oMRI );
    
    handles.AxialSlice = 100;
    
    axes(handles.axe_MRI_Sli);
    img = get_2DImage( oMRI, 2, handles.AxialSlice );
    seg = get_2DSegImage( oMRI, 2, handles.AxialSlice );
    if( isempty(img) || isempty(seg) )
        return;
    end
    segVox = find( seg );
    [ nRows, nCols ] = size( img );
    imgRGB = [ img, img, img ];
    imgRGB( [segVox, segVox+nRows*nCols*2] ) = zeros( 1, numel(segVox)*2, 'uint8' );
    imgRGB( segVox+nRows*nCols ) = 255*ones( 1, numel(segVox), 'uint8' );
    image( reshape(imgRGB,nRows,nCols,3)  ); 
    axis image;
	axis off;
    
    axes(handles.axe_MRI_Row);
    img = get_2DImage( oMRI, 3, 100 );
    seg = get_2DSegImage( oMRI, 3, 100 );
    segVox = find( seg );
    [ nRows, nCols ] = size( img );
    imgRGB = [ img, img, img ];
    imgRGB( [segVox, segVox+nRows*nCols*2] ) = zeros( 1, numel(segVox)*2, 'uint8' );
    imgRGB( segVox+nRows*nCols ) = 255*ones( 1, numel(segVox), 'uint8' );
    image( reshape(imgRGB,nRows,nCols,3)  ); 
    axis image;
	axis off;
    
    axes(handles.axe_MRI_Col);
    img = get_2DImage( oMRI, 1, 100 );
    seg = get_2DSegImage( oMRI, 1, 100 );
    segVox = find( seg );
    [ nRows, nCols ] = size( img );
    imgRGB = [ img, img, img ];
    imgRGB( [segVox, segVox+nRows*nCols*2] ) = zeros( 1, numel(segVox)*2, 'uint8' );
    imgRGB( segVox+nRows*nCols ) = 255*ones( 1, numel(segVox), 'uint8' );
    image( reshape(imgRGB,nRows,nCols,3)  ); 
    axis image;
	axis off;
    
%     set( handles.slider1, 'max', handles.maxIntensity );
%     set( handles.slider1, 'min', 0 );
%     set( handles.slider1, 'Value', handles.maxIntensity );
%     set( handles.txt_IMax, 'String', sprintf( '%d', handles.maxIntensity ) );
    
    
    
    % Update handles structure
    guidata(hObject, handles);
    


% --- Outputs from this function are returned to the command line.
function varargout = IO_Dlg_MRIView_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

    %Waits for OK or Cancel before setting output
    uiwait(hObject);
    
    %If 'hObject' is not a handle anymore, user has exited Dlg with 'X'
    if( ishandle(hObject) )
        handles = guidata( hObject );
        varargout{1} = handles.output;
        close(hObject);
    end
    

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    Int = get(hObject,'Value');
    %disp( Int );
    set(handles.txt_IMax, 'String', sprintf('%d',int16(Int)) );
    axes( handles.axe_MRI_Sli );
    colormap( gray( Int ) );
    axes( handles.axe_MRI_Row );
    colormap( gray( Int ) );
    axes( handles.axe_MRI_Col );
    colormap( gray( Int ) );
    
    
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function txt_XCoord_Callback(hObject, eventdata, handles)
% hObject    handle to txt_XCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_XCoord as text
%        str2double(get(hObject,'String')) returns contents of txt_XCoord as a double


% --- Executes during object creation, after setting all properties.
function txt_XCoord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_XCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_YCoord_Callback(hObject, eventdata, handles)
% hObject    handle to txt_YCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_YCoord as text
%        str2double(get(hObject,'String')) returns contents of txt_YCoord as a double


% --- Executes during object creation, after setting all properties.
function txt_YCoord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_YCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ZCoord_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ZCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ZCoord as text
%        str2double(get(hObject,'String')) returns contents of txt_ZCoord as a double


% --- Executes during object creation, after setting all properties.
function txt_ZCoord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ZCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axe_MRI_Sli_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to axe_MRI_Sli (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: place code in OpeningFcn to populate axe_MRI_Sli
%     disp( 'axe_MRI_Sli_CreateFcn' );
%     SizeX = size( handles.matVoxelsXY, 1 );
%     SizeY = size( handles.matVoxelsXY, 2 );
%     SizeZ = size( handles.matVoxelsXY, 3 );
%     PosX = ceil(SizeX/2);
%     PosY = ceil(SizeY/2);
%     PosZ = ceil(SizeZ/2);
%     if( PosX && PosY && PosZ )
%         axes(hObject);
%         %colormap(gray(1200));
%         image( handles.matVoxelsXY(:,:,PosZ) );
%     end
    
    

% --- Executes during object creation, after setting all properties.
function axe_MRI_Row_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to axe_MRI_Row (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: place code in OpeningFcn to populate axe_MRI_Row
%     disp( 'axe_MRI_Row_CreateFcn' );
%     SizeX = size( handles.matVoxelsXY, 1 );
%     SizeY = size( handles.matVoxelsXY, 2 );
%     SizeZ = size( handles.matVoxelsXY, 3 );
%     PosX = ceil(SizeX/2);
%     PosY = ceil(SizeY/2);
%     PosZ = ceil(SizeZ/2);
%     if( PosX && PosY && PosZ )
%         axes(hObject);
%         %colormap(gray(1200));
%         image( reshape( handles.matVoxelsXY(:,PosY,:), SizeX, SizeZ ) );
%     end

% --- Executes during object creation, after setting all properties.
function axe_MRI_Col_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to axe_MRI_Col (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

%     disp( 'axe_MRI_Col_CreateFcn' );
%     SizeX = size( handles.matVoxelsXY, 1 );
%     SizeY = size( handles.matVoxelsXY, 2 );
%     SizeZ = size( handles.matVoxelsXY, 3 );
%     PosX = ceil(SizeX/2);
%     PosY = ceil(SizeY/2);
%     PosZ = ceil(SizeZ/2);
%     if( PosX && PosY && PosZ )
%         axes(hObject);
%         %colormap(gray(1200));
%         image( reshape( handles.matVoxelsXY(PosX,:,:), SizeY, SizeZ ) );
%     end





function txt_IMax_Callback(hObject, eventdata, handles)
% hObject    handle to txt_IMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_IMax as text
%        str2double(get(hObject,'String')) returns contents of txt_IMax as a double


% --- Executes during object creation, after setting all properties.
function txt_IMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_IMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Axial_Up.
function Axial_Up_Callback(hObject, eventdata, handles)
    % hObject    handle to Axial_Up (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    global PrjData;
    
    IMax = get_IMaxVox( get_MRI_Data( PrjData ) );
    handles.AxialSlice = handles.AxialSlice+1;
    
    axes(handles.axe_MRI_Sli);
    img = get_2DImage( get_MRI_Data( PrjData ), 2, handles.AxialSlice );
    seg = get_2DSegImage( get_MRI_Data( PrjData ), 2, handles.AxialSlice );
    if( isempty(img) || isempty(seg) )
        return;
    end
    segVox = find( seg );
    [ nRows, nCols ] = size( img );
    imgRGB = [ img, img, img ];
    imgRGB( [segVox, segVox+nRows*nCols*2] ) = zeros( 1, numel(segVox)*2, 'uint8' );
    imgRGB( segVox+nRows*nCols ) = 255*ones( 1, numel(segVox), 'uint8' );
    image( reshape(imgRGB,nRows,nCols,3)  );
    axis image;
	axis off;
    guidata( handles.figure1, handles );

% --- Executes on button press in Axial_Down.
function Axial_Down_Callback(hObject, eventdata, handles)
% hObject    handle to Axial_Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global PrjData;
    
    IMax = get_IMaxVox( get_MRI_Data( PrjData ) );
    handles.AxialSlice = handles.AxialSlice-1;
    
    axes(handles.axe_MRI_Sli);
    img = get_2DImage( get_MRI_Data( PrjData ), 2, handles.AxialSlice );
    %seg = get_2DSegImage( get_MRI_Data( PrjData ), 2, handles.AxialSlice );
    seg = get_2DSegImage( get_MRI_Data( PrjData ), 2, handles.AxialSlice );
    if( isempty(img) || isempty(seg) )
        return;
    end
    segVox = find( seg );
    [ nRows, nCols ] = size( img );
    imgRGB = [ img, img, img ];
    imgRGB( [segVox, segVox+nRows*nCols*2] ) = zeros( 1, numel(segVox)*2, 'uint8' );
    imgRGB( segVox+nRows*nCols ) = 255*ones( 1, numel(segVox), 'uint8' );
    image( reshape(imgRGB,nRows,nCols,3)  ); 
    axis image;
	axis off;
    guidata( handles.figure1, handles );

% function [Helm_MRI_Coords, Helm_MRI_Norms ] = Coreg_Helmet_To_MRI(oHelmet, oMRI_Data)
%     
%     vHolesH = get_vHoles( oHelmet );
%     sMtg = get_Mtg( oHelmet );
% 
%     %Transfert des coordonnees sous forme de matrices
%     CoordsH = [vHolesH(find([vHolesH.Type] == 400)).Coord];
%     NormalsH = [vHolesH(find([vHolesH.Type] == 400)).Normal];
%     matCoordsH = [CoordsH.x;CoordsH.y;CoordsH.z]';
%     matNormalsH = [NormalsH.x;NormalsH.y;NormalsH.z]';
% 
%     matHelmetRefs = sMtg.matFiducials;
% 
%     matMRIRefs = get_matFiducials( oMRI_Data );
% 
% 	%Construction de la matrice de coregistration
%     IsHelmetRH = true;
%     IsMRIRH = false;
%     Physical_HelmCoordSystemSize = [1,1,1]; %Le casque provient d'un systeme orthonorme et metrique
%     matCoRegistration = Create_CoRegistration_Matrix( matHelmetRefs,                matMRIRefs, ...
%                                                       Physical_HelmCoordSystemSize, get_PhysicalVoxDim( oMRI_Data ), ...
%                                                       IsHelmetRH,                   IsMRIRH );
%     matCoordsH(:,4) = 1;
%     matNormalsH(:,4) = 1;
% 
%     %Transfert de systeme de coords
%     Helm_MRI_Coords = matCoordsH*matCoRegistration;
%     Helm_MRI_Norms = matNormalsH*matCoRegistration;
% 
%     %Normales unitaires (en voxel-space)
%     for( iNorm=1:size(Helm_MRI_Norms,1) )
%         Helm_MRI_Norms(iNorm,1:3) = Helm_MRI_Norms(iNorm,1:3) / norm(Helm_MRI_Norms(iNorm,1:3));
%     end

% function Fit_SrcDets_OnSkinSegmentation(handles, oHelmet, oMRI_Data)
% 
%     HCoords = handles.Helm_MRI_Coords;
%     HNorms = handles.Helm_MRI_Norms;
%     matSkinSegmentation = get_SkinSegmentation( oMRI_Data );
%     [nVoxX,nVoxY,nVoxZ] = size( matSkinSegmentation );
%     
%     vPhysicalVoxDim = get_PhysicalVoxDim( oMRI_Data );
%     if( numel(find( vPhysicalVoxDim == vPhysicalVoxDim(1) ) ) < 3 )
%         warndlg( 'Voxel size is not cubic: possible errors in distances calculated on skin segmentation (ERR204)');
%     end
%     MaxDepth_m = 0.05; %in meters (5 cm)
%     MaxDepth_vox = ceil(MaxDepth_m/mean(vPhysicalVoxDim));
%     
%     %Structure de memorisation des profondeurs de trous.
%     H_SkinDepth = zeros( 1,size(HCoords,1) );
%     
%     vHoles = get_vHoles(oHelmet);
%     v_iHolesToDisplay = find([vHoles.Type] == 400);
%     
%     hold on;
%     [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI_Data );
%     patch('Vertices',VertexBuffer,'Faces',IndexBuffer, 'EdgeColor', 'none', 'FaceColor', [0.5,0.5,0.5], 'FaceLighting', 'gouraud', 'SpecularColorReflectance', 0.5, 'DiffuseStrength', 0.5, 'SpecularStrength', 1);
%     light( 'Position', [500,500,500] );
%     light( 'Position', [500,-500,500] );
%     light( 'Position', [-500,0,0] );
%     
%     %Pour chaque position de trou ce casque (en voxel-space)
%     for( iElem = 1:size(HCoords,1) )
%         
%         DepthSteps = 0:1:MaxDepth_vox;
%         HCoordsStepsX = int16(HCoords(iElem,1)-HNorms(iElem,1)*DepthSteps);
%         HCoordsStepsY = int16(HCoords(iElem,2)-HNorms(iElem,2)*DepthSteps);
%         HCoordsStepsZ = int16(HCoords(iElem,3)-HNorms(iElem,3)*DepthSteps);
%         
%         SkinReached = false;
%         Depth = 0;
%         
%         %hold on;
%         plot3( HCoordsStepsX, HCoordsStepsY,HCoordsStepsZ, 'r' );
%            
%         %Pour au moins chaque voxel dans la trajectoire de la normale
%         for( iDepth=1:numel(DepthSteps) )
%             
%             %Si a l'interieur de la limite de voxels (-1 et +1 pour cube)
%             if(   ( HCoordsStepsX(iDepth) > 1 && HCoordsStepsX(iDepth) < nVoxX ) ...
%                 &&( HCoordsStepsY(iDepth) > 1 && HCoordsStepsY(iDepth) < nVoxY ) ...
%                 &&( HCoordsStepsZ(iDepth) > 1 && HCoordsStepsZ(iDepth) < nVoxZ ) )
%                 
%                 %Verifier que le nombre de voxels du cube centre en
%                 %HCoordsStepsX/Y/Z faisant partie de la segmentation de la
%                 %peau est d'au moins 13 (13/27 des voxels du cube)
%                 nVoxelsInSeg = numel( find( matSkinSegmentation(HCoordsStepsX(iDepth)-1:HCoordsStepsX(iDepth)+1, ...
%                                                                 HCoordsStepsY(iDepth)-1:HCoordsStepsY(iDepth)+1, ...
%                                                                 HCoordsStepsZ(iDepth)-1:HCoordsStepsZ(iDepth)+1 ) ) );
%                 if( numel( find( matSkinSegmentation(HCoordsStepsX(iDepth)-1:HCoordsStepsX(iDepth)+1, ...
%                                                      HCoordsStepsY(iDepth)-1:HCoordsStepsY(iDepth)+1, ...
%                                                      HCoordsStepsZ(iDepth)-1:HCoordsStepsZ(iDepth)+1 ) ) ) >= 13 )
%                     SkinReached = true;
%                     Depth = (iDepth-1)*HNorms(iElem,1:3)*vPhysicalVoxDim';
%                     plot3( [HCoordsStepsX(iDepth),HCoordsStepsX(iDepth)-1], [HCoordsStepsY(iDepth),HCoordsStepsY(iDepth)-1],[HCoordsStepsZ(iDepth),HCoordsStepsZ(iDepth)-1], 'color', [1,1,1], 'linewidth', 3 );
%                     break;
%                 end
%             end
%         end
%            
%         H_SkinDepth(iElem) = Depth;
%     end
%     
%     cameratoolbar( handles.figure1 );
