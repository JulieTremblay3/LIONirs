function varargout = IO_Dlg_DisplayOptions(varargin)
    % IO_DLG_DISPLAYOPTIONS M-file for IO_Dlg_DisplayOptions.fig
    %      IO_DLG_DISPLAYOPTIONS, by itself, creates a new IO_DLG_DISPLAYOPTIONS or raises the existing
    %      singleton*.
    %
    %      H = IO_DLG_DISPLAYOPTIONS returns the handle to a new IO_DLG_DISPLAYOPTIONS or the handle to
    %      the existing singleton*.
    %
    %      IO_DLG_DISPLAYOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in IO_DLG_DISPLAYOPTIONS.M with the given input arguments.
    %
    %      IO_DLG_DISPLAYOPTIONS('Property','Value',...) creates a new IO_DLG_DISPLAYOPTIONS or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before IO_Dlg_DisplayOptions_OpeningFunction gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to IO_Dlg_DisplayOptions_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help IO_Dlg_DisplayOptions

    % Last Modified by GUIDE v2.5 29-Nov-2017 11:25:25

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @IO_Dlg_DisplayOptions_OpeningFcn, ...
                       'gui_OutputFcn',  @IO_Dlg_DisplayOptions_OutputFcn, ...
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
% --- Executes just before IO_Dlg_DisplayOptions is made visible.
function IO_Dlg_DisplayOptions_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to IO_Dlg_DisplayOptions (see VARARGIN)

    % UIWAIT makes IO_Dlg_DisplayOptions wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
    
    if( ~numel( varargin ) || ~isa( varargin{1}, 'IO_DisplayOptions' ) )
        hDlg = errordlg( 'Error in IO_Dlg_DisplayOptions_OpeningFcn: varargin not a ''IO_DisplayOptions''. Loading Defaults' );
        uiwait(hDlg);
        oDispOpt = IO_DisplayOptions;
    else
        oDispOpt = varargin{1};
    end
    
    % Choose default command line output for IO_Dlg_DisplayOptions
    handles.output = oDispOpt;
    
    %Sauvegarde de l'objet
    handles.m_oDispOpt = oDispOpt;
    
    % Update handles structure
    guidata(hObject, handles);
    
    %vecteur de handles de checkboxes
    v_handlesCheckBoxes = findobj(handles.figure_DisplayOptions,'Style','checkbox');
    prefix = 'chk_';
    
    %Chargement des options d'affichage
    for( iBtn=1:length(v_handlesCheckBoxes) )
        
        %Handle du checkbox
        hChk = v_handlesCheckBoxes(iBtn);
        
        %Tag du checkbox
        BtnTag = get( hChk, 'Tag' );
        
        %String de l'option d'affichage
        strDisplayOption = BtnTag(5:numel(BtnTag));
        
        %booleen : l'option est elle cochee?
        bIsActive = get_DispOptChecked( oDispOpt, strDisplayOption );
        
        %Coche l'option si necessaire
        set( hChk, 'Value', bIsActive );
    end
    
    Diameter = get_SphereDiameter(oDispOpt);
    set(handles.string_diameter,'string',Diameter);    
    atlasnum = get_DispAtlasZone(oDispOpt);
    set(handles.edit_atlasnum,'string',atlasnum);
    Interval = get_HoleIdInterval( oDispOpt );
    popupStrings =  get(handles.popup_HoleIdInterval, 'String' );
    for( iString=1:numel(popupStrings) )
        if( strcmp( popupStrings{iString}, sprintf( '%d', Interval ) ) )
            set( handles.popup_HoleIdInterval, 'Value', iString );
            break;
        end
    end
    
    
    
%--------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = IO_Dlg_DisplayOptions_OutputFcn(hObject, eventdata, handles) 
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
    
    
%--------------------------------------------------------------------------
function btn_Cancel_Callback(hObject, eventdata, handles)
    %Unlocks IO_Dlg_DisplayOptions_OutputFcn uiwait
    uiresume(handles.figure_DisplayOptions);
    
%--------------------------------------------------------------------------
function btn_Ok_Callback(hObject, eventdata, handles)
    
    % Command line output for IO_Dlg_DisplayOptions
    handles.output = SaveChanges( handles.m_oDispOpt, handles );
    
    % Saves output into guidata
    guidata(handles.figure_DisplayOptions, handles);
    
    %Unlocks IO_Dlg_DisplayOptions_OutputFcn uiwait
    uiresume(handles.figure_DisplayOptions);

    
function oDispOpt = SaveChanges( oDispOpt, handles )

    %vecteur de handles de checkboxes
    v_handlesCheckBoxes = findobj( handles.figure_DisplayOptions,'Style','checkbox' );
    prefix = 'chk_';
    
    %Sauvegarde des options d'affichage dans l'object IO_DisplayOptions
    for( iBtn=1:length(v_handlesCheckBoxes) )
        
        %Handle du checkbox
        hChk = v_handlesCheckBoxes(iBtn);
        
        %Tag du checkbox
        BtnTag = get( hChk, 'Tag' );
        
        %String de l'option d'affichage (meme que checkbox sans le 'chk_')
        strDisplayOption = BtnTag(5:numel(BtnTag));
        
        %L'objet prend la valeur du checkbox
        oDispOpt = set_DispOptChecked( oDispOpt, strDisplayOption, get( hChk, 'Value' ) );
    end
    
    %Determination de l'intervalle
    popupStrings = get(handles.popup_HoleIdInterval, 'String' );
    NumSelection = get( handles.popup_HoleIdInterval, 'Value' );
    strSelection = popupStrings{ NumSelection };
    
    %L'objet prend la valeur de l'intervalle
    oDispOpt = set_HoleIdInterval( oDispOpt, strread(strSelection) );
    
    %Ajout diamètre
    name = get(handles.string_diameter,'string');
    oDispOpt = set_SphereDiameter(oDispOpt,name);
    
    name = get(handles.edit_atlasnum,'string');
    oDispOpt = set_DispAtlasZone(oDispOpt,name);
%--------------------------------------------------------------------------
% --- Executes on selection change in popup_HoleIdInterval.
function popup_HoleIdInterval_Callback(hObject, eventdata, handles)
    % hObject    handle to popup_HoleIdInterval (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns popup_HoleIdInterval contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from popup_HoleIdInterval


%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function popup_HoleIdInterval_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to popup_HoleIdInterval (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%Digitalisation
function chk_Dig_HelmetHoles_Callback(                     hObject, eventdata, handles )
function chk_Dig_HelmetHolesNormals_Callback(              hObject, eventdata, handles )
function chk_Dig_HelmetRegistrationPoints_Callback(        hObject, eventdata, handles )
function chk_Dig_SubjectFiducials_Callback(                hObject, eventdata, handles )
function chk_Dig_SubjectRegistrationTestMarkers_Callback(  hObject, eventdata, handles )
function chk_Dig_CoordinateSystemAxes_Callback(            hObject, eventdata, handles )
function chk_SphereHoles_Callback(                         hObject, eventdata, handles )
function chk_line_holes_Callback(                          hObject, eventdata, handles )

%Labels
function chk_Lbl_HelmetHoleId_Callback(                    hObject, eventdata, handles )
function chk_Lbl_OptodesHoleId_Callback(                   hObject, eventdata, handles )
function chk_Lbl_OptodesFibersId_Callback(                 hObject, eventdata, handles )

%Appearance
function chk_App_BoldLabels_Callback(                      hObject, eventdata, handles )
function chk_App_InvertedColorScheme_Callback(             hObject, eventdata, handles )
function chk_App_BoldLines_Callback(                       hObject, eventdata, handles )
function chk_App_ItemsAlwaysVisible_Callback(              hObject, eventdata, handles )
function chk_App_Reflections_Callback(                     hObject, eventdata, handles )
function chk_App_Transparency_Callback(                    hObject, eventdata, handles )

%Optodes
function chk_Opt_FittedOnSkin_Callback(                    hObject, eventdata, handles )
function chk_Opt_CortexIntersectionMarker_Callback(        hObject, eventdata, handles )
function chk_Opt_MuxChannel_Callback(                      hObject, eventdata, handles )

%Electrodes
function chk_Ele_PositionOnHelmet_Callback(                hObject, eventdata, handles )

%MRI
function chk_MRI_SurfaceCortexHiRes_Callback(              hObject, eventdata, handles )
    if( get( hObject,'Value' ) )
        set( handles.chk_MRI_SurfaceCortexLowRes, 'Value', false )
        set( handles.chk_MRI_SurfaceSkin,         'Value', false )
    end
    
function chk_MRI_SurfaceCortexLowRes_Callback(             hObject, eventdata, handles )
    if( get( hObject,'Value' ) )
        set( handles.chk_MRI_SurfaceCortexHiRes,  'Value', false )
        set( handles.chk_MRI_SurfaceSkin,         'Value', false )
    end

function chk_MRI_SurfaceSkin_Callback(                     hObject, eventdata, handles )
    if( get( hObject,'Value' ) )
        set( handles.chk_MRI_SurfaceCortexHiRes,  'Value', false )
        set( handles.chk_MRI_SurfaceCortexLowRes, 'Value', false )
    end

function chk_MRI_ROIMarkers_Callback(                      hObject, eventdata, handles )
function chk_MRI_SurfaceHolesIntersectionMarkers_Callback( hObject, eventdata, handles )
function chk_MRI_SubjectFiducials_Callback(                hObject, eventdata, handles )







% --- Executes on button press in btn_MRI_Image.
function btn_MRI_Image_Callback(hObject, eventdata, handles)
% hObject    handle to btn_MRI_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name path] = uigetfile('.img');
set(handles.popup_MRIImage,'string',[path name])

% --- Executes on selection change in popup_MRIImage.
function popup_MRIImage_Callback(hObject, eventdata, handles)
% hObject    handle to popup_MRIImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_MRIImage contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_MRIImage


% --- Executes during object creation, after setting all properties.
function popup_MRIImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_MRIImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in chk_MRI_Image.
function chk_MRI_Image_Callback(hObject, eventdata, handles)
% hObject    handle to chk_MRI_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_MRI_Image




% --- Executes on button press in btn_clear.
function btn_clear_Callback(hObject, eventdata, handles)
% hObject    handle to btn_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.popup_MRIImage,'string','')



% --- Executes on button press in btn_Talairach.
function btn_Talairach_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Talairach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edit_atlasnum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_atlasnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_atlasnum as text
%        str2double(get(hObject,'String')) returns contents of edit_atlasnum as a double


% --- Executes during object creation, after setting all properties.
function edit_atlasnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_atlasnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in chk_SphereHolesDisplay.
function chk_SphereHolesDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to chk_SphereHolesDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_SphereHolesDisplay



function string_diameter_Callback(hObject, eventdata, handles)
% hObject    handle to string_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of string_diameter as text
%        str2double(get(hObject,'String')) returns contents of string_diameter as a double


% --- Executes during object creation, after setting all properties.
function string_diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to string_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in chk_App_BackgroundColor.
function chk_App_BackgroundColor_Callback(hObject, eventdata, handles)
% hObject    handle to chk_App_BackgroundColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_App_BackgroundColor
