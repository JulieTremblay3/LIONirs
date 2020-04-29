function varargout = IO_Dlg_GenParams(varargin)
    % IO_DLG_GENPARAMS M-file for IO_Dlg_GenParams.fig
    %      IO_DLG_GENPARAMS, by itself, creates a new IO_DLG_GENPARAMS or raises the existing
    %      singleton*.
    %
    %      H = IO_DLG_GENPARAMS returns the handle to a new IO_DLG_GENPARAMS or the handle to
    %      the existing singleton*.
    %
    %      IO_DLG_GENPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in IO_DLG_GENPARAMS.M with the given input arguments.
    %
    %      IO_DLG_GENPARAMS('Property','Value',...) creates a new IO_DLG_GENPARAMS or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before IO_Dlg_GenParams_OpeningFunction gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to IO_Dlg_GenParams_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES


    % Edit the above text to modify the response to help IO_Dlg_GenParams

    % Last Modified by GUIDE v2.5 23-Jan-2013 15:48:06

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @IO_Dlg_GenParams_OpeningFcn, ...
                       'gui_OutputFcn',  @IO_Dlg_GenParams_OutputFcn, ...
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


%------------------------------------------------------------------
function IO_Dlg_GenParams_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to IO_Dlg_GenParams (see VARARGIN)

    if( ~numel( varargin ) || ~isa( varargin{1}, 'Helmet' ) )
        hDlg = errordlg( 'Error in IO_Dlg_GenParams_OpeningFcn: varargin not a ''Helmet''. Loading Defaults' );
        uiwait(hDlg);
        oHelmet = Helmet;
    else
        oHelmet = varargin{1};
    end
    
    % Choose default command line output for IO_Dlg_DisplayOptions
    handles.output = oHelmet;
    
    %Sauvegarde du parametre d'entree
    handles.m_oHelmet = oHelmet;
    
    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes IO_Dlg_GenParams wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
    InitControls(handles);

%------------------------------------------------------------------
function varargout = IO_Dlg_GenParams_OutputFcn(hObject, eventdata, handles) 
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

%------------------------------------------------------------------
function txt_MinDistMux_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------
function txt_MinDistMux_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%------------------------------------------------------------------
function txt_MaxDistMux_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------
function txt_MaxDistMux_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%------------------------------------------------------------------
function txt_DistContamination_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------
function txt_DistContamination_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%------------------------------------------------------------------
function txt_NbSrcDisp_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------
function txt_NbSrcDisp_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
%------------------------------------------------------------------
function txt_NbDetDisp_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------
function txt_NbDetDisp_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%------------------------------------------------------------------
function txt_NbSrcHole_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------
function txt_NbSrcHole_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%------------------------------------------------------------------
function btn_Cancel_Callback(hObject, eventdata, handles)
    
    uiresume(handles.figure1);

%------------------------------------------------------------------
function txt_NbSrcPerBank_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------
function txt_NbSrcPerBank_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%------------------------------------------------------------------
function txt_NbBanks_Callback(hObject, eventdata, handles)

%------------------------------------------------------------------
function txt_NbBanks_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
%------------------------------------------------------------------
function txt_DistCalcThickness_Callback(hObject, eventdata, handles)
    %handles.DistCalcThickness = str2double(get(hObject,'String'));
    
%------------------------------------------------------------------
function txt_DistCalcThickness_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

% -------------------------------------------------------------------------
function InitControls( handles )
    
    oHelmet = handles.m_oHelmet;
    sMtg = get_Mtg(oHelmet);
    
    set( handles.txt_MinDistMux,         'String', sprintf('%.1f', sMtg.Gen_Params.Min_Dist_Mux*100      ) );
    set( handles.txt_MaxDistMux,         'String', sprintf('%.1f', sMtg.Gen_Params.Max_Dist_Mux*100      ) );
    set( handles.txt_DistContamination,  'String', sprintf('%.1f', sMtg.Gen_Params.DistContamination*100 ) );
    set( handles.txt_NbSrcDisp,          'String', sprintf('%d', sMtg.Gen_Params.Nb_Sources_Disp       ) );
    set( handles.txt_NbDetDisp,          'String', sprintf('%d', sMtg.Gen_Params.Nb_Detect_Disp        ) );
    set( handles.txt_NbSrcHole,          'String', sprintf('%d', sMtg.Gen_Params.Nb_Longueur_Onde      ) );
    set( handles.txt_NbSrcPerBank,       'String', sprintf('%d', sMtg.Gen_Params.NbSrcPerBank          ) );
    set( handles.txt_DistCalcThickness,  'String', sprintf('%.1f', sMtg.Gen_Params.UniformSkinDepth*100    ) );
    set( handles.txt_NbBanks,            'String', sprintf('%.1f', sMtg.Gen_Params.NbBanks    ) );
    if ~isfield(sMtg.Gen_Params,'MontageAutomatic')
        sMtg.Gen_Params.MontageAutomatic = 1;
    end
    if ~isfield(sMtg.Gen_Params,'AcqSystem')
        sMtg.Gen_Params.AcqSystem = 'ISS';
    end
    if strcmp(sMtg.Gen_Params.AcqSystem,'ISS')
        set(handles.popup_AcqSystem,'value',1)
    elseif strcmp(sMtg.Gen_Params.AcqSystem,'IMAGINC')
        set(handles.popup_AcqSystem,'value',2)
    elseif strcmp(sMtg.Gen_Params.AcqSystem,'IMAGINCHUBERT')
        set(handles.popup_AcqSystem,'value',3)
    end
        
    set(handles.radio_MontageAutomatic,'value',sMtg.Gen_Params.MontageAutomatic);
    if( sMtg.Gen_Params.bUseUniformSkinDepth )
        set( handles.chk_DistUniform, 'Value', true );
    else
        set( handles.chk_DistOnSkin, 'Value', true );
        set( handles.txt_DistCalcThickness, 'Enable', 'off' );
    end    
    if strcmp(sMtg.Gen_Params.ElectrodeType,'10-10') 
        set(handles.radio_1010,'value',1) 
    elseif strcmp(sMtg.Gen_Params.ElectrodeType,'E128') 
         set(handles.radio_E128,'value',1) 
    elseif strcmp(sMtg.Gen_Params.ElectrodeType,'Custom') 
        set(handles.radio_custom,'value',1) 
    end
    guidata(handles.figure1, handles);
    
   

%------------------------------------------------------------------
function btn_Ok_Callback(hObject, eventdata, handles)

    oHelmet = handles.m_oHelmet;
    sMtg = get_Mtg(oHelmet);

    sOldMtg = sMtg;
    
    %Memorisation des parametres de generation de montage
    sMtg.Gen_Params.Min_Dist_Mux       = str2double(get(handles.txt_MinDistMux, 'String' ))/100;
    sMtg.Gen_Params.Max_Dist_Mux       = str2double(get(handles.txt_MaxDistMux, 'String' ))/100;
    sMtg.Gen_Params.DistContamination  = str2double(get(handles.txt_DistContamination, 'String' ))/100;
    sMtg.Gen_Params.Nb_Sources_Disp    = str2double(get(handles.txt_NbSrcDisp, 'String' ));
    sMtg.Gen_Params.Nb_Detect_Disp     = str2double(get(handles.txt_NbDetDisp, 'String' ));
    sMtg.Gen_Params.Nb_Longueur_Onde   = str2double(get(handles.txt_NbSrcHole, 'String' ));
    sMtg.Gen_Params.NbSrcPerBank       = str2double(get(handles.txt_NbSrcPerBank, 'String' ));
    sMtg.Gen_Params.UniformSkinDepth   = str2double(get(handles.txt_DistCalcThickness, 'String' ))/100;
    sMtg.Gen_Params.NbBanks            = str2double(get(handles.txt_NbBanks, 'String' ));
    sMtg.Gen_Params.bUseUniformSkinDepth = get( handles.chk_DistUniform, 'Value' );
    sMtg.Gen_Params.MontageAutomatic = get( handles.radio_MontageAutomatic, 'Value' );
    if get(handles.radio_1010,'value')
        sMtg.Gen_Params.ElectrodeType = '10-10';
    elseif get(handles.radio_E128,'value')
        sMtg.Gen_Params.ElectrodeType = 'E128';
    elseif  get(handles.radio_custom,'value')
        sMtg.Gen_Params.ElectrodeType = 'Custom';
    end
    
    
    id = get(handles.popup_AcqSystem,'value');
    AcqSystem = get(handles.popup_AcqSystem,'string');
    sMtg.Gen_Params.AcqSystem = AcqSystem{id};
    
    bMustResetMtg = false;
    
    if(    sMtg.Gen_Params.Nb_Sources_Disp < sOldMtg.Gen_Params.Nb_Sources_Disp ...
        || sMtg.Gen_Params.Nb_Detect_Disp  < numel(sOldMtg.v_pDet) ...
        || sMtg.Gen_Params.NbSrcPerBank    < sOldMtg.Gen_Params.NbSrcPerBank  ...
        || get_MtgSrc_Contaminated( set_Mtg( oHelmet, sMtg ) ) )
        bMustResetMtg = true;
    end
        
    if( bMustResetMtg )
        button = questdlg('Are you sure?','Warning','Yes','No','No');
%         button = questdlg('Parameters cannot be changed without clearing montage. Are you sure you want to clear current Montage?','Warning','Yes','No','No');
        if( strcmp( button, 'No' ) )
            return;
        end        
    end

    %Memorisation du montage dans le casque dans le projet
    oHelmet = set_Mtg( oHelmet, sMtg );    
    handles.m_oHelmet = oHelmet;    
    handles.output = oHelmet;
    guidata(handles.figure1, handles);
    uiresume(handles.figure1);

% --- Executes on button press in chk_DistOnSkin.
function chk_DistOnSkin_Callback(hObject, eventdata, handles)
% hObject    handle to chk_DistOnSkin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_DistOnSkin
if( get(hObject,'Value') )
    set( handles.chk_DistUniform, 'Value', false );
    set( handles.txt_DistCalcThickness, 'Enable', 'off' );
else
    set( handles.chk_DistUniform, 'Value', true );
    set( handles.txt_DistCalcThickness, 'Enable', 'on' );
end
    
    

% --- Executes on button press in chk_DistUniform.
function chk_DistUniform_Callback(hObject, eventdata, handles)
% hObject    handle to chk_DistUniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_DistUniform
if( get(hObject,'Value') )
    set( handles.chk_DistOnSkin, 'Value', false );
    set( handles.txt_DistCalcThickness, 'Enable', 'on' );
else
    set( handles.chk_DistOnSkin, 'Value', true );
    set( handles.txt_DistCalcThickness, 'Enable', 'off' );
end



% --- Executes on button press in radio_montageautomatic.
function radio_montageautomatic_Callback(hObject, eventdata, handles)
% hObject    handle to radio_montageautomatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_montageautomatic




% --- Executes on button press in radio_MontageAutomatic.
function radio_MontageAutomatic_Callback(hObject, eventdata, handles)
% hObject    handle to radio_MontageAutomatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_MontageAutomatic




% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radio_E128.
function radio_E128_Callback(hObject, eventdata, handles)
% hObject    handle to radio_E128 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_E128


% --- Executes on button press in radio_custom.
function radio_custom_Callback(hObject, eventdata, handles)
% hObject    handle to radio_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_custom




% --- Executes on selection change in popup_AcqSystem.
function popup_AcqSystem_Callback(hObject, eventdata, handles)
% hObject    handle to popup_AcqSystem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_AcqSystem contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_AcqSystem
if get(handles.popup_AcqSystem,'value')==1
    set(handles.txt_NbSrcDisp,'string','64')
    set(handles.txt_NbSrcDisp,'enable','on')
    set(handles.txt_NbDetDisp,'string','32')
    set(handles.txt_NbDetDisp,'enable','on')
    set(handles.txt_NbBanks,'string','4')
    set(handles.txt_NbBanks,'enable','on')
    set(handles.txt_NbSrcPerBank,'string','16')
    set(handles.txt_NbSrcPerBank,'enable','on')
    set(handles.txt_NbSrcHole,'string','2')
    set(handles.txt_NbSrcHole,'enable','on')
    
elseif get(handles.popup_AcqSystem,'value')==2
    set(handles.txt_NbSrcDisp,'string','64')
    set(handles.txt_NbSrcDisp,'enable','off')
    set(handles.txt_NbDetDisp,'string','32')
    set(handles.txt_NbDetDisp,'enable','off')
    set(handles.txt_NbBanks,'string','4')
    set(handles.txt_NbBanks,'enable','off')
    set(handles.txt_NbSrcPerBank,'string','16')
    set(handles.txt_NbSrcPerBank,'enable','off')
    set(handles.txt_NbSrcHole,'string','2')
    set(handles.txt_NbSrcHole,'enable','off')
    
end

% --- Executes during object creation, after setting all properties.
function popup_AcqSystem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_AcqSystem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


