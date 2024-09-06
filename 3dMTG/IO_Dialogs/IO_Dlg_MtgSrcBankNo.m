function varargout = IO_Dlg_MtgSrcBankNo(varargin)
% IO_DLG_MTGSRCBANKNO M-file for IO_Dlg_MtgSrcBankNo.fig
%      IO_DLG_MTGSRCBANKNO, by itself, creates a new IO_DLG_MTGSRCBANKNO or raises the existing
%      singleton*.
%
%      H = IO_DLG_MTGSRCBANKNO returns the handle to a new IO_DLG_MTGSRCBANKNO or the handle to
%      the existing singleton*.
%
%      IO_DLG_MTGSRCBANKNO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IO_DLG_MTGSRCBANKNO.M with the given input arguments.
%
%      IO_DLG_MTGSRCBANKNO('Property','Value',...) creates a new IO_DLG_MTGSRCBANKNO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IO_Dlg_MtgSrcBankNo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IO_Dlg_MtgSrcBankNo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IO_Dlg_MtgSrcBankNo

% Last Modified by GUIDE v2.5 29-Jun-2015 09:53:32

% Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @IO_Dlg_MtgSrcBankNo_OpeningFcn, ...
                       'gui_OutputFcn',  @IO_Dlg_MtgSrcBankNo_OutputFcn, ...
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


% --- Executes just before IO_Dlg_MtgSrcBankNo is made visible.
function IO_Dlg_MtgSrcBankNo_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to IO_Dlg_MtgSrcBankNo (see VARARGIN)

    if( ~numel( varargin ) || ~isa( varargin{1}, 'Helmet' ) )
        hDlg = errordlg( 'Error in IO_Dlg_MtgDetAlpha_OpeningFcn: varargin not a ''Helmet''. Loading Defaults' );
        uiwait(hDlg);
        oHelmet = Helmet;
    else
        oHelmet = varargin{1};
    end
    
    % Choose default command line output for IO_Dlg_DisplayOptions
    handles.output = oHelmet;
    
    %Sauvegarde de l'objet
    handles.m_oHelmet = oHelmet;
    
    handles.pItem = varargin{2};
    
    % Update handles structure
    guidata(hObject, handles);
    
    sMtg = get_Mtg( oHelmet );

    %Generation du string de popup
    [PhysicalSrcGroupSyncNo,PhysicalSrcCombinations] = get_PhysicalSrcGroups( oHelmet );
    
    if( sMtg.Gen_Params.Nb_Longueur_Onde == 1 )
        SrcCodes = PhysicalSrcCombinations;
    elseif( sMtg.Gen_Params.Nb_Longueur_Onde == 2 )
        SrcCodes = PhysicalSrcCombinations(:,1)+PhysicalSrcCombinations(:,2)*1000;
    else
        disp('Error: IOMtg cannot work with more than 2 wavelength');
    end
    
    popup_String(1,:) = sprintf( '%-10s', get_Mtg_nFib2strFib( oHelmet, SrcCodes(1) ) );
    for( i=1:get_Mtg_MaxSrcHoles(oHelmet) )
        popup_String(i,:) = sprintf( '%-10s', get_Mtg_nFib2strFib( oHelmet, SrcCodes(i) ) );
    end
    
    set( handles.popup_SrcFibers, 'String', popup_String );
    
    
    %Etablissement de la selection par defaut
    if( numel(sMtg.v_HolesMtg) >= handles.pItem )
        
        Selection = find( SrcCodes == sMtg.v_HolesMtg(handles.pItem) );
        
        if( numel(Selection) == 1 )
            %Valider que la selection est possible
            if( Selection >= 1 || Selection <= get_Mtg_MaxSrcHoles(oHelmet) )
                %Mettre la selection comme valeur initiale
                set( handles.popup_SrcFibers, 'Value', Selection );
            end
        end
    end

    % UIWAIT makes IO_Dlg_MtgSrcBankNo wait for user response (see UIRESUME)
    % uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IO_Dlg_MtgSrcBankNo_OutputFcn(hObject, eventdata, handles) 
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
function popup_SrcFibers_Callback(hObject, eventdata, handles)
    % Hints: contents = get(hObject,'String') returns popup_SrcFibers contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from popup_SrcFibers

%------------------------------------------------------------------
function popup_SrcFibers_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%------------------------------------------------------------------
function btn_Ok_Callback(hObject, eventdata, handles)
        
    oHelmet = handles.m_oHelmet;
    sMtg = get_Mtg( oHelmet );
    
    [PhysicalSrcGroupSyncNo,PhysicalSrcCombinations] = get_PhysicalSrcGroups( oHelmet );
    
    if( sMtg.Gen_Params.Nb_Longueur_Onde == 1 )
        SrcCodes = PhysicalSrcCombinations;
    elseif( sMtg.Gen_Params.Nb_Longueur_Onde == 2 )
        SrcCodes = PhysicalSrcCombinations(:,1)+PhysicalSrcCombinations(:,2)*1000;
    else
        disp('Error: IOMtg cannot work with more than 2 wavelength');
    end
    
    SrcFibersSelected = SrcCodes( get( handles.popup_SrcFibers, 'Value' ) );
    
    %Verifier si la fibre est deja utilisee (retirer les redondances)
    if( ~isempty( find( sMtg.v_HolesMtg == SrcFibersSelected ) ) )
        sMtg.v_HolesMtg( find( sMtg.v_HolesMtg == SrcFibersSelected ) ) = 0;
    end
    
    sMtg.v_HolesMtg( handles.pItem ) = SrcFibersSelected;
    
    %Verifier si le montage a du bon sens
    if( get_MtgSrc_Contaminated( set_Mtg(oHelmet, sMtg ) ) )
        h = warndlg( 'Contamination detected: this source cannot be used', 'Warning' );
        waitfor(h);
    else
        handles.output = set_Mtg(oHelmet, sMtg );
        guidata( handles.figure1, handles );
        uiresume( handles.figure1 );
    end

%------------------------------------------------------------------
function btn_Cancel_Callback(hObject, eventdata, handles)
  	uiresume( handles.figure1 );


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
