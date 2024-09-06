function varargout = IO_Dlg_MtgDetAlpha(varargin)
% IO_DLG_MTGDETALPHA M-file for IO_Dlg_MtgDetAlpha.fig
%      IO_DLG_MTGDETALPHA, by itself, creates a new IO_DLG_MTGDETALPHA or raises the existing
%      singleton*.
%
%      H = IO_DLG_MTGDETALPHA returns the handle to a new IO_DLG_MTGDETALPHA or the handle to
%      the existing singleton*.
%
%      IO_DLG_MTGDETALPHA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IO_DLG_MTGDETALPHA.M with the given input arguments.
%
%      IO_DLG_MTGDETALPHA('Property','Value',...) creates a new IO_DLG_MTGDETALPHA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IO_Dlg_MtgDetAlpha_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IO_Dlg_MtgDetAlpha_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IO_Dlg_MtgDetAlpha

% Last Modified by GUIDE v2.5 19-Jul-2006 16:24:50

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @IO_Dlg_MtgDetAlpha_OpeningFcn, ...
                       'gui_OutputFcn',  @IO_Dlg_MtgDetAlpha_OutputFcn, ...
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


% --- Executes just before IO_Dlg_MtgDetAlpha is made visible.
function IO_Dlg_MtgDetAlpha_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to IO_Dlg_MtgDetAlpha (see
    % VARARGIN)
    
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
    if ~isfield(sMtg.Gen_Params,'AcqSystem')
         sMtg.Gen_Params.AcqSystem = 'ISS'
    end
    if strcmp(sMtg.Gen_Params.AcqSystem, 'ISS')
        popup_String(1,1) = 'A';
        for( i=1:sMtg.Gen_Params.Nb_Detect_Disp )
            popup_String(i,1) = char( int16('A')+(i-1) );
        end
    else
         popup_String=[] ;
        for i=1:sMtg.Gen_Params.Nb_Detect_Disp 
            popup_String = [popup_String;'D',sprintf('%02.0f',i)];%,str2num(i);
        end
    end 
    set( handles.popup_DetAlpha, 'String', popup_String );
    
    if( numel(sMtg.v_HolesMtg) >= handles.pItem )
        Selection = sMtg.v_HolesMtg(handles.pItem)/1000000;
        
        %Valider que la selection est possible
        if( Selection < 1 || Selection > sMtg.Gen_Params.Nb_Detect_Disp )
            Selection = 1;
        end
        
        %Mettre la selection comme valeur initiale
        set( handles.popup_DetAlpha, 'Value', Selection );
    end

    % UIWAIT makes IO_Dlg_MtgDetAlpha wait for user response (see UIRESUME)
    % uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IO_Dlg_MtgDetAlpha_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on selection change in popup_DetAlpha.
function popup_DetAlpha_Callback(hObject, eventdata, handles)
    % hObject    handle to popup_DetAlpha (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns popup_DetAlpha contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from popup_DetAlpha


%------------------------------------------------------------------
function popup_DetAlpha_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to popup_DetAlpha (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%------------------------------------------------------------------
function btn_Ok_Callback(hObject, eventdata, handles)

    oHelmet = handles.m_oHelmet;
    sMtg = get_Mtg( oHelmet );
    
    %Verifier la possibilite du montage
    DetFiberSelected = get( handles.popup_DetAlpha, 'Value' )*1000000;
    
    %Verifier si la fibre est deja utilisee (retirer les redondances)
    if( ~isempty( find( sMtg.v_HolesMtg == DetFiberSelected ) ) )
        sMtg.v_HolesMtg( find( sMtg.v_HolesMtg == DetFiberSelected ) ) = 0;
    end
        
    sMtg.v_HolesMtg( handles.pItem ) = DetFiberSelected;
    
    %Verifier si le montage a du bon sens
    if( get_MtgSrc_Contaminated( set_Mtg(oHelmet, sMtg ) ) )
        h = warndlg( 'This Detector cannot be used: too many sources are placed in proximity', 'Warning' );
        waitfor(h);
    else
        handles.output = set_Mtg(oHelmet, sMtg );
        guidata( handles.figure1, handles );
        uiresume( handles.figure1 );
    end
    


%------------------------------------------------------------------
function btn_Cancel_Callback(hObject, eventdata, handles)
    uiresume( handles.figure1 );
