function varargout = IO_Dlg_MtgElectrode(varargin)
% IO_DLG_MTGELECTRODE M-file for IO_Dlg_MtgElectrode.fig
%      IO_DLG_MTGELECTRODE, by itself, creates a new IO_DLG_MTGELECTRODE or raises the existing
%      singleton*.
%
%      H = IO_DLG_MTGELECTRODE returns the handle to a new IO_DLG_MTGELECTRODE or the handle to
%      the existing singleton*.
%
%      IO_DLG_MTGELECTRODE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IO_DLG_MTGELECTRODE.M with the given input arguments.
%
%      IO_DLG_MTGELECTRODE('Property','Value',...) creates a new IO_DLG_MTGELECTRODE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IO_Dlg_MtgElectrode_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IO_Dlg_MtgElectrode_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IO_Dlg_MtgElectrode

% Last Modified by GUIDE v2.5 23-Aug-2007 12:41:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IO_Dlg_MtgElectrode_OpeningFcn, ...
                   'gui_OutputFcn',  @IO_Dlg_MtgElectrode_OutputFcn, ...
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


% --- Executes just before IO_Dlg_MtgElectrode is made visible.
function IO_Dlg_MtgElectrode_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IO_Dlg_MtgElectrode (see VARARGIN)

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
    
    if strcmp(sMtg.Gen_Params.ElectrodeType,'10-10') 
        list = liste_electrode(1);
    elseif strcmp(sMtg.Gen_Params.ElectrodeType,'E128') 
        list = liste_electrode(2);
    elseif strcmp(sMtg.Gen_Params.ElectrodeType,'Custom') 
        list =liste_electrode(3);
        if ~iscell(list)
         list ={list};
        end
    end
    set( handles.popup_Electrode, 'String', list);
    if( numel(sMtg.v_HolesMtg) >= handles.pItem )
        Selection = sMtg.v_HolesEle(handles.pItem);
        if isempty(Selection)
            ind = strmatch( Selection,list);
        else
            ind = 1;
        end
        if isempty(ind)
            ind = 1;
        end
        %Mettre la selection comme valeur initiale
        set( handles.popup_Electrode, 'Value', ind );
    end

    % UIWAIT makes IO_Dlg_MtgDetAlpha wait for user response (see UIRESUME)
    % uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = IO_Dlg_MtgElectrode_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in btn_Ok.
function btn_Ok_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    oHelmet = handles.m_oHelmet;
    sMtg = get_Mtg( oHelmet );
    
    %Verifier la possibilite du montage
    EleFiberSelected = get( handles.popup_Electrode, 'Value' );
    if strcmp(sMtg.Gen_Params.ElectrodeType,'10-10') 
        list = liste_electrode(1);
    elseif strcmp(sMtg.Gen_Params.ElectrodeType,'E128') 
        list = liste_electrode(2);
    elseif strcmp(sMtg.Gen_Params.ElectrodeType,'Custom') 
        list =  get( handles.popup_Electrode,'string');
        sMtg.v_HolesEle( handles.pItem ) =  list;
        handles.output = set_Mtg(oHelmet, sMtg );
        guidata( handles.figure1, handles );
        uiresume( handles.figure1 );
        return
    end
    liste_ele = list;
    label_ele = list(EleFiberSelected );
    %Verifier si la fibre est deja utilisee (retirer les redondances)
    for i = 1:numel(sMtg.v_HolesEle)
        if iscellstr(sMtg.v_HolesEle(i))
            if strmatch(label_ele, sMtg.v_HolesEle(i)) 
                sMtg.v_HolesEle{i} = [];
            end
        end
    end
    sMtg.v_HolesEle( handles.pItem ) =  label_ele;
    handles.output = set_Mtg(oHelmet, sMtg );
    guidata( handles.figure1, handles );
    uiresume( handles.figure1 );
  
    
% --- Executes on button press in btn_cancel.
function btn_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to btn_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 uiresume( handles.figure1 );
 
% --- Executes on selection change in popup_Electrode.
function popup_Electrode_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_Electrode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Electrode


% --- Executes during object creation, after setting all properties.
function popup_Electrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


