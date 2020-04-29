function varargout = warndlg_modal(varargin)
% WARNDLG_MODAL M-file for warndlg_modal.fig
%      WARNDLG_MODAL, by itself, creates a new WARNDLG_MODAL or raises the existing
%      singleton*.
%
%      H = WARNDLG_MODAL returns the handle to a new WARNDLG_MODAL or the handle to
%      the existing singleton*.
%
%      WARNDLG_MODAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WARNDLG_MODAL.M with the given input arguments.
%
%      WARNDLG_MODAL('Property','Value',...) creates a new WARNDLG_MODAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before warndlg_modal_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to warndlg_modal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help warndlg_modal

% Last Modified by GUIDE v2.5 24-Apr-2007 09:27:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @warndlg_modal_OpeningFcn, ...
                   'gui_OutputFcn',  @warndlg_modal_OutputFcn, ...
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


% --- Executes just before warndlg_modal is made visible.
function warndlg_modal_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to warndlg_modal (see VARARGIN)

    % Choose default command line output for warndlg_modal
    handles.output = hObject;

    if( numel( varargin ) >= 1 )
        set( handles.text_Message, 'String', varargin{1} );
        if( numel( varargin ) >= 2 )
            set( handles.warndlg_modal, 'Name', varargin{2} );
        end
    end

    axes(handles.axe_Icon);
    matpix = imread( 'warning_Icon.bmp' );
    image( matpix  ); 
    axis image;
    axis off;
    
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes warndlg_modal wait for user response (see UIRESUME)
    % uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = warndlg_modal_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

    uiwait( hObject );
    
    if( ishandle(hObject) )
        close(hObject);
    end

% --- Executes on button press in warndlg_modal.
%function warndlg_modal_Callback(hObject, eventdata, handles)
% hObject    handle to warndlg_modal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_ok.
function btn_ok_Callback(hObject, eventdata, handles)
    % hObject    handle to btn_ok (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    uiresume(handles.warndlg_modal);

