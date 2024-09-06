function varargout = IO_Dlg_Export(varargin)
% IO_DLG_EXPORT M-file for IO_Dlg_Export.fig
%      IO_DLG_EXPORT, by itself, creates a new IO_DLG_EXPORT or raises the existing
%      singleton*.
%
%      H = IO_DLG_EXPORT returns the handle to a new IO_DLG_EXPORT or the handle to
%      the existing singleton*.
%
%      IO_DLG_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IO_DLG_EXPORT.M with the given input arguments.
%
%      IO_DLG_EXPORT('Property','Value',...) creates a new IO_DLG_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IO_Dlg_Export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IO_Dlg_Export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IO_Dlg_Export

% Last Modified by GUIDE v2.5 27-Oct-2014 10:49:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IO_Dlg_Export_OpeningFcn, ...
                   'gui_OutputFcn',  @IO_Dlg_Export_OutputFcn, ...
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


% --- Executes just before IO_Dlg_Export is made visible.
function IO_Dlg_Export_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to IO_Dlg_Export (see VARARGIN)

    % Choose default command line output for IO_Dlg_Export
    handles.output = hObject;

    handles.m_oPrjData = varargin{1};
    
    % Update handles structure
    guidata(hObject, handles);
    
    axes(handles.axe_Export_Mtg_Icon);
    matpix = imread( 'Export_Mtg_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_Export_MtgSEQ_Icon);
    matpix = imread( 'Export_MtgSEQ_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_Export_MtgELP_Icon);
    matpix = imread( 'Export_MtgELP_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_Export_CompleteELP_Icon);
    matpix = imread( 'Export_CompleteELP_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_Export_BrainsightTool_Icon);
    matpix = imread( 'Export_BrainsightTool_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    axes(handles.axe_Export_iMagicELE_Icon);
    matpix = imread( 'Export_iMagicELE_Icon.bmp' );
    image( matpix  ); 
    axis image;
	axis off;
    
    
% UIWAIT makes IO_Dlg_Export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IO_Dlg_Export_OutputFcn(hObject, eventdata, handles) 

    % Get default command line output from handles structure
    varargout{1} = handles.output;    
    uiwait( hObject );    
    if( ishandle(hObject) )
        varargout{1} = handles.output;
        close( hObject );
    end

function radio_Export_MtgSEQ_Callback(hObject, eventdata, handles)
    SelectRadio(hObject);
    
function radio_Export_MtgELP_Callback(hObject, eventdata, handles)
    SelectRadio(hObject);
    
function radio_Export_BrainsightTool_Callback(hObject, eventdata, handles)
    SelectRadio(hObject);

function radio_Export_iMagicELE_Callback(hObject, eventdata, handles)
    SelectRadio(hObject);

function radio_Export_Mtg_Callback(hObject, eventdata, handles)
   SelectRadio(hObject);

function radio_Export_LocELP_Callback(hObject, eventdata, handles)
    SelectRadio(hObject);

function SelectRadio( hObject )
    set( findobj('Tag','radio_Export_MtgSEQ'), 'Value', false );
    set( findobj('Tag','radio_Export_MtgELP'), 'Value', false );
    set( findobj('Tag','radio_Export_LocELP'), 'Value', false );
    set( findobj('Tag','radio_Export_BrainsightTool'), 'Value', false );
    set( findobj('Tag','radio_Export_iMagicELE'), 'Value', false );
    set( findobj('Tag','radio_Export_Mtg'),'value',false);
    set( hObject, 'Value', true );
  
function chk_SrcDetOnly_Callback(hObject, eventdata, handles)

function btn_Cancel_Callback(hObject, eventdata, handles)
    uiresume( handles.figure_Export );

function btn_Ok_Callback(hObject, eventdata, handles)
    oDigCompleteHelmet = get_Dig_CompleteHelmet(handles.m_oPrjData);
    oHelmet            = get_Helmet(handles.m_oPrjData);
    oMRI_Data          = get_MRI_Data(handles.m_oPrjData);
    if (get(handles.radio_Export_Mtg,'value'))
        [FileName,PathName] = uiputfile({'*.mtg','Reuse montage to create an other project(*.mtg)'},'Save As...');
                if( isempty(strfind(FileName, '.mtg')) )
                FileName = [FileName, '.mtg'];
                end
            save( [PathName,FileName ],'oHelmet');
    end

    if(     get( handles.radio_Export_MtgSEQ, 'Value' ) )
        [FileName,PathName] = uiputfile({'*.seq','Brainsight Sequence File (*.seq)'},'Save As...');
        if( ~isempty(FileName) )
            if( isempty(strfind(FileName, '.seq')) )
                FileName = [FileName, '.seq'];
            end
            Write_Brainsight_SEQ( oHelmet, fullfile(PathName,FileName) );
        else
            return;
        end
    elseif( get( handles.radio_Export_MtgELP, 'Value' ) )
        [FileName,PathName] = uiputfile({'*.elp','Opt3D .elp File (*.elp)'},'Save As...');
        if( ~isempty(FileName) )
                  %MATTRANSFORM = 
            IRMAdjustment.matTransformManual = get_matTransformManual(oMRI_Data)
            IRMAdjustment.matTransform = get_matTransform(oMRI_Data)
            save(fullfile(PathName,[FileName,'.mat']),'IRMAdjustment')
           
            if( isempty(strfind(FileName, '.elp')) ) 
                FileName = [FileName, '.elp'];
            end        
            Write_OPT3D_Elp( oHelmet, fullfile(PathName,FileName),get(handles.radio_allholeelp,'value'), get(handles.radio_fitonskin,'value'),get(handles.radio_MRIfid,'value'));
      
        else
            return;
        end
    elseif( get( handles.radio_Export_LocELP, 'Value' ) )
        if( isempty( oDigCompleteHelmet ) )
            warndlg_modal( 'Cannot create file: Complete helmet LOCATOR .elp file not found in project' );
            return;
        else
            [FileName,PathName] = uiputfile({'*.elp','LOCATOR File (*.elp)'},'Save As...');
            if( ~isempty(FileName) )
                if( isempty(strfind(FileName, '.elp')) )
                    FileName = [FileName, '.elp'];
                end
                Write_LOCATOR_Elp( oDigCompleteHelmet, fullfile(PathName,FileName) );
            else
                return;
            end
        end
        
    elseif( get( handles.radio_Export_BrainsightTool, 'Value' ) )
        if( isempty( oDigCompleteHelmet ) )
            warndlg_modal( 'Cannot create file: Complete helmet LOCATOR .elp file not found in project' );
            return;
        else
            [FileName,PathName] = uiputfile({'*.txt','BrainSight Tool File (*.txt)'},'Save As...');
            bExportSrcDetOnly = get( handles.chk_SrcDetOnly, 'Value' );
            if( ~isempty(FileName) )
                if( isempty(strfind(FileName, '.txt')) )
                    FileName = [FileName, '.txt'];
                end
                Write_Brainsight_Tool( oHelmet, oDigCompleteHelmet, bExportSrcDetOnly, fullfile(PathName,FileName) );
            else
                return;
            end
        end
        
    elseif( get( handles.radio_Export_iMagicELE, 'Value' ) )
        if( isempty( find( get_matFiducials( oMRI_Data ) ) ) )
            warndlg_modal( 'Cannot create .ELE file: No fiducial markers found in MRI Data' );
            return;
        else
            [FileName,PathName] = uiputfile({'*.ele','iMagic Electrode File (*.ele)'},'Save As...');
            if( ~isempty(FileName) )
                if( isempty(strfind(FileName, '.ele')) )
                    FileName = [FileName, '.ele'];
                end
                if get(handles.radio_Helmet,'value');projection=1;
                elseif get(handles.radio_Skin,'value');projection=2;
                elseif get(handles.radio_Cortex,'value');projection = 3;
                end
                Write_iMagic_ELE( oHelmet, oMRI_Data, fullfile(PathName,FileName),projection );
            else
                return;
            end
        end
    end
    uiresume( handles.figure_Export );


% --- Executes on button press in radiobutton10.
function radiobutton10_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton10




% --- Executes on button press in radio_allholeelp.
function radio_allholeelp_Callback(hObject, eventdata, handles)
% hObject    handle to radio_allholeelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_allholeelp




% --- Executes on button press in radio_fitonskin.
function radio_fitonskin_Callback(hObject, eventdata, handles)
% hObject    handle to radio_fitonskin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_fitonskin




% --- Executes on button press in radio_MRIfid.
function radio_MRIfid_Callback(hObject, eventdata, handles)
% hObject    handle to radio_MRIfid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_MRIfid


