%Description: fonction principale de l'intarface du casque.
% Communique avec Homer pour les information du .DOT
% Utilise un .prj prédéfini ou à définir
% Actualise l'affichage des graphiques d'Homer et GLM lors d'un changement
function varargout = IO_HelmetMTG_Display(varargin)
% IO_HELMETMTG_DISPLAY M-file for IO_HelmetMTG_Display.fig
% Fonction permettant d'afficher le casque avec homer
% On l'utilise pour selectionner les canaux desirés
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @IO_HelmetMTG_Display_OpeningFcn, ...
    'gui_OutputFcn',  @IO_HelmetMTG_Display_OutputFcn, ...
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


% --- Executes just before IO_HelmetMTG_Display is made visible.
function IO_HelmetMTG_Display_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IO_HelmetMTG_Display (see VARARGIN)

% Choose default command line output for IO_HelmetMTG_Display

setappdata(0,'guiHelmet',handles.IO_HelmetMTG); 
setappdata(0,'guiExportFigureOption',[]);
setappdata(handles.IO_HelmetMTG,'fhresetview',@resetview);
set(handles.IO_HelmetMTG,'CloseRequestFcn',@IO_Helmetclosereq);
handles.output = hObject;
%if homer or SPM gui or SPM video call the program
if numel(varargin)==2
   option = varargin{2};
else
   option = [];
end

figure(handles.IO_HelmetMTG)
if ~isfield(handles,'radio_guiSPMnirsHSJ')
        handles.radio_guiSPMnirsHSJ = uicontrol('style','popupmenu','string',[{'Homer'},{'SPMGui'},{'SPMvideo'}],'tag','radio_guiSPMnirsHSJ','unit','normalized','position',[0.05,0.1,0.12,0.03]);
end
if ~isempty(getappdata(0,'guiHOMER'));
    set(handles.radio_guiSPMnirsHSJ,'visible','on','value',1); %Homer
end
if isfield(option,'SPMgui')
     set(handles.radio_guiSPMnirsHSJ,'visible','on','value',2); %Gui SPM
     set(handles.popupmenu_display,'string',{'Current module','Mean start stop time','Projection Channel .cp','PARAFAC','PCA','GLM','COMPONENT'});
end
if isfield(option,'video')
      set(handles.radio_guiSPMnirsHSJ,'visible','on','value',3); %Gui video
end

set(handles.radio_guiSPMnirsHSJ, 'visible', 'off'); % Set the view selector to invisible as it is not needed.
if isfield(option,'restorestetting')
else%~isempty(varargin)
    set(handles.radio_Channel,'value',0);
    set(handles.radio_holes,'value',0);
    set(handles.radio_ViewMeasListAct,'value',0);
    set(handles.radio_viewskin,'value',1);
    if isfield(option,'scalemin')
        set(handles.edit_climmin,'string',option.scalemin);
    else
        set(handles.edit_climmin,'string','-1');
    end
    if isfield(option,'scalemax')
        set(handles.edit_climmax,'string',option.scalemax);
    else
        set(handles.edit_climmax,'string','1');
    end
    if isfield(option,'scaletr')
        if option.scaletr < option.scalemax*0.5
            set(handles.edit_cminthreshold,'string',option.scaletr)
        else
            set(handles.edit_cminthreshold,'string',option.scalemax/10);
        end
    else
       set(handles.edit_cminthreshold,'string',0.1);
    end
    if isfield(option,'cover') %show the uncover region in a different color
      set(handles.radio_show_cover,'value',option.cover)
    end
    if isfield(option,'skintype')
        if option.skintype==0
            set(handles.radio_viewskin,'value',1)
            set(handles.radio_viewcortex,'value',0)
        else
            set(handles.radio_viewskin,'value',0)
            set(handles.radio_viewcortex,'value',1)
        end
    else
        set(handles.radio_viewskin,'value',1)
        set(handles.radio_viewcortex,'value',0)
    end

end


handles.selectsrs = 0; % aucune source est utilisé
handles.selectsrsMarker = 0;
cmin = str2num(get(handles.edit_climmin,'string'));
cmax = str2num(get(handles.edit_climmax,'string'));
caxis([cmin,cmax]);
axes(handles.axes_Mtg2);
c = colorbar;
set(c,'fontsize',14)
guidata(hObject, handles);
global currentsub
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
    try %Vérifier l'ouverture du projet
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    catch
        return
    end
dcm_obj = datacursormode(handles.IO_HelmetMTG);
set(dcm_obj,'UpdateFcn',@myupdatefcn)
set(dcm_obj,'SnapToDataVertex','on')
set(dcm_obj,'DisplayStyle','window')

try %Vérifier la validiter du casque de données
    PrjData_file = PMI{currentsub}.prj_name;
    set(handles.IO_HelmetMTG,'name',['IO_Helmet  : ' , PrjData_file])
    LoadedStruct = load(PrjData_file,'-mat');
    PrjStruct = Load_PrjStruct(LoadedStruct,false);
catch %Proposer d'ouvrir un autre casque
    [filename, filepath] = uigetfile('.prj','Please select a project');
    PrjData_file = [filepath filename];
    if filename == 0
        close(gcf)
        return
    end
    LoadedStruct = load(PrjData_file,'-mat');
    PrjStruct = Load_PrjStruct(LoadedStruct,false);
    PMI{currentsub}.prj_name = PrjData_file;
    set(guiHOMER,'UserData',PMI);
end
setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
guidata(handles.IO_HelmetMTG, handles);
cameratoolbar(handles.IO_HelmetMTG);
resetview(handles);

p = mfilename('fullpath');
[tok,rem]= strtok(fliplr(p),'\');
path = fliplr(rem);
strBmp =[path,'AngleRotationCote.jpg'];
matpix = imread( strBmp,'jpg' );
set(handles.text_angle_rotation_cote, 'CData', matpix );
strBmp =[path,'AngleRotationFace.jpg'];
matpix = imread( strBmp,'jpg' );
set(handles.text_angle_rotation_face, 'CData', matpix );
strBmp =[path,'Play.jpg'];
matpix = imread( strBmp,'jpg' );
set(handles.btn_play, 'CData', matpix );
strBmp =[path,'Pause.jpg'];
matpix = imread( strBmp,'jpg' );
set(handles.btn_stop, 'CData', matpix );
colormap jet;


% --- Outputs from this function are returned to the command line.
function varargout = IO_HelmetMTG_Display_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function menu_IO_Openprj_Callback(hObject, eventdata, handles)
% Description : Verifier la validiter du projet a ouvrir et offrir a
% l'utilisateur d'ouvrir un nouveau projet.
% hObject    handle to menu_IO_Openprj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentsub

if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
try %Vérifier l'ouverture du projet
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
catch
    return
end

[filename, filepath] = uigetfile('.prj','Please select a project');
PrjData_file = [filepath filename];
if filename == 0
    return
end
LoadedStruct = load(PrjData_file,'-mat');
PrjStruct = Load_PrjStruct(LoadedStruct,false);
PMI{currentsub}.prj_name = PrjData_file;
set(handles.IO_HelmetMTG,'name',['IO_Helmet  : ' , PrjData_file]);
set(guiHOMER,'UserData',PMI);
setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
guidata(handles.IO_HelmetMTG, handles);
resetview(handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function IO_HelmetMTG_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to IO_HelmetMTG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
persistent currentsrs
if isempty(currentsrs)
    currentsrs.SDPairs = 0;
    currentsrs.type = 0;
end
global currentsub;
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
%HOMERhandles = guihandles(guiHOMER);
HOMERhandles = guidata(guiHOMER);
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
MeasList = PMI{currentsub}.data(cf).MeasList;
idxLambda = 1; %get(HOMERhandles.HOMer_popupmenu_Display2,'value');
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
DispHelm=get_Helmet(PrjStruct);
vHoles = get_vHoles( DispHelm );
sMtg = get_Mtg( DispHelm );

if size(PMI{currentsub}.plotLst,1)>1
    PMI{currentsub}.plotLst = PMI{currentsub}.plotLst';
end
% selection d'un canal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.radio_select_one_pairs,'value');
    if currentsrs.type == 0 %si on a rien on prend une source ou un detecteur
        pClosest = get_SrsDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint' ), get( handles.axes_Mtg2,  'CameraPosition' ));
        if ~(pClosest == 0)
            DetSrs = sMtg.v_HolesMtg(pClosest);
            [idxMin,is_srs] = Srs_n2SDpairs(DetSrs);
            currentsrs.SDPairs = idxMin;
            currentsrs.type = is_srs;
            plot3(vHoles(pClosest).Coord.x,vHoles(pClosest).Coord.y,vHoles(pClosest).Coord.z,'+','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',3);
            DispParameter.reset = 0;
        end
    elseif currentsrs.type == 1 %si on a un source on cherche un detecteur pres et on cree une paire
        Det_list = MeasList([find( MeasList(:,1)==currentsrs.SDPairs & MeasList(:,4)==idxLambda)],2);
        if ~isempty(Det_list)
            pClosest = get_pDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint'), Det_list);
            if ~(pClosest == 0)
                det = sMtg.v_HolesMtg(pClosest)/1000000;
                plotLstc = [find( MeasList(:,1)==currentsrs.SDPairs & MeasList(:,2)== det & MeasList(:,4)==idxLambda)];
                PMI{currentsub}.plotLst = plotLstc ;
                PMI{currentsub}.plot = [MeasList(plotLstc,1),MeasList(plotLstc,2)];
            end
        end
        currentsrs.type = 0;
        currentsrs.SDPairs = 0 ;
        DispParameter.reset = 1;
    elseif currentsrs.type == 2 %si on a un detecteur on cherche une source pres et on cree une paire
        Srs_list = MeasList([find( MeasList(:,2)==currentsrs.SDPairs & MeasList(:,4)==idxLambda)],1);
        pClosest = get_pSrsClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint'), Srs_list);
        if  ~(pClosest == 0)
            [srs,is_srs] = Srs_n2SDpairs(sMtg.v_HolesMtg(pClosest));
            plotLstc = [find( MeasList(:,1)== srs & MeasList(:,2)== currentsrs.SDPairs & MeasList(:,4)==idxLambda)];
            PMI{currentsub}.plotLst = plotLstc ;
            PMI{currentsub}.plot = [MeasList(plotLstc,1),MeasList(plotLstc,2)];
        end
        currentsrs.type = 0;
        currentsrs.SDPairs = 0 ;
        DispParameter.reset = 1;
    end
end
% selection d'un mux channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.radio_select_mux,'value');
    pClosest = get_SrsDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint' ), get( handles.axes_Mtg2,  'CameraPosition' ));
    if ~(pClosest == 0)
        DetSrs = sMtg.v_HolesMtg(pClosest);
        [idxMin,is_srs] = Srs_n2SDpairs(DetSrs);
        if is_srs == 1;
            lst = find( MeasList(:,1)==idxMin & MeasList(:,4)==idxLambda );
            PMI{currentsub}.plotLst = lst;
            PMI{currentsub}.plot = [idxMin*ones(length(lst),1) MeasList(lst,2)];
        elseif is_srs == 2;
            lst = find( MeasList(:,2)==idxMin & MeasList(:,4)==idxLambda );
            PMI{currentsub}.plotLst = lst;
            PMI{currentsub}.plot = [MeasList(lst,1) idxMin*ones(length(lst),1)];
        end
    end
    currentsrs.type = 0;
    currentsrs.SDPairs = 0 ;
    DispParameter.reset = 1;
end

if get(handles.radio_restore_mux,'value');
    pClosest = get_SrsDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint' ), get( handles.axes_Mtg2,  'CameraPosition' ));
    if ~(pClosest == 0)
        DetSrs = sMtg.v_HolesMtg(pClosest);
        [idxMin,is_srs] = Srs_n2SDpairs(DetSrs);
        if is_srs == 1;
            lst = find( MeasList(:,1)==idxMin)% & MeasList(:,4)==idxLambda );

        elseif is_srs == 2;
            lst = find( MeasList(:,2)==idxMin)% & MeasList(:,4)==idxLambda );

        end

        PMI{currentsub}.data(cf).MeasListAct(lst) = 1;
    end
    currentsrs.type = 0;
    currentsrs.SDPairs = 0 ;
    DispParameter.reset = 1;
end

% Dot out the selected mux.
if get(handles.radio_remove_mux,'value');
    pClosest = get_SrsDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint' ), get( handles.axes_Mtg2,  'CameraPosition' ));
    if ~(pClosest == 0)
        DetSrs = sMtg.v_HolesMtg(pClosest);
        [idxMin,is_srs] = Srs_n2SDpairs(DetSrs);
        if is_srs == 1;
            lst = find( MeasList(:,1)==idxMin & MeasList(:,4)==idxLambda );
            %PMI{currentsub}.plotLst = lst;
            %PMI{currentsub}.plot = [idxMin*ones(length(lst),1) MeasList(lst,2)];
        elseif is_srs == 2;
            lst = find( MeasList(:,2)==idxMin & MeasList(:,4)==idxLambda );
            %PMI{currentsub}.plotLst = lst;
            %PMI{currentsub}.plot = [MeasList(lst,1) idxMin*ones(length(lst),1)];
        end
        PMI{currentsub}.data(cf).MeasListAct(lst) = 0;
    end
    currentsrs.type = 0;
    currentsrs.SDPairs = 0 ;
    DispParameter.reset = 1;
end

% Remove the mux from the selection.
if get(handles.radio_delete_mux, 'value')
    pClosest = get_SrsDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint' ),get( handles.axes_Mtg2,  'CameraPosition' ) );
    if ~(pClosest == 0)
        DetSrs = sMtg.v_HolesMtg(pClosest);
        [idxMin,is_srs] = Srs_n2SDpairs(DetSrs);
        assert( (is_srs == 1 || is_srs == 2 || is_srs == 0), 'Error. The node type is not valid.'); % Verify if the node has a valid type.
        lst = []; % Empty the list.
        if is_srs == 1 % If source.
            lst = find( MeasList(:,1)==idxMin & MeasList(:,4)==idxLambda )'; % Get the list of channels next to the selected source.
        elseif is_srs == 2 % If detector.
            lst = find( MeasList(:,2)==idxMin & MeasList(:,4)==idxLambda )'; % Get the list of channels next to the selected detector.
        end
        if PMI{currentsub}.plot(1)==0
            PMI{currentsub}.plot=[];
            PMI{currentsub}.plotLst=[];
        end
        lst = setdiff(PMI{currentsub}.plotLst, lst); % Substract the list of currently selected channels with the list of tool-selected channels.
        PMI{currentsub}.plotLst = lst;
        PMI{currentsub}.plot = [MeasList(lst,1) MeasList(lst,2)];
        if length(PMI{currentsub}.plotLst) < 1 % If there is no channel left...
            PMI{currentsub}.plotLst=[1]; % Add a default one to avoid an exception.
            PMI{currentsub}.plot=[0,0];
        end
    end
    currentsrs.type = 0;
    currentsrs.SDPairs = 0 ;
    DispParameter.reset = 1;
end

% selection d'un plusieurs mux channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.radio_manymuxs,'value');
    pClosest = get_SrsDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint' ),get( handles.axes_Mtg2,  'CameraPosition' ) );
    if ~(pClosest == 0)
        DetSrs = sMtg.v_HolesMtg(pClosest);
        [idxMin,is_srs] = Srs_n2SDpairs(DetSrs);
        if is_srs == 1;
            lst = find( MeasList(:,1)==idxMin & MeasList(:,4)==idxLambda )';
            ind = [];
            for i = 1:numel(lst)
                indbad = find(lst(i)==PMI{currentsub}.plotLst);
                if ~isempty(indbad)
                    ind = [ind,i];
                end
            end
            if ~isempty(ind)
                lst(ind) = [];
            end
            if PMI{currentsub}.plot(1)==0
                PMI{currentsub}.plot=[];
                PMI{currentsub}.plotLst=[];
            end
            PMI{currentsub}.plotLst = [PMI{currentsub}.plotLst, lst];
            plot_lst = PMI{currentsub}.plot
            PMI{currentsub}.plot = [ plot_lst ; [idxMin*ones(length(lst),1) MeasList(lst,2)]];
        elseif is_srs == 2;
            lst = find( MeasList(:,2)==idxMin & MeasList(:,4)==idxLambda )';
            ind = [];
            for i = 1:numel(lst)
                indbad = find(lst(i)==PMI{currentsub}.plotLst);
                if ~isempty(indbad)
                    ind = [ind,i];
                end
            end
            if ~isempty(ind)
                lst(ind) = [];
            end
            if PMI{currentsub}.plot(1)==0
                PMI{currentsub}.plot=[];
                PMI{currentsub}.plotLst=[];
            end
            PMI{currentsub}.plotLst =  [PMI{currentsub}.plotLst,lst];
            plot_lst = PMI{currentsub}.plot;
            PMI{currentsub}.plot = [ plot_lst ;[  MeasList(lst,1) idxMin*ones(length(lst),1)]];
        end
    end
    currentsrs.type = 0;
    currentsrs.SDPairs = 0 ;
    DispParameter.reset = 1;
end

% selection de plusieurs canaux % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.radio_multiples_pairs,'value');
    if currentsrs.type == 0 %si on a rien on prend une source ou un detecteur
        pClosest = get_SrsDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint' ),get( handles.axes_Mtg2,  'CameraPosition' ) );
        if ~(pClosest == 0)
            DetSrs = sMtg.v_HolesMtg(pClosest);
            [idxMin,is_srs] = Srs_n2SDpairs(DetSrs);
            currentsrs.SDPairs = idxMin;
            currentsrs.type = is_srs;
            plot3(vHoles(pClosest).Coord.x,vHoles(pClosest).Coord.y,vHoles(pClosest).Coord.z,'+','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',3);
            DispParameter.reset = 0;
        end
    elseif currentsrs.type == 1 %si on a un source on cherche un detecteur pres et on cree une paire
        Det_list = MeasList([find( MeasList(:,1)==currentsrs.SDPairs & MeasList(:,4)==idxLambda)],2);
        pClosest = get_pDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint'), Det_list);
        if ~(pClosest == 0)
            det = sMtg.v_HolesMtg(pClosest)/1000000;
            plotLst = PMI{currentsub}.plotLst;
            plotdata  = PMI{currentsub}.plot;
            if plotdata == 0
                plotLst = [];
                plotdata = [];
            end
            plotLstc = find( MeasList(:,1)== currentsrs.SDPairs & MeasList(:,2)== det & MeasList(:,4)==idxLambda);
            if ~isempty(plotLstc)
                if isempty(find(plotLstc==plotLst))
                    plotLst = [plotLst ,  plotLstc ];
                    PMI{currentsub}.plotLst = plotLst;
                    PMI{currentsub}.plot = [MeasList(plotLst,1),MeasList(plotLst,2)];
                end
            end
        end
        currentsrs.type = 0;
        currentsrs.SDPairs = 0 ;
        DispParameter.reset = 1;
    elseif currentsrs.type == 2 %si on a un detecteur on cherche une source pres et on cree une paire
        Srs_list = MeasList([find( MeasList(:,2)==currentsrs.SDPairs & MeasList(:,4)==idxLambda)],1);
        pClosest = get_pSrsClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint'), Srs_list);
        if ~(pClosest == 0)
            [srs,is_srs] = Srs_n2SDpairs(sMtg.v_HolesMtg(pClosest));
            plotLst = PMI{currentsub}.plotLst;
            plotdata  = PMI{currentsub}.plot;
            if plotdata == 0
                plotLst = [];
                plotdata = [];
            end
            plotLstc = (find( MeasList(:,1)== srs & MeasList(:,2)== currentsrs.SDPairs & MeasList(:,4)==idxLambda));
            if ~isempty(plotLstc)
                if isempty(find(plotLstc==plotLst))
                    plotLst = [plotLst ,  plotLstc ];
                    PMI{currentsub}.plotLst = plotLst ;
                    PMI{currentsub}.plot = [MeasList(plotLst,1),MeasList(plotLst,2)];
                end
            end
        end
        currentsrs.type = 0;
        currentsrs.SDPairs = 0 ;
        DispParameter.reset = 1;
    end
end
% selection de tous les canaux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.radio_selectall,'value')
    lst = find(MeasList(:,4)==idxLambda );
    PMI{currentsub}.plotLst = lst;
    PMI{currentsub}.plot = [MeasList(lst,1) MeasList(lst,2)];
    currentsrs.type = 0;
    currentsrs.SDPairs = 0 ;
    DispParameter.reset = 1;
end

%deselection d'un canal MeasListAct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.radio_remove,'value')|get(handles.radio_restore,'value')|get(handles.radio_removeselection,'value');
    if currentsrs.type == 0 %si on a rien on prend une source ou un detecteur
        pClosest = get_SrsDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint' ),get( handles.axes_Mtg2,  'CameraPosition' ) );
        if ~(pClosest == 0)
            DetSrs = sMtg.v_HolesMtg(pClosest);
            [idxMin,is_srs] = Srs_n2SDpairs(DetSrs);
            currentsrs.SDPairs = idxMin;
            currentsrs.type = is_srs;
            plot3(vHoles(pClosest).Coord.x,vHoles(pClosest).Coord.y,vHoles(pClosest).Coord.z,'+','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',3);
            DispParameter.reset = 0;
        end
    elseif currentsrs.type == 1 %si on a un source on cherche un detecteur pres et on cree une paire
        Det_list = MeasList([find( MeasList(:,1)==currentsrs.SDPairs & MeasList(:,4)==idxLambda)],2);
        Det_listind = sMtg.v_HolesMtg(sMtg.v_pDet)==Det_list(1)*1000000
        if ~isempty(Det_list)
            pClosest = get_pDetClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint'), Det_list);
            if ~(pClosest == 0)
                 sMtg.v_HolesMtg(pClosest)
                 %sMtg.v_HolesMtg(sMtg.v_pDet)/1000000
                det = sMtg.v_HolesMtg(pClosest)/1000000;
                plotLstc = [find( MeasList(:,1)==currentsrs.SDPairs & MeasList(:,2)== det & MeasList(:,4)==idxLambda)];
                plotLstc1 = [find( MeasList(:,1)==currentsrs.SDPairs & MeasList(:,2)== det & MeasList(:,4)==1)];
                plotLstc2 = [find( MeasList(:,1)==currentsrs.SDPairs & MeasList(:,2)== det & MeasList(:,4)==2)];
                if get(handles.radio_remove,'value')
                    PMI{currentsub}.data(cf).MeasListAct(plotLstc1) = 0;
                    PMI{currentsub}.data(cf).MeasListAct(plotLstc2) = 0;
                elseif get(handles.radio_restore,'value')
                    PMI{currentsub}.data(cf).MeasListAct(plotLstc1) = 1;
                    PMI{currentsub}.data(cf).MeasListAct(plotLstc2) = 1;
                elseif get(handles.radio_removeselection,'value')
                    ind = find(PMI{currentsub}.plotLst == plotLstc);
                    %lst = PMI{currentsub}.plotLst(ind);
                    PMI{currentsub}.plotLst(ind) =[];
                    PMI{currentsub}.plot(ind,:) = [];
                    if length(PMI{currentsub}.plotLst) < 1 % If there is no channel left...
                        PMI{currentsub}.plotLst=[1]; % Add a default one to avoid an exception.
                        PMI{currentsub}.plot=[0,0];
                    end
                end
            end
        end
        currentsrs.type = 0;
        currentsrs.SDPairs = 0 ;
        DispParameter.reset = 1;
    elseif currentsrs.type == 2 %si on a un detecteur on cherche une source pres et on cree une paire
        Srs_list = MeasList([find( MeasList(:,2)==currentsrs.SDPairs & MeasList(:,4)==idxLambda)],1);
        pClosest = get_pSrsClicked( DispHelm, get( handles.axes_Mtg2,  'CurrentPoint'), Srs_list);
        if  ~(pClosest == 0)
            [srs,is_srs] = Srs_n2SDpairs(sMtg.v_HolesMtg(pClosest));
            plotLstc1 = [find( MeasList(:,1)== srs & MeasList(:,2)== currentsrs.SDPairs & MeasList(:,4)==1)];
            plotLstc =  plotLstc1
            plotLstc2 = [find( MeasList(:,1)== srs & MeasList(:,2)== currentsrs.SDPairs & MeasList(:,4)==2)];
            if get(handles.radio_remove,'value')
                PMI{currentsub}.data(cf).MeasListAct(plotLstc1) = 0;
                PMI{currentsub}.data(cf).MeasListAct(plotLstc2) = 0;
            elseif get(handles.radio_restore,'value')
                PMI{currentsub}.data(cf).MeasListAct(plotLstc1) = 1;
                PMI{currentsub}.data(cf).MeasListAct(plotLstc2) = 1;
            elseif get(handles.radio_removeselection,'value')
                ind = find(PMI{currentsub}.plotLst == plotLstc);
                %lst = PMI{currentsub}.plotLst(ind);
                PMI{currentsub}.plotLst(ind) =[];
                PMI{currentsub}.plot(ind,:) = [];
            end
        end
        currentsrs.type = 0;
        currentsrs.SDPairs = 0 ;
        DispParameter.reset = 1;
    end
end

set(guiHOMER,'UserData',PMI)
%Actualisation des fenêtres utilisés
if DispParameter.reset == 1
    if get(handles.radio_guiSPMnirsHSJ,'value')==1 %HOMER
        if get(HOMERhandles.HOMer_Tab,'value')
            plotAxes_SDG( HOMERhandles, 0 );
            plotAxes3(HOMERhandles, 0 );
            plotAxes1( HOMERhandles, 0, 0, 0 );
        end
        if get(HOMERhandles.LISA_tab,'value')
            plotAxes_LisaA(HOMERhandles,0);
            plotAxes_LisaB(HOMERhandles,0);
            plotAxes_LisaAmoinsB(HOMERhandles,0);
        end
        if get(HOMERhandles.MMARGE_Tab,'value')
            plot_MMARGE_axes_avg(HOMERhandles,0)
            plot_MMARGE_axes_avgHistogramme(HOMERhandles,0)
        end
        if ~isempty(getappdata(0,'guiGML'))& isfield(PMI{currentsub}.data(cf),'GML')
            guiGML = getappdata(0,'guiGML');
            GMLhandles = guihandles(guiGML);
            plot_axes_donnees(GMLhandles,0);
            plot_axes_corr(GMLhandles,0);
            plot_axes_yfrh(GMLhandles,0);
        end
    elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
        plotAxes_d1(HOMERhandles);
    elseif get(handles.radio_guiSPMnirsHSJ,'value')==3

    end
    resetview(handles)
end

% --- Executes on button press in radio_select_mul_pairs.
function radio_select_mul_pairs_Callback(hObject, eventdata, handles)
% hObject    handle to radio_select_mul_pairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in radio_viewfront.
function radio_viewfront_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewfront (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetview(handles)

% --- Executes on button press in radio_viewfiducies.
function radio_viewfiducies_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewfiducies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetview(handles)

% --- Executes on button press in radio_viewchannels.
function radio_viewchannels_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% resetview(handles)

% --- Executes on button press in radio_viewdistance.
function radio_viewdistance_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewdistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_viewdistance
% set(handles.radio_IDnb,'value',0)
% set(handles.radio_nbavg,'value',0)
resetview(handles)



% --- Executes on button press in btn_selectview.
function btn_selectview_Callback(hObject, eventdata, handles)
% hObject    handle to btn_selectview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Donnees de homer et de iomtg
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==4
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
HOMERhandles = guihandles(guiHOMER);
PMI = get(guiHOMER,'UserData');
global currentsub;
PMI = get(HOMERhandles.figure1,'UserData');
cf = PMI{currentsub}.currentFile;
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
DispHelm = get_Helmet(PrjStruct);
sMtg = get_Mtg( DispHelm );
vHoles = get_vHoles( DispHelm );

%Matrice d'aligement du plan par rapport à la camera
vAffichageInitial = [0,0,1]';
vAffichageChoisi = get(handles.axes_Mtg2,'CameraPosition')';
if vAffichageChoisi(1)==0 & vAffichageChoisi(2)==0
    vAffichageChoisi = [0.24202;-0.027849;2.1865];
end
matrixAlignement = AlignementVecteurs(vAffichageChoisi,vAffichageInitial );
%Matrice de transformation rotation dans le plan 2D en par rapport au vecteur up de la camera
vup = get(handles.axes_Mtg2,'CameraUpVector');
vup_r = [vup,1]*matrixAlignement;
theta = -atan2(vup_r(1),vup_r(2));
if vup(2)==0 & vup(3)==0
    theta = 1.5;
end
matrixRotation = makehgtform('zrotate',theta);
PMI{currentsub}.matTransform = matrixAlignement*matrixRotation;

%Transformation des positions SD.srs et SD.det
matSrcPos = [];
matDetPos = [];
oMRI = get_MRI_Data(PrjStruct);

if get(handles.radio_viewskin,'value')
    [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
    if ~isempty(VertexBuffer)
        dist_skin = questdlg('Do you want optode position with projection distance on skin');
    end
elseif get(handles.radio_viewcortex,'value')
    [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes(oMRI );
    if ~isempty(VertexBuffer)
        dist_skin = questdlg('Do you want optode position with projection distance on cortex');
    end
else
    dist_skin = 'No';
end
if strcmp(dist_skin,'Cancel');
    return
end
for srs = 1:numel(PMI{currentsub}.SD.SrcPos(:,1))
    Srs_n = SDpairs2Srs_n(srs,1);
    Srs_n2 = mod(Srs_n,1000)*1000 + floor(Srs_n/1000);
    p_srs = find(sMtg.v_HolesMtg == Srs_n |  sMtg.v_HolesMtg == Srs_n2);
    if ~isempty(p_srs)
        if strcmp(dist_skin,'No')
            srsx = vHoles(p_srs).Coord.x;
            srsy = vHoles(p_srs).Coord.y;
            srsz = vHoles(p_srs).Coord.z;
        elseif get(handles.radio_viewskin,'value')
            srsx = vHoles(p_srs).Coord.x-vHoles(p_srs).Normal.x*vHoles(p_srs).SkinDepth;
            srsy = vHoles(p_srs).Coord.y-vHoles(p_srs).Normal.y*vHoles(p_srs).SkinDepth;
            srsz = vHoles(p_srs).Coord.z-vHoles(p_srs).Normal.z*vHoles(p_srs).SkinDepth;
        elseif get(handles.radio_viewcortex,'value')
            srsx = vHoles(p_srs).Coord.x-vHoles(p_srs).Normal.x*vHoles(p_srs).CortexDepth;
            srsy = vHoles(p_srs).Coord.y-vHoles(p_srs).Normal.y*vHoles(p_srs).CortexDepth;
            srsz = vHoles(p_srs).Coord.z-vHoles(p_srs).Normal.z*vHoles(p_srs).CortexDepth;
        end
    else
        srsx = 0;
        srsy = 0;
        srsz = 0;
    end
    matSrcPos = [matSrcPos;srsx,srsy,srsz,1];

end
for det = 1:numel(PMI{currentsub}.SD.DetPos(:,1));
    Det_n = SDpairs2Srs_n(det,2);
    p_det = find(sMtg.v_HolesMtg == Det_n);
    if ~isempty(p_det)
        if strcmp(dist_skin,'No')
            detx = vHoles(p_det).Coord.x;
            dety = vHoles(p_det).Coord.y;
            detz = vHoles(p_det).Coord.z;
        elseif get(handles.radio_viewskin,'value')
            detx = vHoles(p_det).Coord.x-vHoles(p_det).Normal.x*vHoles(p_det).SkinDepth;
            dety = vHoles(p_det).Coord.y-vHoles(p_det).Normal.y*vHoles(p_det).SkinDepth;
            detz = vHoles(p_det).Coord.z-vHoles(p_det).Normal.z*vHoles(p_det).SkinDepth;
        elseif get(handles.radio_viewcortex,'value')
            detx = vHoles(p_det).Coord.x-vHoles(p_det).Normal.x*vHoles(p_det).CortexDepth;
            dety = vHoles(p_det).Coord.y-vHoles(p_det).Normal.y*vHoles(p_det).CortexDepth;
            detz = vHoles(p_det).Coord.z-vHoles(p_det).Normal.z*vHoles(p_det).CortexDepth;
        end
    else
        srsx = 0;
        srsy = 0;
        srsz = 0;
    end
    matDetPos = [matDetPos;detx,dety,detz,1];

end
matSrcPostr = matSrcPos * matrixAlignement*matrixRotation;
PMI{currentsub}.SD.SrcPos = matSrcPostr(:,1:3)*100;
matDetPostr = matDetPos * matrixAlignement*matrixRotation;
PMI{currentsub}.SD.DetPos = matDetPostr(:,1:3)*100;
matcoord = [PMI{currentsub}.SD.SrcPos; PMI{currentsub}.SD.DetPos];
%Selection des trous sur la surface entre z_max et z_max - profondeurz
PMI{currentsub}.data(cf).MeasList = PMI{currentsub}.data(cf).MeasListini;

%PMI{currentsub}.data(cf).MeasListAct =
if get(handles.radio_profondeurz,'value')
    set(handles.radio_profondeurz,'value',0);
    if get(handles.radio_dottedmemory,'value')
        PMI{currentsub}.data(cf).MeasListActNoisy = PMI{currentsub}.data(cf).MeasListAct;
    end
    set(handles.radio_dottedmemory,'value',0)
    profondeur_z = str2num(get(handles.edit_profondeurz,'string'));
    maxz = max(matcoord(:,3));
    minz = max(matcoord(:,3))-profondeur_z;
    for srs = 1:numel(PMI{currentsub}.SD.SrcPos(:,1))
        if find(PMI{currentsub}.SD.SrcPos(srs,3) <= maxz & PMI{currentsub}.SD.SrcPos(srs,3) >= minz);
        else
            PMI{currentsub}.SD.SrcPos(srs,:) = [0,0,0];
            Meas_delete = find(PMI{currentsub}.data(cf).MeasList(:,1)== srs);
            PMI{currentsub}.data(cf).MeasList(Meas_delete,:) = 0;
            PMI{currentsub}.data(cf).MeasListAct(Meas_delete) = 0;
            %sauvegarder l'ancien
        end
    end
    for det = 1:numel(PMI{currentsub}.SD.DetPos(:,1))
        if find(PMI{currentsub}.SD.DetPos(det,3) <= maxz & PMI{currentsub}.SD.DetPos(det,3) >= minz)
        else
            PMI{currentsub}.SD.DetPos(det,:) = [0,0,0];
            Meas_delete = find(PMI{currentsub}.data(cf).MeasList(:,2)== det);
            PMI{currentsub}.data(cf).MeasList(Meas_delete,:) = 0;
            PMI{currentsub}.data(cf).MeasListAct(Meas_delete) = 0;
        end
    end
else
    set(handles.radio_profondeurz,'value',1);
    set(handles.radio_dottedmemory,'value',1)
    if isfield(PMI{currentsub}.data(cf),'MeasListActNoisy')
        PMI{currentsub}.data(cf).MeasListAct = PMI{currentsub}.data(cf).MeasListActNoisy;
    else
        if get(handles.radio_dottedmemory,'value')
            PMI{currentsub}.data(cf).MeasListActNoisy = PMI{currentsub}.data(cf).MeasListAct;
        end
    end
end

%use to homer
%Ajustement des échelles x et y pour l'affichage de l'image projeté, on
minx = min(matcoord(:,1));
maxx = max(matcoord(:,1));
miny = min(matcoord(:,2));
maxy = max(matcoord(:,2));
%diminue légerement la précision pour que les matrices soient de la
%même  23*23
stepx =  (round((maxx-minx)/22 * 1000))/1000;
maxx =  minx + stepx*22;
stepy = (round((maxy-miny)/22 * 1000))/1000;
maxy =  miny + stepy*22;
Xvox = sprintf('%f:%f:%f',minx,stepx,maxx);
Yvox = sprintf('%f:%f:%f',miny,stepy,maxy);
set(HOMERhandles.HOMer_edit_Xvox,'string',Xvox);
set(HOMERhandles.HOMer_edit_Yvox,'string',Yvox);

% Sauvegarde des parametre PMI le MeasList, SD.SrcPos, SD.DetPos,
% PMI{currentsub}.matTransform on ete transformé, actualisation des images
% de homer
set(HOMERhandles.figure1,'UserData',PMI);
plotAxes_SDG( HOMERhandles, 0 );
plotAxes3( HOMERhandles, 0 );
plotAxes1( HOMERhandles, 0, 0, 0 );
resetview(handles)



function edit_profondeurz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_profondeurz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_profondeurz as text
%        str2double(get(hObject,'String')) returns contents of edit_profondeurz as a double


% --- Executes during object creation, after setting all properties.
function edit_profondeurz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_profondeurz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_profondeurz.
function radio_profondeurz_Callback(hObject, eventdata, handles)
% hObject    handle to radio_profondeurz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_profondeurz




% --- Executes on button press in radio_viewskin.
function radio_viewskin_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewskin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_viewskin
if get(handles.radio_viewskin,'value')
    set(handles.radio_viewcortex,'value',0)
    set(handles.radio_viewatlas,'value',0) 
    datacursormode off
end
    resetview(handles)



% --- Executes on slider movement.
function slider_time_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==4
    1
end
HOMERhandles = guihandles(guiHOMER);
PMI = get(guiHOMER,'UserData');
global currentsub
cf = PMI{currentsub}.currentFile;
scroll_ratio = get(handles.slider_time,'value');
if scroll_ratio == 0
    echantillon_time = 1;
elseif isfield(PMI{currentsub},'data') &  isfield(PMI{currentsub}.data(cf).HRF,'tHRF')
    echantillon_time = floor(scroll_ratio*size(PMI{currentsub}.data(cf).HRF.tHRF,2));
    %val = find(str2num(get(handles.text_Time,'string'))>=PMI{currentsub}.data(cf).HRF.tHRF);
    %echantillon_time =val(end);
else
    msgbox('Please average the data before');
    set(handles.slider_time,'value',0);
    guidata(hObject,handles);
    return
end

if get(handles.radio_HBO,'value'), dconc = 1;,...
elseif get(handles.radio_HBR, 'value'),dconc = 2;,...
elseif get(handles.radio_HBT,'value'),dconc = 3;,...
else
dconc = 1;
end

type=get(handles.popupmenu_display,'value');
if get(handles.radio_guiSPMnirsHSJ,'value')==1 %Homer old version
    d1=selectimagetype(echantillon_time,type,dconc);
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    %actually used 
    type=get(handles.popupmenu_display,'value');
    typelabel = get(handles.popupmenu_display,'string');
    nbhalf = size(PMI{currentsub}.data(cf).HRF.AvgC,2)/2;
    if strcmp(typelabel{type}, 'Current module') %type==1 %Current module value 
        if dconc ==1 %HbO ou 830
            d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf);
        elseif dconc==2 %HbR ou 690
            d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,nbhalf+1:end);
        elseif dconc==3 %HbT HbO+HbR ou 830+690
            d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf)+...
                PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,nbhalf+1:end);
        end
    elseif type==0 %TOPO do not support topo for release no in doc
         [name,path]= uigetfile({'*.vcolor';'*.img'});
         if path==0
             d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf);
             set(handles.popupmenu_display,'string',{'Current module',['Topo ']});
         else
             [path1, name1, ext] = fileparts([path name]);
             d1= opentopo([path name]);
             set(handles.popupmenu_display,'string',{'Current module',['Topo ',name1]});
         end
    elseif strcmp(typelabel{type}, 'Mean start stop time')  %Mean time start time stop 
    echantillon_timestart = find(str2num(get(HOMERhandles.edit_time_start,'string'))<=PMI{currentsub}.data(cf).HRF.tHRF);
    echantillon_timestop = find(str2num(get(HOMERhandles.edit_time_stop,'string'))<=PMI{currentsub}.data(cf).HRF.tHRF);
    echantillon_time = echantillon_timestart(1): echantillon_timestop(1) ;
       if dconc ==1 %HbO ou 830
            d1 = mean(PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf));
        elseif dconc==2 %HbR ou 690
            d1 = mean(PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,nbhalf+1:end));
        elseif dconc==3 %HbT HbO+HbR ou 830+690
            d1 = mean(PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf))+...
                mean(PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,nbhalf+1:end));
        end
        echantillon_time = echantillon_time(1);
    elseif strcmp(typelabel{type}, 'Projection Channel .cp')  %load D1 matrix 'Projection Channel .cp',
        [name,path]= uigetfile({'*.mat';'Select channel projection *.cp'});
        d1 = load('-mat',[path name]);
        if isfield(d1,'zonelist')
        zonelist = d1.zonelist;
        val = d1.A;
        d1 = zeros( size(PMI{1}.data.MeasList,1)/2,1);     
        for i=1:numel(zonelist)
                       
            [DetL,SrsL]=strtok(zonelist{i},' ')
            if numel(DetL)>=2 %nirx           
                SDdetL = StrBoxy2SDDet(DetL);
                SDsrsL = str2num(SrsL(3:end));
            else
                SDdetL = StrBoxy2SDDet_ISS(DetL);
                SDsrsL = StrBoxy2SDPairs(SrsL(2:end));
            end
            ML_new = PMI{1,1}.data.MeasList;
            L1 = find(ML_new(:,1)==SDsrsL & ML_new(:,2)==SDdetL & ML_new(:,4)==1);
            %L2 = find(ML_new(:,1)==SDsrsL & ML_new(:,2)==SDdetL & ML_new(:,4)==2);
            if isempty(L1)
                sprintf(['check ', DetL,' ', SrsL]);
                listname{i,1} = [DetL ' ' SrsL];
            else
                listHBOch(i,1)= L1;           
                listname{i,1} = [DetL ' ' SrsL];    
            end
        end
         d1(listHBOch) = val
        else
             d1 = d1.A;            
            
        end            
      
        set(handles.edit_D1matrix,'string',[path,name]) 
        get(handles.edit_D1matrix,'position') 
         set(handles.edit_D1matrix,'value',1)
         guidata(hObject,handles);
    elseif strcmp(typelabel{type}, 'PARAFAC')  %PARAFAC FACTOR
         FacSpatial = PMI{currentsub}.tmpPARAFAC.Factors{2}
         selected = PMI{currentsub}.tmpPARAFAC.selected;
         d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1)
         if dconc==1
            d1(PMI{1}.tmpPARAFAC.listgood,1) =FacSpatial(:,selected);
         elseif dconc==2
             d1(PMI{1}.tmpPARAFAC.listgood,1) =-FacSpatial(:,selected);
         end
         %d1 =  FacSpatial(:,selected);
    elseif strcmp(typelabel{type}, 'ICA') %ICA FACTOR
         FacSpatial = PMI{currentsub}.tmpICA.Factors{2}
         selected =PMI{currentsub}.tmpICA.selected;
         if dconc==1
         d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1)
         d1(PMI{1}.tmpICA.listgood(1:end/2),1) =FacSpatial(selected,1:end/2)';
         elseif dconc==2
            d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1)
            d1(PMI{1}.tmpICA.listgood(1:end/2),1) =FacSpatial(selected,end/2+1:end)
         end
         %d1 =  FacSpatial(:,selected);
     elseif strcmp(typelabel{type}, 'PCA') %PCA FACTOR
         s = PMI{currentsub}.tmpPCA.s;
         u = PMI{currentsub}.tmpPCA.u;
         v = PMI{currentsub}.tmpPCA.v;
         listgood = PMI{currentsub}.tmpPCA.listgood;
         lstSV =PMI{currentsub}.tmpPCA.selected; %composante 1
         temp = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
%          figure;plot(u(:,lstSV))
%           figure;plot(temp)
         ampl= s(lstSV,lstSV)*v(:,lstSV)'
         id830 = find(listgood<=size(PMI{1}.data.MeasListAct,1)/2)
         id690 = find(listgood>(size(PMI{1}.data.MeasListAct,1)/2))
         if dconc==1
             d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1);
            d1(PMI{1}.tmpPCA.listgood(id830),1) =ampl(id830);
         elseif dconc==2
            d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1)
            d1(PMI{1}.tmpPCA.listgood(id690)-(size(PMI{1}.data.MeasListAct,1)/2),1)=ampl(id690);
         end
         %d1 =  FacSpatial(:,selected);
    elseif type==0 %Detrend
        if 1
            if PMI{currentsub}.tmpDetrend.display == 1
                B = PMI{currentsub}.tmpDetrend.beta(:,1);
            elseif PMI{currentsub}.tmpDetrend.display == 2  %Show Slope B2
                 B = PMI{currentsub}.tmpDetrend.beta(:,2);
            elseif PMI{currentsub}.tmpDetrend.display == 3
                B = PMI{currentsub}.tmpDetrend.Rsquare;
            end
         listgood = PMI{currentsub}.tmpDetrend.listgood;
         id830 = find(listgood<=size(PMI{1}.data.MeasListAct,1)/2)
         id690 = find(listgood>(size(PMI{1}.data.MeasListAct,1)/2))
         d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1);
         if dconc==1

            d1(PMI{1}.tmpDetrend.listgood(id830),1) =B(id830);
            elseif dconc==2
            d1(PMI{1}.tmpDetrend.listgood(id830),1) =B(id690);
         end
        end
    elseif type==0 %Current module set tstart to zero.
        idstart = find(PMI{currentsub}.tstart>PMI{currentsub}.data(cf).HRF.tHRF);
        idstart = idstart(end);
        PMI{currentsub}.data(cf).HRF.tHRF(echantillon_time);
        if dconc ==1 %HbO ou 830
            d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf)- PMI{currentsub}.data(cf).HRF.AvgC(idstart,1:nbhalf);
        elseif dconc==2 %HbR ou 690
            d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,nbhalf+1:end)- PMI{currentsub}.data(cf).HRF.AvgC(idstart,1:nbhalf);
        elseif dconc==3 %HbT HbO+HbR ou 830+690
            d1 =( PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf)+...
                PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,nbhalf+1:end))- (PMI{currentsub}.data(cf).HRF.AvgC(idstart,1:nbhalf)+...
                PMI{currentsub}.data(cf).HRF.AvgC(idstart,nbhalf+1:end))
        end
      elseif  type==0 %PARAFAC wavelet
         FacSpatial = PMI{currentsub}.tmpPARAFACwavelet.Factors{3}
         selected = PMI{currentsub}.tmpPARAFACwavelet.selected;
         d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1)
         if dconc==1
            d1(PMI{1}.tmpPARAFACwavelet.listgood,1) =FacSpatial(:,selected);
         elseif dconc==2
             d1(PMI{1}.tmpPARAFACwavelet.listgood,1) =-FacSpatial(:,selected);
         end
    elseif strcmp(typelabel{type}, 'GLM')  %GLM
         beta = PMI{currentsub}.tmpGLM.beta;
         selected = PMI{currentsub}.tmpGLM.selected;
         d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1)
         if dconc==1
            d1(PMI{1}.tmpGLM.listgood,1) =beta(selected,:)';
         elseif dconc==2
             d1(PMI{1}.tmpGLM.listgood,1) =-beta(selected,:);
         end
    elseif strcmp(typelabel{type}, 'COMPONENT')%Selected component topo
          d1 = zeros( size(PMI{1}.data.MeasListAct,1)/2,1); 
          %tmp list good could containt both wavelenght of just one if it
          %contain both separate hbo and hbr as first and second half else
          %visualised hbo as positive and hbr as negative
          if max(PMI{currentsub}.tmplistgood) <= numel(d1)
              %on affiche juste du HBO
            listgoodHbO =PMI{currentsub}.tmplistgood;
            topoHbO =PMI{currentsub}.tmptopo;
            topoHbR = -PMI{currentsub}.tmptopo; %inverse pour HbR cas parafac only half of the channel ! 
          elseif  min(PMI{currentsub}.tmplistgood)>numel(d1)
            listgoodHbR =PMI{currentsub}.tmplistgood;
            topoHbR =PMI{currentsub}.tmptopo;
          elseif numel(PMI{currentsub}.tmplistgood)>numel(d1)
             listgoodHbR = PMI{currentsub}.tmplistgood(end/2+1:end);
             listgoodHbO = PMI{currentsub}.tmplistgood(1:end/2);
             topoHbR = PMI{currentsub}.tmptopo(end/2+1:end);
             topoHbO = PMI{currentsub}.tmptopo(1:end/2);
         else
             listgoodHbR = PMI{currentsub}.tmplistgood(1:end/2);
             listgoodHbO = PMI{currentsub}.tmplistgood(1:end/2);
             topoHbR = -PMI{currentsub}.tmptopo(end/2+1:end);
             topoHbO = PMI{currentsub}.tmptopo(listgoodHbR);
         end

         if dconc==1
             [n,m]=size(topoHbO);
             if n<m
                 topoHbO =   topoHbO';
             end
            d1(listgoodHbO,1) =sum(topoHbO,2);
         elseif dconc==2
             [n,m]=size(topoHbR);
             if n<m
                 topoHbR =   topoHbR';
             end
             d1(listgoodHbO,1)=sum(topoHbR,2);
         end
   
    end
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    d1=selectimagetype(echantillon_time,type,dconc,get(handles.radio_guiSPMnirsHSJ,'value'));
end
indna = isnan(d1);
d1(indna)=0;
%d1 = d1(1:end/2)
%d1 = d1((end/2+1):end)
d1temp = d1; %log(d1);
%figure;plot(d1temp)
%d1temp = ones(size(d1));%JTHARDCODE
if numel(d1)==size(PMI{currentsub}.data(cf).MeasList,1)/2
    for numd1 = 1:size(d1temp,1)
        d1 = d1temp(numd1,:);
        PMI{1}.coloramp = (d1-min(d1temp))/(max(d1temp)-min(d1temp));
    end
end
d1 = d1temp;
set(guiHOMER,'UserData',PMI);
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
if get(handles.radio_viewskin,'value')
    type = 0;
elseif get(handles.radio_viewcortex,'value')
    type = 1;
end
if strcmp(get(handles.context_topo_radial,'checked'),'on')
    PrjStruct = display_MRIcolor(PrjStruct,PMI,d1, type);
elseif strcmp(get(handles.context_topo_IWD,'checked'),'on')
    PrjStruct = display_MRIcolor_InverseWeightedDistance(PrjStruct,PMI,d1, type);
end
    
setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);

%handles.PrjStruct = display_MRIcolorHomer(handles)
if isfield(PMI{currentsub}.data(cf).HRF,'tHRF')
    time = PMI{currentsub}.data(cf).HRF.tHRF(1,echantillon_time);
else
    time = 1;
end
set(handles.text_Time,'string',num2str(time));
handles.d1 = d1;
guidata(hObject,handles);
resetview(handles);


% --- Executes during object creation, after setting all properties.
function slider_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function SelectHoles_Callback(hObject, eventdata, handles)
% hObject    handle to SelectHoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function resetview(handles)
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
DispHelm = get_Helmet( PrjStruct );
sMtg = get_Mtg( DispHelm );
vHoles = get_vHoles( DispHelm );
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI =get(guiHOMER,'UserData');
DispParameter.axes1 = handles.axes_Mtg2;
DispParameter.viewhelmet = get(handles.radio_holes,'value');
%     DispParameter.viewhelmet = 1;
DispParameter.reset = 1;
DispParameter.channel = get(handles.radio_Channel,'value');
DispParameter.viewfront = get(handles.radio_viewfront,'value');
if get(handles.radio_viewdistance,'value')& get(handles.popup_labelch,'value')==1
    DispParameter.D1label = 1;
    DispParameter.dist = 0;
    DispParameter.idnb = 0;
    DispParameter.nbavg =0;
    if isfield(handles,'d1')
        DispParameter.d1 = handles.d1;
    else
        DispParameter.D1label = 0;
    end
elseif get(handles.radio_viewdistance,'value')&get(handles.popup_labelch,'value')==2
    DispParameter.D1label = 0;
    DispParameter.dist = 1;
    DispParameter.idnb = 0;
    DispParameter.nbavg = 0;
elseif get(handles.radio_viewdistance,'value')&get(handles.popup_labelch,'value')==3
    DispParameter.D1label = 0;
    DispParameter.dist = 0;
    DispParameter.idnb = 1;
    DispParameter.nbavg = 0; 
elseif get(handles.radio_viewdistance,'value')&get(handles.popup_labelch,'value')==4
    DispParameter.D1label = 0;
    DispParameter.dist = 0;
    DispParameter.idnb = 0;
    DispParameter.nbavg = 1;
else
    DispParameter.D1label = 0;
    DispParameter.dist = 0;
    DispParameter.idnb = 0;
    DispParameter.nbavg = 0;
end
DispParameter.fiducies = get(handles.radio_viewfiducies,'value');
DispParameter.viewskin = get(handles.radio_viewskin,'value');
DispParameter.viewcortex = get(handles.radio_viewcortex,'value');
DispParameter.viewatlas = get(handles.radio_viewatlas,'value');
DispParameter.viewcortexintersection = get(handles.radio_normal_cortex,'value');
DispParameter.cthresh = str2num(get(handles.edit_cminthreshold,'string'));
DispParameter.cmin = str2num(get(handles.edit_climmin,'string'));
DispParameter.cmax = str2num(get(handles.edit_climmax,'string'));
DispParameter.viewscale = get(handles.btn_scale,'value');
DispParameter.MRIviewlight = get(handles.radio_light,'value');
DispParameter.MRIviewtransparent = get(handles.radio_transparent,'value');
DispParameter.viewMeasListAct = get(handles.radio_ViewMeasListAct,'value');
DispParameter.enumaltas = str2num(get(handles.edit_atlasview,'string'));
DispParameter.viewcortexandatlas  = get(handles.radio_cortexandatlas,'value');
DispParameter.HideNotCover = get(handles.radio_show_cover,'value');
DispParameter.LineWidth = str2num(get(handles.edit_channelwidth,'string'));
if get(handles.radio_d1color,'value')==1;
    DispParameter.viewd1chcolor  = 1;
elseif get(handles.radio_viewch,'value')==1;
    DispParameter.viewd1chcolor  = 2;
else
    DispParameter.viewd1chcolor = 0;
end
id=get(handles.popupmenu_backgroundcolor,'value');
str = get(handles.popupmenu_backgroundcolor,'string');
DispParameter.backgroundcolor = str{id};

IO_displayHelmet(PrjStruct,PMI,DispParameter);
guidata(handles.IO_HelmetMTG, handles);



% --------------------------------------------------------------------
function menu_Define_Colormap_Callback(hObject, eventdata, handles)

axes(handles.axes_Mtg2);
get(handles.IO_HelmetMTG);
colormapeditor;
cmin = str2num(get(handles.edit_climmin,'string'));
cmax = str2num(get(handles.edit_climmax,'string'));
caxis([cmin,cmax]);
c = colorbar;
set(c,'fontsize',14)


% --------------------------------------------------------------------
function menu_define_colorchannel_Callback(hObject, eventdata, handles)
% hObject    handle to menu_define_colorchannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI = get(guiHOMER,'UserData');
global currentsub
cf = PMI{currentsub}.currentFile;
PMI{currentsub}.color;
ml  = PMI{currentsub}.data(cf).MeasList;
scrsz = get(0,'ScreenSize');
scrsunit= get(0,'units')
fig = figure;
set(fig,'MenuBar','none','Name','Define color channel','NumberTitle','off','units',scrsunit)%,'windowstyle','modal');
set(fig,'Position',[1 1 scrsz(3)/2 scrsz(4)/2])
posfig = get(fig,'position')
nbchanel = size(ml,1)

pos_x = 0.05;
pos_y = 0.9;
width_x = 0.07;
height_y = 0.03;
% pos_x = posfig(1); %0.05;
% pos_y = posfig(2);
% width_x =  posfig(3)/25;
% height_y = posfig(4)/30;;

for ch = 1:size(ml,1)/2
    if pos_y < 0
        pos_y = 0.9;
        pos_x = pos_x + 2.2*width_x;
    end
    srs = SDPairs2strboxy(ml(ch,1));
    det = SDDet2strboxy(ml(ch,2));
    if PMI{currentsub}.data(cf).MeasListAct(ch)
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, ch,...
            numel(PMI{currentsub}.color)/3);
        line_button(ch) = uicontrol('BackgroundColor',[192/256,192/256,192/256],...
            'units','normalized',...
            'style','text',...
            'fontweight','bold',...
            'horizontalAlignment','center',...
            'string',[srs,'_',det],...
            'position', [pos_x pos_y width_x height_y],...
            'foregroundcolor',PMI{currentsub}.color(idxc,:),...
            'unit','characters',...
            'fontunit','normalized');
    else
        line_button(ch) = uicontrol('BackgroundColor',[192/256,192/256,192/256],...
            'units','normalized',...
            'style','text',...
            'fontangle','italic',...
            'horizontalAlignment','center',...
            'string',[srs,'_',det],...
            'foregroundcolor',PMI{currentsub}.color(idxc,:),...
            'position', [pos_x pos_y width_x height_y],...
            'unit','characters',...
            'fontunit','normalized');
    end
    handles.ch_button(ch) = uicontrol('BackgroundColor',PMI{currentsub}.color(idxc,:),...
        'units','normalized',...
        'style','pushbutton',...
        'Callback',{'selectcolor',handles,ch},...
        'horizontalAlignment','center',...
        'foregroundcolor','white',...
        'position', [pos_x+width_x  pos_y width_x height_y],...
        'unit','characters',...
        'fontunit','normalized');

    pos_y = pos_y - height_y;

end

guidata(handles.IO_HelmetMTG, handles);
for i = 1:numel(PMI{currentsub}.plotLst)
    PMI{currentsub}.color(PMI{currentsub}.plotLst(i),:)= [0,0,1]
end
set(guiHOMER,'userdata',PMI)



function edit_climmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_climmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_climmax as text
%        str2double(get(hObject,'String')) returns contents of edit_climmax as a double
cmin = str2num(get(handles.edit_climmin,'string'));
cmax = str2num(get(handles.edit_climmax,'string'));
caxis([cmin,cmax]);
c = colorbar;
set(c,'fontsize',14)
resetview(handles)
% --- Executes during object creation, after setting all properties.
function edit_climmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_climmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_climmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_climmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_climmin as text
%        str2double(get(hObject,'String')) returns contents of edit_climmin as a double
cmin = str2num(get(handles.edit_climmin,'string'));
cmax = str2num(get(handles.edit_climmax,'string'));
caxis([cmin,cmax]);
c = colorbar;
set(c,'fontsize',14)
resetview(handles)

% --- Executes during object creation, after setting all properties.
function edit_climmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_climmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in radio_viewcortex.
function radio_viewcortex_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewcortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_viewcortex
if get(handles.radio_viewcortex,'value')
    set(handles.radio_viewskin,'value',0)
    set(handles.radio_viewatlas,'value',0) 
    datacursormode off
end
resetview(handles)

% --- Executes on button press in radio_light.
function radio_light_Callback(hObject, eventdata, handles)
% hObject    handle to radio_light (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_light
resetview(handles)

% --- Executes on button press in radio_transparent.
function radio_transparent_Callback(hObject, eventdata, handles)
% hObject    handle to radio_transparent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_transparent

if get(handles.radio_transparent,'value')
    set(handles.IO_HelmetMTG,'renderer','openGL')
else
    set(handles.IO_HelmetMTG,'renderer','zbuffer')
end
resetview(handles)
guidata(handles.IO_HelmetMTG, handles);

% --- Executes on button press in radio_HBO.
function radio_HBO_Callback(hObject, eventdata, handles)
% hObject    handle to radio_HBO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_HBO

slider_time_Callback(hObject, eventdata, handles)


% --- Executes on button press in radio_HBR.
function radio_HBR_Callback(hObject, eventdata, handles)
% hObject    handle to radio_HBR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_HBR


slider_time_Callback(hObject, eventdata, handles)


% --- Executes on button press in radio_HBT.
function radio_HBT_Callback(hObject, eventdata, handles)
% hObject    handle to radio_HBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_HBT


slider_time_Callback(hObject, eventdata, handles)



function IO_Helmetclosereq(src,evnt)
handles = guidata(gcf); 
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    handlesHSJ = guidata(guiHOMER);
    set( handlesHSJ.radio_3dmontageupdate,'value',0)
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end


setappdata(0,'guiHelmet',[])
closereq


% --- Executes on button press in btn_HomerImageSkin.
function btn_HomerImageSkin_Callback(hObject, eventdata, handles)
% hObject    handle to btn_HomerImageSkin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
DOT = get(guiHOMER,'UserData');
global currentsub
cf =DOT{currentsub}.currentFile;
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
oMRI = get_MRI_Data(PrjStruct);

if get(handles.radio_viewskin,'value')
    [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
    vColor = get_SkinVcolor( oMRI );
    vColor = zeros(size(VertexBuffer,1),1);
elseif get(handles.radio_viewcortex,'value')
    %[VertexBuffer, IndexBuffer] = get_CortexMeshHiRes( oMRI );
    [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI );
    vColor = get_CortexLowResVcolor( oMRI );
    vColor = zeros(size(VertexBuffer,1),1);
else
    msgbox('Please select skin or cortex to display')
    return
end
VertexBuffer(:,4) = 1;

if( ~isempty( VertexBuffer) && ~isempty( IndexBuffer) )
    VertexBuffer(:,4) = 1;
    T = get_matTransform(oMRI);
    Vertex_tmp = VertexBuffer*T;
    if ~isfield(DOT{currentsub},'matTransform')
        msgbox('Sorry, you must select a plane and create an image to activate this function')
        return
    end
    matTransform = DOT{currentsub}.matTransform;
    if ~isfield(DOT{currentsub},'IMG')
        msgbox('Homer image concentration must be created')
        return
    end
    bebe=1 %FLAG GAUSSIAN BÉBÉ
    if get(handles.radio_HBO,'value') & isfield(DOT{currentsub}.IMG,'HbT')
        if bebe
            %GAUSSIAN POUR LES BEBE
            x = -2:0.15:2;
            y = -2:0.15:2;
            for i=1:numel(x)
                for j = 1:numel(y)
                    xval = x(i);
                    yval = y(j);
                    img(i,j)=2.71828183^(-xval^2/2+-yval^2/2);
                end
            end
            Image = img *0.19;
        else
            Image = DOT{currentsub}.IMG.HbO;
        end

    elseif get(handles.radio_HBR,'value') & isfield(DOT{currentsub}.IMG,'HbT')
        Image = DOT{currentsub}.IMG.HbR;
    elseif get(handles.radio_HBT,'value')& isfield(DOT{currentsub}.IMG,'HbT')
        Image = DOT{currentsub}.IMG.HbT;
    else
        msgbox('Homer image concentration must be created')
        return
    end
    %                 v_center = (Vertex_tmp(IndexBuffer(:,1),1:4) + Vertex_tmp(IndexBuffer(:,2),1:4)+ Vertex_tmp(IndexBuffer(:,3),1:4))/3;
    v_center = (Vertex_tmp(:,1:4) + Vertex_tmp(:,1:4)+ Vertex_tmp(:,1:4))/3;
    v_center_trp = v_center * matTransform;
    %cadrage en cm
    scale = 100;
    if 0% bebe
        x_min =0
        x_max = 0.05
        y_min = 0
        y_max = 0.05
        z_min = 0
        z_max = 0.05
%         mat_SrcDet = [DOT{currentsub}.SD.SrcPos(:,3);DOT{1}.SD.DetPos(:,3)];
%         ind = find(mat_SrcDet(:));
    else
        x_min =  DOT{currentsub}.IMG.Medium.CompVol.X(1)/scale%-0.005;
        x_max =  DOT{currentsub}.IMG.Medium.CompVol.X(end)/scale%+0.005;
        y_min = DOT{currentsub}.IMG.Medium.CompVol.Y(1)/scale%-0.005;
        y_max =DOT{currentsub}.IMG.Medium.CompVol.Y(end)/scale%+0.005;
        mat_SrcDet = [DOT{currentsub}.SD.SrcPos(:,3);DOT{1}.SD.DetPos(:,3)];
        ind = find(mat_SrcDet(:));
        z_min = min(mat_SrcDet(ind))/scale - 0.05;
        z_max = max(mat_SrcDet(ind))/scale + 0.05;
    end
    %   par pixel d'affichage
    h = waitbar(0,'Please wait...');
    [num_y num_x]= size(Image);
    for x = 1:(num_x)
        waitbar(x/num_x)
        x_min1 = x_min + (x-1)*(x_max-x_min)/num_x;
        x_max1 = x_min + (x)*(x_max-x_min)/num_x;
        for y = 1:(num_y)
            y_min1 = y_min + (y-1)*(y_max-y_min)/num_y;
            y_max1 = y_min +(y)*(y_max-y_min)/num_y;
            p_in = find(v_center_trp(:,1) > x_min1 & v_center_trp(:,1) < x_max1 & v_center_trp(:,2) > y_min1 & v_center_trp(:,2) < y_max1 & v_center_trp(:,3) > z_min & v_center_trp(:,3) < z_max);
            for i = 1 : numel(p_in)
                vColor(p_in(i),1) = Image(y,x);
            end
        end
    end
end
close(h)
if get(handles.radio_viewskin,'value')
    oMRI = set_SkinVcolor(oMRI,vColor);
elseif get(handles.radio_viewcortex,'value')
    oMRI = set_CortexLowResVcolor(oMRI,vColor);
end

PrjStruct = set_MRI_Data(PrjStruct, oMRI);
setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
resetview(handles);

% --- Executes on button press in btn_ResetImage.
function btn_ResetImage_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ResetImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
oMRI = get_MRI_Data(PrjStruct);
vColor = get_SkinVcolor( oMRI );
vColor = vColor*0;
oMRI = set_SkinVcolor(oMRI,vColor);
vColor = get_CortexLowResVcolor(oMRI);
vColor = vColor*0;
oMRI = set_CortexLowResVcolor(oMRI,vColor);
PrjStruct = set_MRI_Data(PrjStruct, oMRI);
setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
guidata(hObject,handles);
resetview(handles);

% --- Executes on button press in radio_holes.
function radio_holes_Callback(hObject, eventdata, handles)
% hObject    handle to radio_holes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_holes
resetview(handles);



% --- Executes on button press in btn_GML_mean.
function btn_GML_mean_Callback(hObject, eventdata, handles)
% hObject    handle to btn_GML_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI = get(guiHOMER,'UserData');
global currentsub
cf = PMI{currentsub}.currentFile;

if get(handles.radio_HBO,'value'), dconc = 1;,...
elseif get(handles.radio_HBR, 'value'),dconc = 2;,...
elseif get(handles.radio_HBT,'value'),dconc = 3;,...
else
dconc = 1
end
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
if ~isfield(PMI{currentsub}.data(cf),'GLMstimnum')
    msgbox('You must complete glm analysis before')
    return
end

num_stim = PMI{currentsub}.data(cf).GLMstimnum
d1 = zeros(size(PMI{currentsub}.data(cf).GML,1),1);
isallGML = 1;

for ch = 1:size(PMI{currentsub}.data(cf).GML,1)
    d1(ch)= PMI{currentsub}.data(cf).GML{ch,dconc}.maxcoor;
    %     if isfield(PMI{currentsub}.data(cf).GML{ch,dconc},'mh')
    %         d1(ch) = PMI{currentsub}.data(cf).GML{ch,dconc}.mh(num_stim);
    %     else
    %
    %         d1(ch) = 0;
    %         isallGML = 0;
    %     end
end
% if isallGML == 0
%     h = msgbox('The image don''t represente all channels, only calculated GLM');
% end
if get(handles.radio_viewskin,'value')
    type = 0;
elseif get(handles.radio_viewcortex,'value')
    type = 1;
end
PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type);
setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
guidata(hObject,handles);
resetview(handles);


% --- Executes on button press in radio_Channel.
function radio_Channel_Callback(hObject, eventdata, handles)
% hObject    handle to radio_Channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_Channel
resetview(handles);



% --- Executes on button press in btn_play.
function btn_play_Callback(hObject, eventdata, handles)
% hObject    handle to btn_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.stop_status, 'value', 0);

if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI = get(guiHOMER,'UserData');
HOMERhandles = guihandles(guiHOMER);
global currentsub
cf = PMI{currentsub}.currentFile;
if ~isfield(PMI{currentsub}.data(cf).HRF,'tHRF')
    msgbox('Please average data before');
    return
end
fin =  numel(PMI{currentsub}.data(cf).HRF.tHRF);
time_step = str2num(get(handles.edit_step_movie,'string'));
debut = str2num(get(handles.text_Time,'string'));
if isempty(debut)
    set(handles.text_Time,'string', 1)
    debut = 1;
end
t = PMI{currentsub}.data(cf).HRF.tHRF;
inddebut = find(debut>t);
if ~isempty(inddebut) % This conditionnal is to make sure that at least one value satisfy the condition...
    inddebut = inddebut(end); % If there is at least one value, we can start from the latest one (closest from the chosen start time)...
else
    inddebut = 1; % Otherwise, start at the first value of the time vector.
end
step = round(time_step/(PMI{currentsub}.data(cf).HRF.tHRF(2)-PMI{currentsub}.data(cf).HRF.tHRF(1)));
for echantillon_time = inddebut:step:fin
    if get(handles.stop_status, 'value') % Verify if the stop button has been clicked.
        return
    end
    if get(handles.radio_HBO,'value'), dconc = 1;,...
    elseif get(handles.radio_HBR, 'value'),dconc = 2;,...
    elseif get(handles.radio_HBT,'value'),dconc = 3;,...
    end

if get(handles.radio_guiSPMnirsHSJ,'value')==1
        d1=selectimagetype(echantillon_time,get(handles.popupmenu_display,'value'),dconc,get(handles.radio_guiSPMnirsHSJ,'value'));
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    nbhalf = size(PMI{currentsub}.data(cf).HRF.AvgC,2)/2;
    if dconc ==1 %HbO ou 830
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf);
    elseif dconc==2 %HbR ou 690
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,nbhalf+1:end);
    elseif dconc==3 %HbT HbO+HbR ou 830+690
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,1:nbhalf)+...
            PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,nbhalf+1:end);
    end
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    d1=selectimagetype(echantillon_time,get(handles.popupmenu_display,'value'),dconc,get(handles.radio_guiSPMnirsHSJ,'value'));
end


PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
if get(handles.radio_viewskin,'value')
    type = 0;
elseif get(handles.radio_viewcortex,'value')
    type = 1;
else
    msgbox('Use skin or cortex surface to display video.')
    return
end

if type ==0 |type == 1
    PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type);
    setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
    %handles.PrjStruct = display_MRIcolorHomer(handles)
    time = PMI{currentsub}.data(cf).HRF.tHRF(1,echantillon_time);
    set(handles.text_Time,'string',num2str(time));
    guidata(hObject,handles);
    resetview(handles);
    uiwait(gcf,1);
end
end

% --- Executes on button press in btn_autoscale.
function btn_autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to btn_autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
DispHelm = get_Helmet(PrjStruct);
oMRI = get_MRI_Data(PrjStruct);
if get(handles.radio_viewskin,'value')
    vColor = get_SkinVcolor(oMRI);
elseif get(handles.radio_viewcortex,'value')
    vColor = get_CortexLowResVcolor(oMRI);
elseif get(handles.radio_viewatlas,'value')
    vColor = get_CortexHiResVcolor(oMRI);
end
max_color = max(abs(vColor));
min_color = min(vColor);
if get(handles.radio_viewatlas,'value')
    caxis([min_color, max_color]);
    min_color = sprintf('%0.2g',min_color);
    max_color = sprintf('%0.2g',max_color);
    set(handles.edit_climmin,'string',[min_color]);
    set(handles.edit_climmax,'string',[max_color]);
else
    caxis([-max_color, max_color]);
    max_color = sprintf('%0.2g',max_color);
    set(handles.edit_climmin,'string',['-', max_color]);
    set(handles.edit_climmax,'string',[max_color]);
end
cthresh = str2num(max_color) /10;
set(handles.edit_cminthreshold,'string',num2str(cthresh));
c = colorbar;
set(c,'fontsize',14)



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_step_movie_Callback(hObject, eventdata, handles)
% hObject    handle to edit_step_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_step_movie as text
%        str2double(get(hObject,'String')) returns contents of edit_step_movie as a double


% --- Executes during object creation, after setting all properties.
function edit_step_movie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_step_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function text_Time_Callback(hObject, eventdata, handles)
% hObject    handle to text_Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_Time as text
%        str2double(get(hObject,'String')) returns contents of text_Time as a double
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI = get(guiHOMER,'UserData');
global currentsub
cf = PMI{currentsub}.currentFile;



echantillon_time = find(str2num(get(handles.text_Time,'string'))<=PMI{currentsub}.data(cf).HRF.tHRF);
if isempty(echantillon_time)
    msgbox('No time exist');
    return;
end
scroll_ratio = 1-numel(echantillon_time)/numel(PMI{currentsub}.data(cf).HRF.tHRF);
set(handles.slider_time,'value',scroll_ratio);

slider_time_Callback(hObject, eventdata, handles);


% --- Executes on button press in btn_clearchannel.
function btn_clearchannel_Callback(hObject, eventdata, handles)
% hObject    handle to btn_clearchannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI = get(guiHOMER,'UserData');
global currentsub
cf =PMI{currentsub}.currentFile;
PMI{currentsub}.plotLst=[1];
PMI{currentsub}.plot=[0,0];
set(guiHOMER,'UserData',PMI)
HOMERhandles = guidata(guiHOMER);
plotAxes_d1(HOMERhandles)
resetview(handles)


% --- Executes on button press in radio_ViewMeasListAct.
function radio_ViewMeasListAct_Callback(hObject, eventdata, handles)
% hObject    handle to radio_ViewMeasListAct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_ViewMeasListAct

resetview(handles)


% --- Executes on button press in btn_exportHelmetFigure.
function btn_exportHelmetFigure_Callback(hObject, eventdata, handles,titlefigure)
% hObject    handle to btn_exportHelmetFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiEFO = getappdata(0,'guiExportFigureOption');
if ~isempty(guiEFO)
    EFOhandles = guihandles(guiEFO);
    iscolorbar = get(EFOhandles.radio_colorbar,'value');
    istitle = get(EFOhandles.radio_title,'value');
    pixelx = str2num(get(EFOhandles.edit_xpixel,'string'));
    pixely = str2num(get(EFOhandles.edit_ypixel,'string'));;
    xmax = str2num(get(EFOhandles.edit_xzoom,'string'));
    ymax = str2num(get(EFOhandles.edit_yzoom,'string'));
    zmax = str2num(get(EFOhandles.edit_zzoom,'string'));
else
    iscolorbar = 1;
    istitle = 0;
    pixelx = 1000;
    pixely = 1000;
    xmax = 0.13;
    ymax = 0.13;
    zmax = 0.13;
end
h = figure;
if isfield(handles,'exportview')
    exportview = handles.exportview;
else
    exportview.taille = 0;
    exportview.titlefigure = 1;
    exportview.autotitle = 1;
    exportview.setcolor = [1 1 1]
end
Colormapfig = get(handles.IO_HelmetMTG,'colormap');
set(h,'colormap',Colormapfig);
set(h, 'renderer','zbuffer');
size_screen=get(0,'ScreenSize');
set(h,'unit','pixel','position',[30,30,pixelx,pixely])
newaxesh = gca
prop = get(handles.axes_Mtg2,'CameraTarget');
set(newaxesh,'CameraTarget',prop);
clear prop
prop = get(handles.axes_Mtg2,'CameraPosition');
set(newaxesh,'CameraPosition',prop);
clear prop
prop = get(handles.axes_Mtg2,'CameraUpVector');
set(newaxesh,'CameraUpVector',prop);
clear prop
prop = get(handles.axes_Mtg2,'CameraViewAngle');
set(newaxesh,'CameraViewAngle',prop);
clear prop
prop = get(handles.axes_Mtg2,'CLim');
set(newaxesh,'CLim',prop);
clear prop
prop = get(handles.axes_Mtg2,'CLimMode');
set(newaxesh,'CLimMode',prop);
clear prop
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
DispHelm = get_Helmet( PrjStruct );
sMtg = get_Mtg( DispHelm );
vHoles = get_vHoles( DispHelm );
oMRI = get_MRI_Data(PrjStruct);
vColor = get_SkinVcolor(oMRI);
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI =get(guiHOMER,'UserData');
global currentsub
cf = PMI{currentsub}.currentFile;
DispParameter.axes1 = newaxesh;
DispParameter.viewhelmet = get(handles.radio_holes,'value');
DispParameter.reset = 1;
DispParameter.channel = get(handles.radio_Channel,'value');
DispParameter.viewfront = 1 %get(handles.radio_viewfront,'value');
if get(handles.radio_viewdistance,'value')& get(handles.popup_labelch,'value')==1
    DispParameter.D1label = 1;
    DispParameter.dist = 0;
    DispParameter.idnb = 0;
    DispParameter.nbavg =0;
    DispParameter.d1 = handles.d1;
elseif get(handles.radio_viewdistance,'value')&get(handles.popup_labelch,'value')==2
    DispParameter.D1label = 0;
    DispParameter.dist = 1;
    DispParameter.idnb = 0;
    DispParameter.nbavg = 0;
elseif get(handles.radio_viewdistance,'value')&get(handles.popup_labelch,'value')==3
    DispParameter.D1label = 0;
    DispParameter.dist = 0;
    DispParameter.idnb = 1;
    DispParameter.nbavg = 0; 
elseif get(handles.radio_viewdistance,'value')&get(handles.popup_labelch,'value')==4
    DispParameter.D1label = 0;
    DispParameter.dist = 0;
    DispParameter.idnb = 0;
    DispParameter.nbavg = 1;
else
    DispParameter.D1label = 0;
    DispParameter.dist = 0;
    DispParameter.idnb = 0;
    DispParameter.nbavg = 0;
end
DispParameter.fiducies = get(handles.radio_viewfiducies,'value');
DispParameter.viewskin = get(handles.radio_viewskin,'value');
DispParameter.viewcortex = get(handles.radio_viewcortex,'value');
DispParameter.viewatlas = get(handles.radio_viewatlas,'value');
DispParameter.cthresh = str2num(get(handles.edit_cminthreshold,'string'));
DispParameter.cmin = str2num(get(handles.edit_climmin,'string'));
DispParameter.cmax = str2num(get(handles.edit_climmax,'string'));
DispParameter.enumaltas = str2num(get(handles.edit_atlasview,'string'));
DispParameter.viewscale = get(handles.btn_scale,'value');
DispParameter.viewcortexintersection = get(handles.radio_normal_cortex,'value');
DispParameter.MRIviewlight = get(handles.radio_light,'value');
DispParameter.MRIviewtransparent = get(handles.radio_transparent,'value');
DispParameter.viewMeasListAct = get(handles.radio_ViewMeasListAct,'value');
DispParameter.viewcortexandatlas   = get(handles.radio_cortexandatlas,'value');
DispParameter.HideNotCover = get(handles.radio_show_cover,'value');
DispParameter.LineWidth = str2num(get(handles.edit_channelwidth,'string'));
if get(handles.radio_d1color,'value')==1;
    DispParameter.viewd1chcolor  = 1;
elseif get(handles.radio_viewch,'value')==1;
    DispParameter.viewd1chcolor  = 2;
else
    DispParameter.viewd1chcolor = 0;
end
id=get(handles.popupmenu_backgroundcolor,'value');
str = get(handles.popupmenu_backgroundcolor,'string');
DispParameter.backgroundcolor = str{id};
%     xmax = max(abs(PMI{1}.SD.SrcPos(:,1)))*0.01+0.02;
%     ymax = max(abs(PMI{1}.SD.SrcPos(:,2)))*0.01+0.02;
%     zmax = max(abs(PMI{1}.SD.SrcPos(:,3)))*0.01+0.02;

set(newaxesh,'unit','normalize')
%     set(newaxesh,'Position',[0,0,0.8,0.8]);
%     set(newaxesh,'OuterPosition',[0,0,0.8,0.8])
set(newaxesh,'Position',[0,0,1,1]);
set(newaxesh,'OuterPosition',[0,0,1,1])
set(newaxesh,'Xlim',[-xmax,xmax]);
set(newaxesh,'Ylim',[-ymax,ymax]);
set(newaxesh,'Zlim',[-zmax/3, zmax]);
set(newaxesh,'xtick',[])
set(newaxesh,'ytick',[])
set(newaxesh,'ztick',[])
set(h,'color',exportview.setcolor)
set(newaxesh,'Color',exportview.setcolor)

IO_displayHelmet(PrjStruct,PMI,DispParameter);
guidata(handles.IO_HelmetMTG, handles);

if ~isfield(PMI{currentsub},'subjectNum')
        exportview.autotitle = 0;
end

if  exportview.autotitle == 1
    if get(handles.radio_HBO,'value')
        Conc = 'HbO';
    elseif get(handles.radio_HBR,'value')
        Conc = 'HbR';
    elseif get(handles.radio_HBT,'value')
        Conc = 'HbT';
    end

    suj = PMI{currentsub}.subjectNum;

    if iscell(suj)
        a = suj{1}
        suj = a;
    end
    cond = PMI{currentsub}.data(cf).filenm;
    type= get(handles.popupmenu_display,'string');
    id = get(handles.popupmenu_display,'value');

    titlefigure = [suj, cond, 'Time : ', sprintf('%2.1f',str2num(get(handles.text_Time,'string'))),' ',type{id}, ' : ', Conc ] ;

    if istitle
        title(titlefigure);
    end
    % saveas(h,[titlefigure,'.fig'],'fig')
    set(h, 'name', [titlefigure])
    if iscolorbar
        colorbar
    end
    angle =  get(handles.angle_rotation_cote,'string')
    anglef =  get(handles.angle_rotation_face,'string')

    [name,path] = uiputfile([angle,anglef , suj, cond,sprintf('%2.1f',str2num(get(handles.text_Time,'string'))),type{id}, Conc,'.tif'])
    saveas(h,[path,name(1:end-4),'.tif'],'tif');
    saveas(h,[path,name(1:end-4),'.fig'],'fig');
    save([path,name(1:end-4),'.vcolor'], 'vColor','-mat');
else
   [name,path]= uiputfile('.img');
   if path==0
       return
   end
    savetopo([path,name],vColor,1)
    saveas(h,[path,name(1:end-4),'.tif'],'tif');
    saveas(h,[path,name(1:end-4),'.fig'],'fig');
    close(h);
end
% set(newaxesh,'xcolor',exportview.setcolor)
% set(newaxesh,'ycolor',exportview.setcolor)
% set(newaxesh,'zcolor',exportview.setcolor)




% --- Executes on button press in btn_openimg.
function btn_openimg_Callback(hObject, eventdata, handles)
% hObject    handle to btn_openimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('.img');
if filename == 0
    return
end
load([pathname,filename],'-mat');

PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
oMRI = get_MRI_Data(PrjStruct);

if get(handles.radio_viewskin,'value')
    [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
    vColor = get_SkinVcolor( oMRI );
    vColor = zeros(size(IndexBuffer,1),1);
elseif get(handles.radio_viewcortex,'value')
    %[VertexBuffer, IndexBuffer] = get_CortexMeshHiRes( oMRI );
    [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI );
    vColor = get_CortexLowResVcolor( oMRI );
    vColor = zeros(size(IndexBuffer,1),1);
else
    msgbox('Please select skin or cortex to display')
    return
end
VertexBuffer(:,4) = 1;

if( ~isempty( VertexBuffer) && ~isempty( IndexBuffer) )
    VertexBuffer(:,4) = 1;
    T = get_matTransform(oMRI);
    Vertex_tmp = VertexBuffer*T;
    if ~isfield(IMG,'matTransform')
        msgbox('Sorry, you must select a plane and create an image to activate this function')
        return
    end
    matTransform = IMG.matTransform;

    if get(handles.radio_HBO,'value') & isfield(IMG,'HbO')
        Image = IMG.HbO;
    elseif get(handles.radio_HBR,'value') & isfield(IMG,'HbR')
        Image = IMG.HbR;
    elseif get(handles.radio_HBT,'value')& isfield(IMG,'HbT')
        Image = IMG.HbT;
    else
        msgbox('Homer image concentration must be created')
        return
    end
    v_center = (Vertex_tmp(IndexBuffer(:,1),1:4) + Vertex_tmp(IndexBuffer(:,2),1:4)+ Vertex_tmp(IndexBuffer(:,3),1:4))/3;
    v_center_trp = v_center * matTransform;
    x_min =  IMG.xaxis(1)/100;
    x_max =  IMG.xaxis(end)/100;
    y_min =  IMG.yaxis(1)/100;
    y_max =  IMG.yaxis(end)/100;
    mat_SrcDet = [ IMG.SD.SrcPos(:,3);IMG.SD.DetPos(:,3)];
    ind = find(mat_SrcDet(:));
    z_min = min(mat_SrcDet(ind))/100 - 0.01;
    z_max = max(mat_SrcDet(ind))/100 + 0.01;
    %   par pixel d'affichage
    [num_y num_x]= size(Image);
    for x = 1:(num_x)
        x_min1 = x_min + (x-1)*(x_max-x_min)/num_x;
        x_max1 = x_min + (x)*(x_max-x_min)/num_x;
        for y = 1:(num_y)
            y_min1 = y_min + (y-1)*(y_max-y_min)/num_y;
            y_max1 = y_min +(y)*(y_max-y_min)/num_y;
            p_in = find(v_center_trp(:,1) > x_min1 & v_center_trp(:,1) < x_max1 & v_center_trp(:,2) > y_min1 & v_center_trp(:,2) < y_max1 & v_center_trp(:,3) > z_min & v_center_trp(:,3) < z_max);
            for i = 1 : numel(p_in)
                vColor(p_in(i),1) = Image(y,x);
            end

        end
    end
end

if get(handles.radio_viewskin,'value')
    oMRI = set_SkinVcolor(oMRI,vColor);
elseif get(handles.radio_viewcortex,'value')
    oMRI = set_CortexLowResVcolor(oMRI,vColor);
end

PrjStruct = set_MRI_Data(PrjStruct, oMRI);
setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
resetview(handles);


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in view_back.
function view_back_Callback(hObject, eventdata, handles)
% hObject    handle to view_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Mtg2)
set(handles.axes_Mtg2,'CameraTarget',[0,0,0]);
set(handles.axes_Mtg2,'CameraPosition',[-2.2,0,0]);
set(handles.axes_Mtg2,'CameraUpVector',[0,0,1]);
set(handles.angle_rotation_cote,'string','180');
set(handles.angle_rotation_face,'string','0');

% --- Executes on button press in view_top.
function view_top_Callback(hObject, eventdata, handles)
% hObject    handle to view_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Mtg2)
set(handles.axes_Mtg2,'CameraTarget',[0,0,0]);
set(handles.axes_Mtg2,'CameraPosition',[0,0,2.2]);
set(handles.axes_Mtg2,'CameraUpVector',[-1,0,0]);
set(handles.angle_rotation_cote,'string','0');
set(handles.angle_rotation_face,'string','89');

% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in view_front.
function view_front_Callback(hObject, eventdata, handles)
% hObject    handle to view_front (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Mtg2)
% center = get_Center( ans )
set(handles.axes_Mtg2,'CameraTarget',[0,0,0]);
set(handles.axes_Mtg2,'CameraPosition',[2.2,0,0]);
set(handles.axes_Mtg2,'CameraUpVector',[0,0,1]);
set(handles.angle_rotation_cote,'string','0');
set(handles.angle_rotation_face,'string','0');

% --- Executes on button press in view_right.
function view_right_Callback(hObject, eventdata, handles)
% hObject    handle to view_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_Mtg2)
set(handles.axes_Mtg2,'CameraTarget',[0,0,0]);
set(handles.axes_Mtg2,'CameraPosition',[0,-2.2,0]);
set(handles.axes_Mtg2,'CameraUpVector',[0,0,1]);
set(handles.angle_rotation_cote,'string','-90');
set(handles.angle_rotation_face,'string','0');


% --- Executes on button press in view_left.
function view_left_Callback(hObject, eventdata, handles)
% hObject    handle to view_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_Mtg2)
set(handles.axes_Mtg2,'CameraTarget',[0,0,0]);
set(handles.axes_Mtg2,'CameraPosition',[0,2.2,0]);
set(handles.axes_Mtg2,'CameraUpVector',[0,0,1]);
set(handles.angle_rotation_cote,'string','90');
set(handles.angle_rotation_face,'string','0');

% --- Executes on button press in radio_viewatlas.
function radio_viewatlas_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewatlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_viewatlas
if get(handles.radio_viewatlas,'value')
    set(handles.radio_viewcortex,'value',0)
    set(handles.radio_viewskin,'value',0) 
    dcm_obj = datacursormode(handles.IO_HelmetMTG);
    set(dcm_obj,'UpdateFcn',@myupdatefcn)
    set(dcm_obj,'SnapToDataVertex','on')
    set(dcm_obj,'DisplayStyle','window')
    datacursormode on
end
resetview(handles)


function txt = myupdatefcn(empt,event_obj)
 
pos = get(event_obj,'Position');
axeh = get(empt,'Parent');
figh = get(axeh,'Parent');
handles = guihandles(figh);
cd = get(event_obj.Target,'Vertices');
ind_cursor = find(cd(:,1)==pos(1)& cd(:,2)==pos(2) & cd(:,3)==pos(3));
cd_atlas = get(event_obj.Target,'FaceVertexCData');
fileTOPOmat =[];%get(handles.edit_TOPOmatfile,'string');
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
DispHelm = get_Helmet(PrjStruct);
oMRI = get_MRI_Data(PrjStruct);

 PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct'); 
% oMRI = get_MRI_Data(PrjStruct); 
try
if ~isempty(fileTOPOmat)
    oMRI = get_MRI_Data(PrjStruct);
    [pathtopomat,file,ext]=fileparts(fileTOPOmat);
    if isempty(get(handles.edit_hupdate,'string'))
        h=figure;
        hold on
        set(handles.edit_hupdate,'string',num2str(h));
    else
        h = str2num(get(handles.edit_hupdate,'string'));
        if ishandle(h)
            figure(h);
        else
            h=figure;
            set(handles.edit_hupdate,'string',num2str(h));
        end
        hold on
        cla
    end
    if get(handles.radio_viewskin,'value')
        [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
    else
        [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI );
    end
    set(gca,'fontsize',14)
    ylabel('Amplitude')
    xlabel('TOPO image')
    alldot=[];
    allname = [];
    TOPO= load(fileTOPOmat);
    TOPO = TOPO.TOPO;
    module = numel(TOPO.pp);
    colorbar = ['r', 'g', 'b', 'c', 'm'];
   icolor = 1
   datay = [];
    for iTOPO = 1:numel(TOPO.pp(1).p)
        fileTOPOvcolor = TOPO.pp(1).p{iTOPO};
        [pathTOPO,fileTOPO,ext]=fileparts(fileTOPOvcolor);
        vColor = opentopo(fileTOPOvcolor);
        alldot = [alldot,vColor(ind_cursor)];
        allname = [allname,{fileTOPO}];
        if ind_cursor > numel(vColor)
            h = bar(0,0);
           fprintf('WRONG FORMAT')
        else
            h = bar(iTOPO,vColor(ind_cursor),colorbar(icolor));
            datay =[datay, vColor(ind_cursor)];
        end
        icolor= icolor+1;
         if icolor == 6
             icolor = 1;
         end
        set(h,'displayname',fileTOPO,'linewidth',3);

    end

   sprintf('%d\t', datay)

    title(['Coordinate',' x= ',round(num2str(VertexBuffer(ind_cursor,1))),' y= ',...
        round(num2str(VertexBuffer(ind_cursor,2))),...
        ' z= ',round(num2str(VertexBuffer(ind_cursor,3))), ' pos=', num2str(ind_cursor)],'fontsize',14)
    set(h,'name',pathtopomat)
end



catch

end

if cd_atlas == 0
    return
end
nb_atlas= max(cd_atlas);
if get(handles.radio_viewatlas,'value')
caseatlas = get(handles.popup_atlastype,'value');
c = cd_atlas(ind_cursor);
if caseatlas == 4 %69 region
    txt = {['Atlas number : ',num2str(c),' ', findAtlas69(c)]};
elseif caseatlas == 1
    txt = {['Atlas number : ',num2str(c),' ', findAtlas116(c)]};
elseif caseatlas == 3
    txt={['Atlas number : ',num2str(c),' ', findAtlas90(c)]};
elseif caseatlas == 2
    txt = {['Atlas number : ',num2str(c),' ', findAtlasBrodmann48(c)]};
end
else
    txt = ['Amplitude :',num2str(cd_atlas(ind_cursor))];
end
warning off
datacursormode off
%msgbox(txt)


% --- Executes on button press in radio_normal_cortex.
function radio_normal_cortex_Callback(hObject, eventdata, handles)
% hObject    handle to radio_normal_cortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_normal_cortex
resetview(handles)



% --- Executes on button press in radio_dottedmemory.
function radio_dottedmemory_Callback(hObject, eventdata, handles)
% hObject    handle to radio_dottedmemory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_dottedmemory





function edit_cminthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cminthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cminthreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_cminthreshold as a double
resetview(handles)
% cmap = jet(50);
% cmin = str2num(get(handles.edit_climmin,'string'));
% cmax = str2num(get(handles.edit_climmax,'string'));
% nb_color = size(cmap,1);
% delta_map = (cmax-cmin)/nb_color;
% cthresh = str2num(get(handles.edit_cminthreshold,'string'));
% nb_map_min = floor(abs(cmin+cthresh)/delta_map);
% nb_map_max = ceil(abs(cmin-cthresh)/delta_map);
% nb_map_list = nb_map_min:nb_map_max;
% nb_map = size(nb_map_list,2);
% for i_map = 1:nb_map;
%     if get(handles.radio_viewskin,'value')
%         cmap( nb_map_list(i_map),:) = [255/256 255/256 179/256];
%     elseif  get(handles.radio_viewcortex,'value')
%         cmap( nb_map_list(i_map),:) = [88/256 88/256 88/256];
%     end
% end
% colormap(cmap);
% colorbar;
% --- Executes during object creation, after setting all properties.
function edit_cminthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cminthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

return

% --- Executes on button press in btn_stop.
function btn_stop_Callback(hObject, eventdata, handles)
% hObject    handle to btn_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.stop_status, 'value', 1);
    % This is a stop gap solution as a hidden radio button is used to store
    % the play state of the UI. It is set as zero when the play button is
    % clicked and set to 1 when the stop button is pressed. The playback
    % loop will play as long as the value is zero.
return

% --- Executes on selection change in popupmenu_display.
function popupmenu_display_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_display contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_display
slider_time_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit_atlasview_Callback(hObject, eventdata, handles)
% hObject    handle to edit_atlasview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_atlasview as text
%        str2double(get(hObject,'String')) returns contents of edit_atlasview as a double
resetview(handles)



% --- Executes during object creation, after setting all properties.
function radio_viewatlas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radio_viewatlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called





function angle_rotation_Cote_Callback(hObject, eventdata, handles)
% hObject    handle to angle_rotation_Cote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_rotation_Cote as text
%        str2double(get(hObject,'String')) returns contents of angle_rotation_Cote as a double


% --- Executes during object creation, after setting all properties.
function angle_rotation_Cote_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_rotation_Cote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle_rotation_cote_Callback(hObject, eventdata, handles)
% hObject    handle to angle_rotation_cote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_rotation_cote as text
%        str2double(get(hObject,'String')) returns contents of angle_rotation_cote as a double
angle_cote = str2num(get(handles.angle_rotation_cote,'string'));
angle_face = str2num(get(handles.angle_rotation_face,'string'));
axes(handles.axes_Mtg2)
cp_x = cos(angle_cote*pi/180) * cos(angle_face*pi/180)*2.2;
cp_y = sin(angle_cote*pi/180)* cos(angle_face*pi/180)*2.2;
cp_z = sin(angle_face*pi/180)*2.2;
set(handles.axes_Mtg2,'CameraTarget',[0,0,0]);
set(handles.axes_Mtg2,'CameraPosition',[cp_x,cp_y,cp_z]);
set(handles.axes_Mtg2,'CameraUpVector',[0,0,1]);


% --- Executes during object creation, after setting all properties.
function angle_rotation_cote_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_rotation_cote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle_rotation_face_Callback(hObject, eventdata, handles)
% hObject    handle to angle_rotation_face (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_rotation_face as text
%        str2double(get(hObject,'String')) returns contents of angle_rotation_face as a double




% --- Executes during object creation, after setting all properties.
function axes6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes6


% --- Executes during object creation, after setting all properties.
function axes_angle_rotation_cote_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_angle_rotation_cote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_angle_rotation_cote
axes(hObject)




% --- Executes during object creation, after setting all properties.
function text_angle_rotation_cote_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_angle_rotation_cote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function edit_atlasview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_atlasview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function angle_rotation_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_rotation_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popup_atlastype.
function popup_atlastype_Callback(hObject, eventdata, handles)
% hObject    handle to popup_atlastype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_atlastype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_atlastype


% --- Executes during object creation, after setting all properties.
function popup_atlastype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_atlastype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_cortexandatlas.
function radio_cortexandatlas_Callback(hObject, eventdata, handles)
% hObject    handle to radio_cortexandatlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_cortexandatlas
resetview(handles)





% --- Executes on button press in btn_viewmany.
function btn_viewmany_Callback(hObject, eventdata, handles)
% hObject    handle to btn_viewmany (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)[name,path]= uigetfile('.mat');


if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI = get(guiHOMER,'UserData');

guiHelmet = getappdata(0,'guiHelmet');
PrjStruct = getappdata(guiHelmet,'PrjStruct');
type = 0;
[filename,pathname] = uigetfile('.mat','MultiSelect','on')
% if filename == 0
%     return
% end
if ~iscell(filename)
    filename = {filename}
end
[path] = uigetdir()
for i = 1:numel(filename)
    d1 = load('-mat',[pathname filename{i}]);
    d1temp= d1.A;
    for id1temp = 1:size(d1temp,1)
        d1 = d1temp(id1temp,:);
        PrjStruct = display_MRIcolor(PrjStruct,PMI,d1, type);
        setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
        resetview(handles);
        btn_exportHelmetFigure_Callback(hObject, eventdata, handles,[path,'\',num2str(id1temp),filename{i}]);
    end
end






% --- Executes on button press in radio_d1color.
function radio_d1color_Callback(hObject, eventdata, handles)
% hObject    handle to radio_d1color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_d1color
set(handles.radio_viewch, 'value', 0);
resetview(handles);


% --------------------------------------------------------------------
function menu_definecolordet_Callback(hObject, eventdata, handles)
% hObject    handle to menu_definecolordet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI = get(guiHOMER,'UserData');
global currentsub;
cf = PMI{currentsub}.currentFile;
ml = PMI{currentsub}.data(cf).MeasList;
nb_ch = size(ml,1);
half_ch = nb_ch/2;

colordef = [   0.7500    0.7500    0.7500
    1.0000         0         0
    0    1.0000         0
    0    1.0000    1.0000
    0         0    1.0000
    1.0000         0    1.0000
    0         0         0
    0.5000    0.5000    0.5000
    0.5000         0         0
    0.5000    0.5000         0
    0    0.5000         0
    0         0    0.5000
    0.5000         0    0.5000
    30/256          100/256     100/256
    0               100/256     25/256
    179/256             50/256      50/256 ]; %orange
% figure
% hold on
% for i = 1:16
%     plot(1,i,'*','color',colordef(i,:))
% end
button = questdlg('Do you want to attribute same color for each detector','Color Display','Yes','No','Yes');
if strcmp(button,'Yes')
    for idet = 1:16
        ind = find(ml(1:half_ch,2)==idet);
        for i = 1:numel(ind)
            PMI{currentsub}.color(ind(i),:)= colordef(idet,:);
        end
    end
else
    i.color =colorcube(half_ch);
    for ind = 1:size(i.color,1)
        sum(i.color(ind,:));
        if sum(i.color(ind,:))>2.3
            i.color(ind,:) = [0,0,0];
        end
        if i.color(ind,1)==255 & i.color(ind,1)
            i.color(ind,:) = [0,0,0];
        end
    end
    PMI{currentsub}.color = i.color;
end
set(guiHOMER,'UserData',PMI);
updataallview();


% --- Executes on button press in btn_viewzone.
function btn_viewzone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_viewzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentsub;
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
DispHelm = get_Helmet( PrjStruct );
sMtg = get_Mtg( DispHelm );
vHoles = get_vHoles( DispHelm );
matSrcPos = [];
matDetPos = [];
for srs = 1:numel(DOT{currentsub}.SD.SrcPos(:,1))
    Srs_n = SDpairs2Srs_n(srs,1);
    Srs_n2 = mod(Srs_n,1000)*1000 + floor(Srs_n/1000);
    p_srs = find(sMtg.v_HolesMtg == Srs_n |  sMtg.v_HolesMtg == Srs_n2);
    if ~isempty(p_srs)
        srsx = vHoles(p_srs).Coord.x*100;
        srsy = vHoles(p_srs).Coord.y*100;
        srsz = vHoles(p_srs).Coord.z*100;
        matSrcPos = [matSrcPos;srsx,srsy,srsz];
    else
        matSrcPos = [matSrcPos;0,0,0];
    end
end
for det = 1:numel(DOT{currentsub}.SD.DetPos(:,1));
    Det_n = SDpairs2Srs_n(det,2);
    p_det = find(sMtg.v_HolesMtg == Det_n);
    if ~isempty(p_det)
        detx = vHoles(p_det).Coord.x*100;
        dety = vHoles(p_det).Coord.y*100;
        detz = vHoles(p_det).Coord.z*100;
        matDetPos = [matDetPos;detx,dety,detz];
    else
        matDetPos = [matDetPos;0,0,0];
    end

end
[path] = uigetdir()
for num = 1:numel(DOT{currentsub}.zone.plotLst)
    tic
    DOT{currentsub}.plot =  DOT{currentsub}.zone.plot{num};
    DOT{currentsub}.plotLst = DOT{currentsub}.zone.plotLst{num};
    set(guiHOMER,'UserData',DOT);
    updataallview
    coord = [matSrcPos(DOT{currentsub}.plot(:,1),:);matDetPos( DOT{currentsub}.plot(:,2),:)];
    PositionFuite = mean(coord)*0.2
    axes(handles.axes_Mtg2)
    % center = get_Center( ans )
    set(handles.axes_Mtg2,'CameraTarget',[0,0,0]);
    set(handles.axes_Mtg2,'CameraPosition',PositionFuite);
    set(handles.axes_Mtg2,'CameraUpVector',[0,0,1]);
    titlefigure = [path,DOT{currentsub}.zone.label{num}]
    btn_exportHelmetFigure_Callback([], [],handles,titlefigure)
    toc
end


% --- Executes on button press in btn_saveview.
function btn_saveview_Callback(hObject, eventdata, handles)
% hObject    handle to btn_saveview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

angle = get(handles.edit_orientation,'string');
numwindows  = size(str2num(angle),1);
a_cote = [];
a_face = [];
h = handles.IO_HelmetMTG
rem=angle;
[name,path] = uiputfile();
if name==0
    return
end
for i = 1:numwindows
    [token, rem] = strtok(rem,',');
    if isempty(rem)
        a_cote =[a_cote 90];
        a_face = [a_face,0];
    else
        a_cote = [a_cote,str2num(token)];
        [token, rem] = strtok(rem(2:end),';');
        a_face = [a_face,str2num(token)];
    end
end
dist =3;
axisfactor = 0.15;
for i=0:numwindows-1
    tic
    %     h = figure;
    titlefigure = [path,name,'view',sprintf('%02.0f',a_cote(i+1)),sprintf('%02.0f',a_face(i+1))];
    angle_cote = a_cote(i+1);
    angle_face = a_face(i+1);
    cp_x = cos(angle_cote*pi/180) * cos(angle_face*pi/180)*dist;
    cp_y = sin(angle_cote*pi/180)* cos(angle_face*pi/180)*dist;
    cp_z = sin(angle_face*pi/180)*dist;
    set(handles.axes_Mtg2,'CameraTarget',[0,0,0]);
    set(handles.axes_Mtg2,'CameraPosition',[cp_x,cp_y,cp_z]);
    set(handles.axes_Mtg2,'CameraUpVector',[0,0,1]);
    axis([-1,1, -1,1,-1,1].*axisfactor)
    btn_exportHelmetFigure_Callback([], [],handles,titlefigure)
    time=toc;
    disp(['Computation time = ' num2str(time) 'seconds'])
end

function edit_orientation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_orientation as text
%        str2double(get(hObject,'String')) returns contents of edit_orientation as a double


% --- Executes during object creation, after setting all properties.
function edit_orientation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in radio_IDnb.
function radio_IDnb_Callback(hObject, eventdata, handles)
% hObject    handle to radio_IDnb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_IDnb
set(handles.radio_viewdistance,'value',0)
set(handles.radio_nbavg,'value',0')
resetview(handles)

% --- Executes during object creation, after setting all properties.
function axes_Mtg2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_Mtg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_Mtg2




% --- Executes on button press in radio_nbavg.
function radio_nbavg_Callback(hObject, eventdata, handles)
% hObject    handle to radio_nbavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_nbavg
set(handles.radio_IDnb,'value',0)
set(handles.radio_viewdistance,'value',0)
resetview(handles)


% --- Executes on button press in radio_viewch.
function radio_viewch_Callback(hObject, eventdata, handles)
% hObject    handle to radio_viewch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_viewch
set(handles.radio_d1color, 'value', 0);
resetview(handles)



% --- Executes on button press in btn_scale.
function btn_scale_Callback(hObject, eventdata, handles)
% hObject    handle to btn_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of btn_scale

resetview(handles)


% --------------------------------------------------------------------
function menu_exportfigureoption_Callback(hObject, eventdata, handles)
% hObject    handle to menu_exportfigureoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ExportFigureOption


% --------------------------------------------------------------------
function menu_multidisplaymovie_Callback(hObject, eventdata, handles)
% hObject    handle to menu_multidisplaymovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_figure = figure
b = uicontrol('Style','pushbutton','String','Play','units','normalize',...
    'Position',[0.05 0.9 0.1 0.1],'Tag','pushbutton','Callback',{@Multidisplay_IOhelmet,h_figure,1});
%Bouton pour le play image
%Label
b = uicontrol('Style','text','string','Start','units','normalize',...
    'Position',[0.05,0.8,0.1, 0.05],'HorizontalAlignment','left');
b = uicontrol('Style','text','string','Stop','units','normalize',...
    'Position',[0.05,0.75,0.1, 0.05],'HorizontalAlignment','left');
b = uicontrol('Style','text','string','Step','units','normalize',...
    'Position',[0.05,0.7,0.1, 0.05],'tag','text_windowstime','HorizontalAlignment','left');

b = uicontrol('Style','text','string','Nb view','units','normalize',...
    'Position',[0.05,0.6,0.3, 0.05],'HorizontalAlignment','left');
b = uicontrol('Style','text','string','90,30 => Left side with 30 degree','units','normalize',...
    'Position',[0.05,0.55,0.3, 0.05],'HorizontalAlignment','left');
b = uicontrol('Style','text','string','View angle','units','normalize',...
    'Position',[0.05,0.50,0.3, 0.05],'HorizontalAlignment','left');
b = uicontrol('Style','text','string','Axiszoom','units','normalize',...
    'Position',[0.05,0.45,0.1,0.05],'HorizontalAlignment','left');
b = uicontrol('Style','text','string','Space','units','normalize',...
    'Position',[0.05,0.40,0.1,0.05],'HorizontalAlignment','left');
b = uicontrol('Style','text','string','Speed','units','normalize',...
    'Position',[0.05,0.35,0.1,0.05],'HorizontalAlignment','left');
b=uicontrol('Style','radiobutton','string','use avg vcolor','tag','radio_avgvcolor',...
    'units','normalize', 'Position',[0.05,0.30,0.2,0.05],'HorizontalAlignment','left')

b=uicontrol('Style','edit','string','HbOAvg','tag','edit_NameFile',...
    'units','normalize', 'Position',[0.30,0.30,0.4,0.05],'HorizontalAlignment','left')
b=uicontrol('Style','radiobutton','string','use .mat topographic file','tag','radio_avgmat',...
    'units','normalize', 'Position',[0.05,0.25,0.2,0.05],'HorizontalAlignment','left')
%Edit
b = uicontrol('Style','edit','string','0','units','normalize',...
    'Position',[0.15,0.8,0.1, 0.05],'tag','edittimestart');
b = uicontrol('Style','edit','string','1','units','normalize',...
    'Position',[0.15,0.75,0.1, 0.05],'tag','edittimestop');
b = uicontrol('Style','edit','string','1','units','normalize',....
    'Position',[0.15,0.7,0.1, 0.05],'tag','edit_steptime');

b = uicontrol('Style','edit','string','2','units','normalize',...
    'Position',[0.15,0.6,0.1, 0.05],'tag','editwindows');
b = uicontrol('Style','edit','string','90,30;-90,30;0,50;180,0','units','normalize',...
    'Position',[0.15,0.50,0.3, 0.05],'tag','editviewangle');
b = uicontrol('Style','edit','string','0.13','units','normalize',... %avant arriere
    'Position',[0.15,0.45,0.1, 0.05],'tag','edit_xaxiszoom');
b = uicontrol('Style','edit','string','0.13','units','normalize',... %coté
    'Position',[0.25,0.45,0.1, 0.05],'tag','edit_yaxiszoom');
b = uicontrol('Style','edit','string','0.13','units','normalize',... %hauteur
    'Position',[0.35,0.45,0.1, 0.05],'tag','edit_zaxiszoom');

b = uicontrol('Style','edit','string','0.2','units','normalize',...
    'Position',[0.15,0.40,0.1, 0.05],'tag','editspace');
b = uicontrol('Style','edit','string','0.5','units','normalize',...
    'Position',[0.15,0.35,0.1, 0.05],'tag','editspeed');



function ExportMRIShape_Callback(hObject, eventdata, handles)

PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
oMRI = get_MRI_Data(PrjStruct);
if get(handles.radio_guiSPMnirsHSJ,'value')==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end
PMI = get(guiHOMER,'UserData');
HOMERhandles = guihandles(guiHOMER);
DispHelm = get_Helmet( PrjStruct );
vHoles = get_vHoles(DispHelm );
sMtg = get_Mtg( DispHelm );
p = [sMtg.v_pSrc, sMtg.v_pDet];
IRM = logical(zeros(181,217,181));
%p = [5,11];
%    [90,90,108];
%  for ip = 1:numel(p)
%      for jp = 1:numel(p)
%                xi = round([-vHoles(p(ip)).Coord.y,  vHoles(p(ip)).Coord.x,  vHoles(p(ip)).Coord.z] * 1000+[90,60,108]);
%                xj = round([-vHoles(p(jp)).Coord.y,  vHoles(p(jp)).Coord.x,  vHoles(p(jp)).Coord.z] * 1000+[90,60,108]);
%                if xi(1)<xj(1);x = xi(1):xj(1);else x=xi(1):-1:xj(1);end
%                if xi(2)<xj(2); y = xi(2):xj(2);else y=xi(2):-1:xj(2);end
%                if xi(3)<xj(3); z = xi(3):xj(3);else z = xi(3):-1:xj(3);end
%                    for i =x
%                        for j=y
%                            for k=z
%                                     IRM(i,j,k)=1;
%                            end
%                        end
%                    end
%        end
%  end

for ip = 1:numel(p)
    xi(ip,:) = round([-vHoles(p(ip)).Coord.y,  vHoles(p(ip)).Coord.x,  vHoles(p(ip)).Coord.z] * 1000);
end

max(xi(:,2))

fid = fopen('test2.img','w')
fwrite(fid,IRM,'uint8')
fclose(fid)
%
%   MRI_Data
%      fid = [91,35,6,
%      174, 27,114,
%      5,27,114]
%
% AP = fid(2,3)-fid(1,3);
% BH = max(xi(:,3));   %fid(2,2) % centre
% LR = round((fid(2,1)-fid(3,1))/2);
% MRIx= fid(3,1)+LR
% MRIy= fid(3,3)
% MRIz=fid(2,2)
% [x,y,z] =ellipsoid(MRIx,MRIy,MRIz,LR,AP,BH,40);
% plot3(x,y,z)
%                    for i = round(x)
%                        for j=round(y)
%                            for k=round(z)
%                                     IRM(i,j,k)=1;
%                            end
%                        end
%                    end
%  fid = fopen('test2.img','w')
%  fwrite(fid,IRM,'uint8')
%  fclose(fid)


% --- Executes on button press in radio_show_cover.
function radio_show_cover_Callback(hObject, eventdata, handles)
% hObject    handle to radio_show_cover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_show_cover
resetview(handles)

% --------------------------------------------------------------------
function menu_open_Topomat_Callback(hObject, eventdata, handles)
% hObject    handle to menu_open_Topomat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]= uigetfile('*Topo.mat');
if path==0
    set(handles.edit_TOPOmatfile,'string',['']);
else
    set(handles.edit_TOPOmatfile,'string',[path,file]);
end



function edit_TOPOmatfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TOPOmatfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TOPOmatfile as text
%        str2double(get(hObject,'String')) returns contents of edit_TOPOmatfile as a double


% --- Executes during object creation, after setting all properties.
function edit_TOPOmatfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TOPOmatfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_hupdate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hupdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hupdate as text
%        str2double(get(hObject,'String')) returns contents of edit_hupdate as a double


% --- Executes during object creation, after setting all properties.
function edit_hupdate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hupdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_IDnumber_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IDnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IDnumber as text
%        str2double(get(hObject,'String')) returns contents of edit_IDnumber as a double


% --- Executes during object creation, after setting all properties.
function edit_IDnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IDnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D1matrix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_D1matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_D1matrix as text
%        str2double(get(hObject,'String')) returns contents of
%        edit_D1matrix as a double


% --- Executes during object creation, after setting all properties.
function edit_D1matrix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_D1matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popupmenu_backgroundcolor.
function popupmenu_backgroundcolor_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_backgroundcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_backgroundcolor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_backgroundcolor
resetview(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_backgroundcolor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_backgroundcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

    % This panel contains the selection tools.
    cameratoolbar('SetMode', 'nomode'); % Deselect the camera controls.


% --- Executes on selection change in popup_channel_style.
function popup_channel_style_Callback(hObject, eventdata, handles)
% hObject    handle to popup_channel_style (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_channel_style contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_channel_style

% This callback is a shortcut used to keep the previously used radio
% buttons. Their value is set here and then read by other parts of the
% program.

    switch get(handles.popup_channel_style, 'value')
        case 1 % Default
            set(handles.radio_viewch, 'value', 0)
            set(handles.radio_d1color, 'value', 0)
        case 2 % Black
            set(handles.radio_viewch, 'value', 1)
            set(handles.radio_d1color, 'value', 0)
        case 3 % Heatmap
            set(handles.radio_viewch, 'value', 0)
            set(handles.radio_d1color, 'value', 1)
    end
    resetview(handles)


% --- Executes during object creation, after setting all properties.
function popup_channel_style_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_channel_style (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_labelch.
function popup_labelch_Callback(hObject, eventdata, handles)
% hObject    handle to popup_labelch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_labelch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_labelch

resetview(handles)
% --- Executes during object creation, after setting all properties.
function popup_labelch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_labelch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uipanel10_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

resetview(handles)



function edit_channelwidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_channelwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_channelwidth as text
%        str2double(get(hObject,'String')) returns contents of edit_channelwidth as a double
resetview(handles)


% --- Executes during object creation, after setting all properties.
function edit_channelwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_channelwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function context_topo_radial_Callback(hObject, eventdata, handles)
% hObject    handle to context_topo_radial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.context_topo_IWD,'checked','off')
set(handles.context_topo_radial,'checked','on')

% --------------------------------------------------------------------
function context_topo_IWD_Callback(hObject, eventdata, handles)
% hObject    handle to context_topo_IWD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.context_topo_IWD,'checked','on')
set(handles.context_topo_radial,'checked','off')
% --------------------------------------------------------------------
function Topographymode_Callback(hObject, eventdata, handles)
% hObject    handle to Topographymode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
