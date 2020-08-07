function varargout = plot_sessions_GUI(varargin)
% PLOT_SESSIONS_GUI M-file for plot_sessions_GUI.fig
%      PLOT_SESSIONS_GUI, by itself, creates a new PLOT_SESSIONS_GUI or raises the existing
%      singleton*.
%
%      H = PLOT_SESSIONS_GUI returns the handle to a new PLOT_SESSIONS_GUI or the handle to
%      the existing singleton*.
%
%      PLOT_SESSIONS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_SESSIONS_GUI.M with the given input arguments.
%
%      PLOT_SESSIONS_GUI('Property','Value',...) creates a new PLOT_SESSIONS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plot_sessions_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plot_sessions_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plot_sessions_GUI

% Last Modified by GUIDE v2.5 24-Jul-2020 13:51:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plot_sessions_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @plot_sessions_GUI_OutputFcn, ...
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

% --- Executes just before plot_sessions_GUI is made visible.
function handles = plot_sessions_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plot_sessions_GUI (see VARARGIN)

% Choose default command line output for plot_sessions_GUI
handles.output = hObject;
try
handles.NIRSall = varargin{1};
handles.subjectnb = varargin{2};
job = varargin{3};
set(handles.edit_nirsmat, 'String', handles.NIRSall); %Set the file pop-up menu
%job.NIRSmat{filenb,1}
handles.NIRSpath = handles.NIRSall;
load(handles.NIRSpath{handles.subjectnb,1});
handles.NIRS = NIRS;
handles.oldfile = 1;
handles.oldmodule = 1;
%Get the files' names to make the list for the pop-up menu
nb_file = length(handles.NIRS.Dt.fir.pp(end).p);
handles.filepath = [];
module_list = [];
for i = 1:nb_file
    [dir1,fil1,ext1]= fileparts(NIRS.Dt.fir.pp(end).p{i,1});
    handles.filepath = [handles.filepath {['Blocks' sprintf('%03.0f',i),' ',fil1]}];
end
handles.module_list = [];
handles.module_path = [];
for i = 1:length(NIRS.Dt.fir.pp)
    if isempty(NIRS.Dt.fir.pp(i).pre)
        handles.module_list{i} = 'null';
    else
        handles.module_list{i} = [NIRS.Dt.fir.pp(i).pre, ' off'];
    end
    handles.module_path{i} = NIRS.Dt.fir.pp(i).p{1};
end


[s,~] = listdlg('PromptString','Select the module where you will apply manual modification :',...
    'SelectionMode','single',...
    'ListString',handles.module_list,...
    'ListSize',[400,200],...
    'InitialValue',numel(handles.module_list));
lst = numel(handles.module_list);
%NewDirCopyNIRS = 0
[dir1,fil1,ext1] = fileparts(NIRS.Dt.fir.pp(s).p{1});
[dir2,fil2,ext2] = fileparts(handles.NIRSpath{handles.subjectnb,1});
        
if isempty(strfind(NIRS.Dt.fir.pp(s).pre,'Manual Gui'))|~strcmp(upper(dir1),upper(dir2)) %create a new manual field of data where modification will be saved
    NC = NIRS.Cf.H.C.N;
    rDtp = NIRS.Dt.fir.pp(s).p;
    prefix = 'm';
    for f = 1:numel(rDtp)
        try
            d = fopen_NIR(rDtp{f},NC);
        catch
            msgbox(['ERROR Please check file ',rDtp{f} ])
            disp(['ERROR Please check file ',rDtp{f} ,' you can used Utility / folder ajustement to modify the NIRS.mat at the new location.'])
            return
        end
        [dir1,fil1,ext1] = fileparts(rDtp{f});
        [dir2,fil2,ext2] = fileparts(handles.NIRSpath{handles.subjectnb,1});
        if isfield(job,'DelPreviousData')
            DelPreviousData  = job.DelPreviousData;
        else
            DelPreviousData = 0 ;
        end
        %         if isfield(job,'NewDirCopyNIRS')
        %             if isfield(job.NewDirCopyNIRS,'CreateNIRSCopy')
        %                 NewNIRSdir = job.NewDirCopyNIRS.CreateNIRSCopy.NewNIRSdir;
        %                 disp(['Create directory for condition ',NewNIRSdir])
        %                 NewDirCopyNIRS = 1;
        %             else
        %                 NewDirCopyNIRS = 0;
        %             end
        %         else
        %             NewDirCopyNIRS = 0;
        %         end
        infilenir = fullfile(dir1,[fil1 '.nir']);
        infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
        infilevhdr = fullfile(dir1,[fil1 '.vhdr']);
        %         if NewDirCopyNIRS
        %              dir2 = [dir1 filesep NewNIRSdir];
        if ~exist(dir2,'dir'), mkdir(dir2); end;
        outfile = fullfile(dir2,[prefix fil1 ext1]);
        outfilevmrk = fullfile(dir2,[prefix fil1  '.vmrk']);
        outfilevhdr = fullfile(dir2,[prefix fil1  '.vhdr']);
        %             else
        %                 outfile = fullfile(dir1,[prefix fil1 ext1]);
        %                 outfilevmrk = fullfile(dir1,[prefix fil1  '.vmrk']);
        %                 outfilevhdr = fullfile(dir1,[prefix fil1  '.vhdr']);
        %             end
        
        NIRS.Dt.fir.pp(lst+1).pre = ['Manual Gui ', NIRS.Dt.fir.pp(s).pre];
        NIRS.Dt.fir.pp(lst+1).job =  NIRS.Dt.fir.pp(s).job;
        NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
        fwrite_NIR(outfile,d);
        handles.file_vmrk =  outfilevmrk;
        [label,ind_dur_ch]=read_vmrk_all(infilevmrk);
        write_vmrk_all(outfilevmrk,ind_dur_ch,label);
        %vhdr copy
            %write outvhdr file
           ChannelLabels= ConvertmlIDsrs2label(NIRS);
%     for id=1:size(NIRS.Cf.H.C.id,2) 
%          if strcmp(NIRS.Cf.dev.n,'NIRx') %NIRx
%             srs = SDPairs2strboxy(NIRS.Cf.H.C.id(1,id));
%             det = SDDet2strboxy(NIRS.Cf.H.C.id(1,id));
%             ChannelLabels{id,1} = [srs, '_', det];
%          else                             %On ajoute la nomenclature Ste-Justine 
%             srs = SDPairs2strboxy_ISS(NIRS.Cf.H.C.id(1,id));
%             det = SDDet2strboxy_ISS(NIRS.Cf.H.C.id(1,id));
%             ChannelLabels{id,1} = [srs, '_', det];
%         end
%     end
    SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
    nirs_boxy_write_vhdr(outfilevhdr,... %Output file
        outfile,... %DataFile
        outfilevmrk,... %MarkerFile,...
        'nirs_Concatenate_File',... %Function that created the header
        '',... %Channel Resolution
        '',... %Channel Units
        ChannelLabels,... %names given as a column of cells
        SamplingInterval,...%SamplingInterval in microseconds
        size(d,2)); %Number Sample in microseconds
     %UPDATE SelectedFactors to the last module
    try
     load(fullfile(dir2,'SelectedFactors.mat'));
     for ifactor=1:numel(PARCOMP)
        PARCOMP(ifactor).module =  PARCOMP(ifactor).module+1;
        PARCOMP(ifactor).modulestr= ['Manual Gui ', NIRS.Dt.fir.pp(s).pre];
     end
     save(fullfile(dir2,'SelectedFactors.mat'),'PARCOMP');
    catch
    end
        
        
        try
            copyfile(infilevhdr,outfilevhdr)
        catch
        end
        if DelPreviousData
            delete(infilenir);
            delete(infilevmrk);
            delete(infilevhdr);
        end
        
    end
    handles.module_list{lst+1} = ['Manual Gui ', NIRS.Dt.fir.pp(s).pre];
end
handles.NIRS = NIRS;
guidata(handles.figure1, handles);
%   if NewDirCopyNIRS
%       [dir1,fil1,ext1] =fileparts(handles.NIRSpath{handles.subjectnb,1})
%       handles.NIRSpath{handles.subjectnb,1}=fullfile([dir1 filesep NewNIRSdir],'NIRS.mat');
%       save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');
%   else
save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');
% end
set(handles.popupmenu_file, 'String', handles.filepath); %Set the file pop-up menu
set(handles.popupmenu_file, 'Value', 1);
set(handles.popupmenu_module, 'String', handles.module_list); %Set the module pop-up menu
set(handles.popupmenu_module, 'value', numel(handles.module_list));
set(handles.popupmenu_module_hold, 'String', handles.module_list); %Set the module pop-up menu
set(handles.popupmenu_module_hold, 'value', numel(handles.module_list));
if isfield(NIRS,'SPM')
    set(handles.radio_GLMresult,'enable','on')
else
    set(handles.radio_GLMresult,'enable','off')
end

%Get triggers
idfile = get(handles.popupmenu_file,'value');
if isfield(handles.NIRS.Dt.fir,'aux5')
    handles.triggers = handles.NIRS.Dt.fir.aux5{idfile};
    
    unique_trig = unique(handles.triggers(:,1));
    %for i=1:numel(unique_trig)
    trig_list = [];
    for i=1:numel(unique_trig)
        trig_list = [trig_list;  {[num2str(unique_trig(i)),'  on']}];
    end
    set(handles.popupmenu8, 'String', trig_list);
    set(handles.popupmenu8, 'value', numel(trig_list));
    
end
%Compatibilite avec le helmet et display IOSTJ
%Utilisation de certain champ de Homer
setappdata(0,'gui_SPMnirsHSJ',handles.figure1);
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
global currentsub;
currentsub = 1;
PMI{currentsub}.currentFile = 1;
cf = 1;
if isfield(NIRS.Cf.H,'p')
    PMI{currentsub}.prj_name = NIRS.Cf.H.p{1};
end
PMI{currentsub}.data(cf).MeasList = [NIRS.Cf.H.C.id(2:3,:)',...
    ones(size(NIRS.Cf.H.C.id,2),1),...
    [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
fileid=get(handles.popupmenu_file,'value');
% if isfield(handles.NIRS.Cf.H.C, 'okavg')
%     PMI{currentsub}.data(cf).MeasListAct = handles.NIRS.Cf.H.C.okavg;
% else 
     PMI{currentsub}.data(cf).MeasListAct = handles.NIRS.Cf.H.C.ok(:, fileid);
% end
PMI{currentsub}.color = PMIcolordef(size(NIRS.Cf.H.C.id,2));
PMI{currentsub}.plotLst = [1];
PMI{currentsub}.plot = [1,1]; 
set(handles.text4,'string',['Channels HbO 1:' num2str(numel(PMI{currentsub}.data(cf).MeasListAct)/2)]);
[pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(1,1).p{1,1});
try
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    load(fullfile(pathstr,'CorrectionApply.mat'));
    for i=1:numel(PARCORR)
        CORRlist{i} = PARCORR(i).label;
    end
    set(handles.listbox_CorrectionDecomposition,'string',CORRlist)
    %   set(handles.listbox_CorrectionDecomposition,'value',numel(CORRlist))
catch
    set(handles.listbox_CorrectionDecomposition,'string','')
    set(handles.listbox_CorrectionDecomposition,'value',1)
end
try
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    load(fullfile(pathstr,'SelectedFactors.mat'))
    for i=1:numel(PARCOMP)
        COMPlist{i} = PARCOMP(i).label;
    end
    set(handles.listbox_Component,'string',COMPlist);
     set(handles.listbox_Component,'value',1);
catch
    set(handles.listbox_Component,'string','')
    set(handles.listbox_Component,'value',1)
end
try %Check if eeg field
    if isfield(NIRS.Dt,'EEG')
        set(handles.gui_EEG,'enable','on')
        flag_stepnormalisedmissing = 1;
        for imodule = 1:numel(NIRS.Dt.fir.pp)
            if strcmp(NIRS.Dt.fir.pp(imodule).pre,'Segmentation')
                flag_stepnormalisedmissing = 0;
            end
        end
        if flag_stepnormalisedmissing
            set(handles.text_WarningAxesPhysiologie,'visible','on')
            set(handles.text_WarningAxesPhysiologie,'string', 'Use segment to align trig using multimodal data')
        else
            set(handles.text_WarningAxesPhysiologie,'visible','off')
            %msgbox('EEG time have not been synchronize using trigger, used the normalization to segment and synchronised NIRS and EEG according to trigger info')
        end
        
    else
        set(handles.gui_EEG,'enable','off')
    end
catch
end
if isfield(NIRS.Dt,'AUX')
    set(handles.listbox_AUXCH,'enable','on')
    flag_stepnormalisedmissing = 1;
    for imodule = 1:numel(NIRS.Dt.fir.pp)
        NIRS.Dt.fir.pp(imodule).pre;
        if strcmp(NIRS.Dt.fir.pp(imodule).pre,'Segmentation')
            flag_stepnormalisedmissing = 0;
        end
    end
    if flag_stepnormalisedmissing
        set(handles.text_WarningAxesPhysiologie,'visible','on')
        set(handles.text_WarningAxesPhysiologie,'string', 'Use segment to align trig using multimodal data')
    else
        set(handles.text_WarningAxesPhysiologie,'visible','off')
        %msgbox('EEG time have not been synchronize using trigger, used the normalization to segment and synchronised NIRS and EEG according to trigger info')
    end
    
    popuplabel = [];
    for iAUX= 1:numel(NIRS.Dt.AUX)
        popuplabel = [ popuplabel;{NIRS.Dt.AUX(iAUX).label}];
    end
    set(handles.popupauxnb,'string',popuplabel);
else
    set(handles.listbox_AUXCH,'enable','off')
end
if isfield(NIRS.Dt,'Video')
    set(handles.btn_Video,'enable','on');
    flag_stepnormalisedmissing = 1;
    for imodule = 1:numel(NIRS.Dt.fir.pp)
        if strcmp(NIRS.Dt.fir.pp(imodule).pre,'Segmentation')
            flag_stepnormalisedmissing = 0;
        end
    end
    if flag_stepnormalisedmissing
        set(handles.text_WarningAxesPhysiologie,'visible','on');
        set(handles.text_WarningAxesPhysiologie,'string', 'Use segment to align trig using multimodal data');
    else
        set(handles.text_WarningAxesPhysiologie,'visible','off');
    end
else
    set(handles.listbox_AUXCH,'enable','off');
    set(handles.btn_Video,'enable','off');
end
if isfield(NIRS.Dt,'Audio')
    set(handles.btn_AUDIO,'enable','on');
else
    set(handles.btn_AUDIO,'enable','off');
end
if isfield(NIRS.Dt,'Zone')%to be implement save in the NIRS.mat
    set(handles.popupmenu_zone,'string',{''});
    set(handles.popupmenu_zone,'value',1);
else
    set(handles.popupmenu_zone,'string',{''});
    set(handles.popupmenu_zone,'value',1);
end

p = mfilename('fullpath');
[rem,tok]=strtok(fliplr(p),{'','\'});
icon= imread(fullfile(fliplr(tok),'LoupeZoom.jpg'));
set(handles.btn_zoom_timestart_timestop,'cdata',icon );
icon = imread(fullfile(fliplr(tok),'SortArrow.jpg'));
set(handles.btn_SortComponent,'cdata',icon);
icon = imread(fullfile(fliplr(tok),'SortArrow.jpg'));
set(handles.btn_sort_correction,'cdata',icon);
icon = imread(fullfile(fliplr(tok),'LoupeZoomout.jpg'));
set(handles.btn_zoomout,'cdata',icon);

set(guiHOMER,'UserData',PMI);
updatedata(handles,1,1);
updatedisplay(handles);
guidata(handles.figure1, handles);


catch
    guidata(hObject, handles);
end

% UIWAIT makes plot_sessions_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function [handles] = set_popupmenus(handles)
% Sets the NIRS.mat path box, the file pop-up menu and the module pop-up
% menu

%Get the files' names to make the list for the pop-up menu
nb_file = length(handles.NIRS.Dt.fir.pp(1).p);
handles.filepath = [];
module_list = [];
for i = 1:nb_file
    handles.filepath = [handles.filepath NIRS.Dt.fir.pp(1).p(i)];
end
handles.module_list = [];
handles.module_path = [];
for i = 1:length(NIRS.Dt.fir.pp)
    handles.module_list{i} = NIRS.Dt.fir.pp(i).pre;
    handles.module_path{i} = NIRS.Dt.fir.pp(i).p{1};
end
set(handles.popupmenu_file, 'Value', nb_file);
set(handles.popupmenu_file, 'String', handles.filepath); %Set the file pop-up menu
set(handles.popupmenu_module, 'String', handles.module_list); %Set the module pop-up menu



% --- Outputs from this function are returned to the command line.
function varargout = plot_sessions_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('NIRS.mat');
[path,file]
fileopen = get(handles.edit_nirsmat,'string');

% varargin{1,1}= fileopen; 
% varargin{1,2} = '1';
% varargin{1,3} = 'job'
handles = plot_sessions_GUI_OpeningFcn(handles.figure1, [], handles, fileopen,  '1', 'job')


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

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

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on selection change in popupmenu_file.
function popupmenu_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.chosenfile = get(handles.popupmenu_file, 'Value');
oldfile   = handles.oldfile;
set(handles.popupmenu_file,'value',oldfile);
handles = btn_savenoise_Callback(0,0,handles);
set(handles.popupmenu_file,'value',handles.chosenfile);


lst = length(handles.NIRS.Dt.fir.pp); %get the number of processing step
%Get all the processing steps (modules)
handles.module_list = [];
handles.module_path = [];
for i = 1:lst
    handles.module_list{i} = handles.NIRS.Dt.fir.pp(i).pre;
    if numel(handles.NIRS.Dt.fir.pp(i).p)>=handles.chosenfile
        handles.module_path{i} = handles.NIRS.Dt.fir.pp(i).p{handles.chosenfile};
    else
        handles.module_path{i} = handles.NIRS.Dt.fir.pp(i).p{1};
    end
end
% set(handles.popupmenu_module,'string',handles.module_list); JT

trig_ind = get(handles.popupmenu8, 'value');
idfile = get(handles.popupmenu_file,'value');
if idfile ~= handles.oldfile
    if isfield(handles.NIRS.Dt.fir,'aux5')
        if numel(handles.NIRS.Dt.fir.aux5)>=idfile
            handles.triggers = handles.NIRS.Dt.fir.aux5{idfile};
            if ~isempty(handles.triggers)
                handles.triggers(find(handles.triggers(:,2)==0),2)=1;
                unique_trig = unique(handles.triggers(:,1));
                trig_list = [];
                for i=1:numel(unique_trig)
                    trig_list = [trig_list;  {[num2str(unique_trig(i)),'  on']}];
                end
                set(handles.popupmenu8, 'String', trig_list);
                if numel(trig_list)< trig_ind
                    set(handles.popupmenu8, 'value', numel(trig_list));
                end
            else
                handles.triggers = [];
            end
        else
            handles.triggers = [];
        end
    end
    try
        tmp = [];
        for itrig = 1:numel(trig_list)
            trigtype =trig_list{itrig};
            numtrig = str2num(trigtype(1:end-4));
            if  strcmp(trigtype(end-3:end),'  on')
                tmp = [tmp;handles.triggers(ismember(handles.triggers(:,1), numtrig),:)];
            end
            
            %handles.selected_trig(itrig,:) = handles.triggers(ismember(handles.triggers(:,1), str2num(trig_name{trig_ind(itrig)})),:);
        end
        handles.selected_trig=tmp;
    catch
    end
end


updatedata(handles,1,1);
handles.newlist = 1; %Reactualiser la liste de canaux de 3DHelmet
updatedisplay(handles);
handles.oldfile = handles.chosenfile; %get(handles.popupmenu_file, 'Value');
guidata(hObject, handles);



% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_file

% --- Executes during object creation, after setting all properties.
function popupmenu_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nirsmat_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nirsmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(hObject,'String');

handles.subjectnb = get(handles.edit_nirsmat,'value');
% var{1} = handles.NIRSpath
% var{2} = handles.subjectnb

handles = plot_sessions_GUI_OpeningFcn(hObject, [], handles, handles.NIRSpath, handles.subjectnb);
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit_nirsmat as text
%        str2double(get(hObject,'String')) returns contents of edit_nirsmat as a double


% --- Executes during object creation, after setting all properties.
function edit_nirsmat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nirsmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_module.
function popupmenu_module_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Ssave the data noise before if last module manual mode!!!
idmodule = get(handles.popupmenu_module,'value'); %the new module to display
if handles.oldmodule==numel(get(handles.popupmenu_module,'string'))
    disp('Save last module')
    set(handles.popupmenu_module,'value',handles.oldmodule)
    handles = btn_savenoise_Callback(0,0,handles);
end
set(handles.popupmenu_module,'value', idmodule ); %the new module to display

%avg option
set(handles.popup_average,'enable','off');
set(handles.popup_average,'value',1);
if numel(handles.NIRS.Dt.fir.pp(idmodule).pre)>15
    if strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre(1:15), 'Epoch averaging')
        set(handles.popup_average,'enable','on');
    end
end

nb_file = length(handles.NIRS.Dt.fir.pp(idmodule).p);
handles.filepath = [];
for i = 1:nb_file
    [dir1,fil1,ext1]= fileparts(handles.NIRS.Dt.fir.pp(idmodule).p{i});
    handles.filepath = [handles.filepath {['Blocks' sprintf('%03.0f',i)]}];
end
set(handles.popupmenu_file,'string',handles.filepath)
if get(handles.popupmenu_file,'value') > numel(handles.filepath)
    set(handles.popupmenu_file,'value',1)
end
guidata(hObject, handles);
updatedata(handles,1,1)
updatedisplay(handles)

handles.oldmodule = get(handles.popupmenu_module,'value');
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_module contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_module


% --- Executes during object creation, after setting all properties.
function popupmenu_module_CreateFcn(hObject, eventdata, ~)
% hObject    handle to popupmenu_module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_plotLst_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plotLst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

list = str2num(get(hObject,'String'));

% Placer la sélection listbox_normlist à la valeur du premier canal
% sélectionné
idfile = get(handles.popupmenu_file,'value');
imodulenorm = [];
for i=1:numel(handles.NIRS.Dt.fir.pp)
    if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
        if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
            imodulenorm  = i;        %Le dernier is plusieur Step Detection
        end
    end
end
%Prendre la valeur 830 pour chercher dans la list normlist
ch1 = list(1);
if ch1 <= handles.NIRS.Cf.H.C.N/2
    ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
else
    ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
end
Ch = sort([ch1,ch2]);
try
    if ~isempty(imodulenorm)
        NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist(idfile,1);
        if iscell(NormList)
            NormList = NormList{1};
        end
        if ~isempty(NormList)
            ind = find(NormList(:,1)==Ch(1));
            set(handles.listbox_normlist,'value',ind(1));
        end
    end
catch
end

%Update de l'interface avec la nouvelle liste
global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf  = 1;
PMI{currentsub}.plotLst = list;
PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(list,1),PMI{currentsub}.data(cf).MeasList(list,2)];
handles.newlist=1; %update helmet channel
set(guiHOMER,'UserData',PMI)
guidata(hObject, handles);
updatedisplay(handles)



% Hints: get(hObject,'String') returns contents of edit_plotLst as text
%        str2double(get(hObject,'String')) returns contents of edit_plotLst as a double


% --- Executes during object creation, after setting all properties.
function edit_plotLst_CreateFcn(hObject, eventdata, ~)
% hObject    handle to edit_plotLst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in radio_triggers.
function radio_triggers_Callback(hObject, eventdata, handles)
% hObject    handle to radio_triggers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatedisplay(handles);


% Hint: get(hObject,'Value') returns toggle state of radio_triggers



function [rho,handle_out] = corr_disp(handles)
% Process the correlation between coupled channels for the ones that are being displayed.
%

handles.plotlist_corr = handles.plotlist(handles.plotlist <= handles.NC/2);

fs = handles.NIRS.Cf.dev.fs;
win = floor(2.5*fs); %2x win = 5s moving window
samp_length = size(handles.d,2);
rho = [];

hwaitbar = waitbar(0);
for i = 1:length(handles.plotlist_corr) %Loop over all selected channel couples
    waitbar(i/length(handles.plotlist_corr),hwaitbar,'Processing correlations...');
    Idx = handles.plotlist_corr(i);
    d1 = handles.d(Idx,:)';
    d2 = handles.d(Idx+handles.NC/2,:)';
    for t = 1:win
        rho(Idx,t) = corr(d1(1:t*2),d2(1:t*2));
    end
    for t = 1+win:samp_length-win
        rho(Idx,t) = corr(d1(t-win:t+win),d2(t-win:t+win));
    end
    for t = size(handles.d,2)-win+1:samp_length
        rho(Idx,t) = corr(d1(samp_length-(samp_length-t)*2:samp_length),d2(samp_length-(samp_length-t)*2:samp_length));
    end
end
close(hwaitbar);
handle_out = handles;


% --- Executes on button press in btn_View.
function btn_View_Callback(hObject, eventdata, handles)
% hObject    handle to btn_View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
option.conc = 1;
option.SPMgui = 1;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
try
    PMI{currentsub}.prj_name=handles.NIRS.Dt.fir.pp(1).job.prjfile{1};
catch
end
set(guiHOMER,'UserData',PMI);
set(handles.radio_3dmontageupdate,'value',1)
IO_HelmetMTG_Display(handles,option)

function updatedisplay(handles)
%try
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
PMI{currentsub}.currentFile = 1;
cf = PMI{currentsub}.currentFile;
set(handles.text_WarningAxesPhysiologie,'visible','off')
%Show rejected channels
try
    plotLst_ok = PMI{currentsub}.plotLst; %Get the actual plotting list
catch
    plotLst_ok = 1;
    PMI{currentsub}.plotLst = 1;
end
if ~isfield(PMI{currentsub}.data(cf),'MeasListAct')
    bad_ch = [];
else
    bad_ch = find(~PMI{currentsub}.data(cf).MeasListAct);
end

%update trig status
trig_ind = get(handles.popupmenu8, 'value');
trig_name = get(handles.popupmenu8, 'string');
if isfield(handles,'triggers')
    if ~isempty(handles.triggers)
        %handles.NIRS.Dt.fir.aux5bloc{idfile};
        handles.selected_trig = handles.triggers(ismember(handles.triggers(:,1), str2num(char(trig_name{trig_ind}))),:);
    else
        set(handles.radio_triggers,'value',0)
    end
end
ind_bad = ismember(plotLst_ok,bad_ch');
plotLst_bad = plotLst_ok(ind_bad);
PMI{currentsub}.plotLst_bad = plotLst_bad;
set(guiHOMER,'UserData',PMI);


% %Show triggers
% if get(handles.radio_triggers,'value')
%     axes(handles.display_axes);
%     trig_ind = get(handles.popupmenu8, 'value');
%     plot_triggers(PMI{currentsub}.data(cf).HRF.AvgC, handles.triggers(trig_ind,:));
% end

%Update data
handles = plotAxes_d1(handles);
try
    handles = plotAxes_physiologie(handles);
catch
end
PMI = get(guiHOMER,'UserData');
if get(handles.radio_3dmontageupdate,'value')
guiHelmet = getappdata(0,'guiHelmet');
if ~isempty(guiHelmet)
    if isfield(handles,'newlist')
        if handles.newlist==1; % Update helmet if change in the list only
            Helmethandles = guihandles(guiHelmet);
            fhresetview = getappdata(guiHelmet,'fhresetview');
            fhresetview(Helmethandles);
            figure(handles.figure1);
        end
    else
        Helmethandles = guihandles(guiHelmet);
        fhresetview = getappdata(guiHelmet,'fhresetview');
        fhresetview(Helmethandles);
    end
end
handles.newlist=0;
end
guidata(handles.figure1,handles);

% --- Executes on selection change in popupmenu_zone.
function popupmenu_zone_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_zone
% contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_zone
global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
if  ~isempty(PMI{currentsub}.zone.label);
    if ~get(handles.radio_multiplezone,'value');
        num = get(handles.popupmenu_zone,'value');
        PMI{currentsub}.plotLst = PMI{currentsub}.zone.plotLst{num};
        list=PMI{currentsub}.plotLst;
        PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(list,1),PMI{currentsub}.data(cf).MeasList(list,2)];;
    else
         PMI{currentsub}.plotLst
         num = get(handles.popupmenu_zone,'value');
         
        PMI{currentsub}.plotLst = [reshape(PMI{currentsub}.plotLst,numel(PMI{currentsub}.plotLst),1);reshape(PMI{currentsub}.zone.plotLst{num},numel(PMI{currentsub}.zone.plotLst{num}),1)];
        list=PMI{currentsub}.plotLst;
        PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(list,1),PMI{currentsub}.data(cf).MeasList(list,2)];;
    end
end
set(guiHOMER,'UserData',PMI);
handles.newlist=1;
guidata(handles.figure1,handles);
updatedisplay(handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_zone_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_CreateZone.
function btn_CreateZone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_CreateZone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global currentsub
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf  = 1;
d=PMI{currentsub}.data(cf).HRF.AvgC;


time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
currentsub=1;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
if isempty(time_start); indstart = 1; else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart  =1;
    end
    indstart = indstart(end);
end

if isempty(time_stop); indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
    indstop = indstop(end);
end

MeasList = PMI{currentsub}.data(cf).MeasList;
MeasListAct = PMI{currentsub}.data(cf).MeasListAct;
intensnorm=d(indstart:indstop,:);
corrcoef = str2num(get(handles.edit_correlation,'string'));


usewv = get(handles.popup_corrwavelenght,'value');
zone = correlation_channel_gen(intensnorm,MeasListAct,MeasList,corrcoef,usewv);
if ~isempty(zone.color)
    PMI{currentsub}.zone = zone;
    izonelist = 1:numel(PMI{currentsub}.zone.plotLst);
    if get(handles.radio_colorzone,'value')==1
        color = zeros(size(handles.NIRS.Cf.H.C.id,2),3);
        for izone=izonelist;
            plotLst= PMI{currentsub}.zone.plotLst{izone};
            for izoneplotLst=1:numel(plotLst);
                color(plotLst(izoneplotLst),:) = PMI{currentsub}.zone.color(izone,:);
            end
        end
        PMI{1}.color = color;
    else
        PMI{1}.color=PMIcolordef(size(handles.NIRS.Cf.H.C.id,2)); %affichage de couleur des canaux = couleur des zones
    end
    set(handles.popupmenu_zone,'string',zone.label);
end

if 1 %numel(PMI{currentsub}.plotLst)==1
    xref = mean(intensnorm(:,PMI{currentsub}.plotLst),2);
    for ich=1:size(intensnorm,2);
        d2ok = intensnorm(:,ich);
        rall(ich) = corr(xref,d2ok);
    end
    plotLst = find(abs(rall)>corrcoef);
    PMI{currentsub}.plotLst = plotLst;
    set(handles.edit_plotLst,'string',mat2str(plotLst))
    
end
%Be sure both wavelength selected
listgood = PMI{currentsub}.plotLst ; % find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));
nc830 = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
all830 = zeros(numel(PMI{currentsub}.data(cf).MeasListAct),1);
all830(listgood)=1;
tmp = reshape(all830,nc830 ,2);
t = sum(tmp');
listgood = find([t,t]);
PMI{currentsub}.plotLst =listgood;
PMI{currentsub}.plot = [MeasList(plotLst,1),MeasList(plotLst,2)];


set(handles.popupmenu_zone,'value',1);
set(guiHOMER,'UserData',PMI);
handles.newlist=1; %Update 3dhelmet channel display
updatedisplay(handles);


function edit_time_start_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time_start as text
%        str2double(get(hObject,'String')) returns contents of edit_time_start as a double
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
str2num(get(handles.edit_time_start,'string'));
global currentsub;
PMI{currentsub}.tstart = str2num(get(handles.edit_time_start,'string'));
set(guiHOMER,'UserData',PMI);

updatedisplay(handles)

% --- Executes during object creation, after setting all properties.
function edit_time_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_view.
function popup_view_Callback(hObject, eventdata, handles)
% hObject    handle to popup_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_view contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_view

%Brushing data
if get(handles.popup_view,'value')==1|get(handles.popup_view,'value')==3
    set(handles.radio_brushingselectdata,'enable','on')
else
    set(handles.radio_brushingselectdata,'enable','off')
end
%Step option
if get(handles.popup_view,'value')==2|get(handles.popup_view,'value')==4|get(handles.popup_view,'value')==5
    set(handles.edit_stepvalue,'enable','on')
else
    set(handles.edit_stepvalue,'enable','on')
end
%Zone list
if get(handles.popup_view,'value')==5|get(handles.popup_view,'value')==2
    set(handles.edit_zonelist,'enable','on')
else
    set(handles.edit_zonelist,'enable','off')
end



updatedisplay(handles)

% --- Executes during object creation, after setting all properties.
function popup_view_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function display_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to display_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popup_average.
function popup_average_Callback(hObject, eventdata, handles)
% hObject    handle to popup_average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_average contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_average
updatedata(handles,1,1)
updatedisplay(handles)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_average_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --------------------------------------------------------------------
function context_Stop_Callback(hObject, eventdata, handles)
% hObject    handle to context_Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = get(handles.display_axes,'CurrentPoint');
posx =pos(1,1);
axes(handles.display_axes);
startini = str2num(get(handles.edit_time_start,'string'));
if isempty(startini)
    set(handles.edit_time_stop,'string', posx);
else
    if posx < startini
        set(handles.edit_time_start,'string',posx);
        set(handles.edit_time_stop,'string',posx);
    elseif posx > startini
        set(handles.edit_time_start,'string',startini);
        set(handles.edit_time_stop,'string',posx);
    end
end

axes(handles.display_axes);
set(handles.edit_time_stop,'string', posx);
guidata(hObject, handles);
updatedisplay(handles)

% --------------------------------------------------------------------
function context_stop_Callback(~, eventdata, handles)
% hObject    handle to context_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Option_Callback(hObject, eventdata, handles)
% hObject    handle to Option (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_time_stop_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time_stop as text
%        str2double(get(hObject,'String')) returns contents of edit_time_stop as a double
updatedisplay(handles)

% --- Executes during object creation, after setting all properties.
function edit_time_stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function context_newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to context_newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAxes_d1(handles,1)


% --------------------------------------------------------------------
function context_Start_Callback(hObject, eventdata, handles)
% hObject    handle to context_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = get(handles.display_axes,'CurrentPoint');
posx =pos(1,1);
axes(handles.display_axes);
stopini = str2num(get(handles.edit_time_stop,'string'));
if isempty(stopini)
    set(handles.edit_time_start,'string', posx);
else
    if posx < stopini
        set(handles.edit_time_start,'string',posx);
        set(handles.edit_time_stop,'string',stopini);
    elseif posx > stopini
        set(handles.edit_time_start,'string',posx);
        set(handles.edit_time_stop,'string',posx);
    end
end

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
str2num(get(handles.edit_time_start,'string'));
global currentsub;
PMI{currentsub}.tstart = str2num(get(handles.edit_time_start,'string'));
set(guiHOMER,'UserData',PMI);


guidata(hObject, handles);
updatedisplay(handles)



function edit_correlation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_correlation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_correlation as text
%        str2double(get(hObject,'String')) returns contents of edit_correlation as a double


% --- Executes during object creation, after setting all properties.
function edit_correlation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_correlation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updatedata(handles,loadnir,loadvmrk)

global currentsub
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf  = 1;
idmodule = get(handles.popupmenu_module, 'value');
idfile = get(handles.popupmenu_file,'value');
idmodulehold = get(handles.popupmenu_module_hold, 'value');
NC = handles.NIRS.Cf.H.C.N;
oldmodule = handles.oldmodule;
% PMI{currentsub}.data(cf).MeasListAct
% handles.NIRS.Cf.H.C.ok(:, oldmodule)=PMI{currentsub}.data(cf).MeasListAct;
% NIRS = handles.NIRS;
% save(handles.NIRSpath{1},'NIRS');
%PMI{currentsub}.data(cf).MeasListAct
%  NIRS = load(handles.NIRSpath{1});
%  handles.NIRS = NIRS.NIRS;
%Verification pour voir si le manual setting à été sauvegarder
if 0 %numel(handles.NIRS.Dt.fir.pp)==oldmodule
    oldfile = handles.oldfile;
    xname = handles.NIRS.Dt.fir.pp(oldmodule).p{oldfile};
    noise = logical(zeros(size(PMI{currentsub}.data(cf).HRF.AvgC)));
    mrk_type_arr = cellstr('bad_step');
    mrks = [];
    ind = [];
    dur = [];
    x = handles.NIRS.Dt.fir.pp(oldmodule).p{oldfile};
    vmrk_path = [x(1:end-3),'vmrk'];
    %         handles.file_vmrk = handles.NIRS.Dt.fir.pp(end).p{idfile}; %used
    %         in save noise
    
    [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
    if ~isempty(ind_dur_ch)
        maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
        badind = find(maxpoint>size(noise,1));
        if ~isempty(badind)
            disp(['Warning file ' vmrk_path ' marker : ' num2str(badind') ' are out of range'])
            ind_dur_ch(badind,2)=size(noise,1)- ind_dur_ch(badind,1);
        end
        for Idx = 1:size(noise,2)
            mrks = find(ind_dur_ch(:,3)==Idx);
            ind = ind_dur_ch(mrks,1);
            indf = ind + ind_dur_ch(mrks,2) - 1;
            if ~isempty(ind)
                try
                    for i = 1:numel(ind)
                        noise(ind(i):indf(i),Idx) = 1;
                    end
                catch
                    msgbox('Noise reading problem')
                end
            end
        end
    end
    if isfield(PMI{currentsub}.data(cf).HRF,'noise')
        if size(PMI{currentsub}.data(cf).HRF.noise,2)~=size(noise,2)
            msgbox('DEBUG NOISE DIMENSION MISMATCH')
        end
        if ~isempty(find(PMI{currentsub}.data(cf).HRF.noise~= noise))
            choice = questdlg(['Have you forget to save the noise before changing the data ? '],'Warning','Yes','No','Yes');
            if strcmp(choice,'Yes')
                set(handles.popupmenu_module, 'value', oldmodule);
                set(handles.popupmenu_file,'value',oldfile);
                guidata(gco, handles);
                return
            end
        end
    end
    
end %fin vérification
try
    if loadnir==1 %Standard
        %LOAD DATA
        % NIRS.Dt.AUX.pp(2).p{ifile,1}=outfileAUX;
        
        epochavg = 0;
        if numel(handles.NIRS.Dt.fir.pp(idmodule).pre)>15
            if  strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre(1:15), 'Epoch averaging');
                epochavg = 1;
            end
        end
        if epochavg
            [pathstr, name, ext] = fileparts(handles.NIRS.Dt.fir.pp(idmodule).p{idfile});
               xnumstim = fullfile(pathstr,[name,'_events',ext]);
                numstim = fopen_NIR(xnumstim,NC)';
                PMI{currentsub}.data(cf).HRF.navg = numstim;            
            if get(handles.popup_average,'value')==1 %mean
                xname =fullfile(pathstr,[name,ext]);             
            elseif get(handles.popup_average,'value')==2 %mean et std
                xname =fullfile(pathstr,[name,ext]);
                xstd = fullfile(pathstr,[name,'_std',ext]);;
                AvgStd = fopen_NIR(xstd,NC)';
                PMI{currentsub}.data(cf).HRF.AvgStdErr = AvgStd./ numstim.^0.5;    
            elseif get(handles.popup_average,'value')==3 %tval
                xname =fullfile(pathstr,[name,'_tval',ext]);
            elseif get(handles.popup_average,'value')==4 %event
                xname =fullfile(pathstr,[name,'_events',ext]);
            elseif 5 %tval & conc
                xname =fullfile(pathstr,[name,ext]);
                fopen_NIR(xname,NC)';
                PMI{currentsub}.data(cf).HRF = 2;
            end
        else
            if strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre,'READ_RAW_NIRS')
                xname = handles.NIRS.Dt.fir.pp(idmodule).p{idfile};
            else
                xname = handles.NIRS.Dt.fir.pp(idmodule).p{idfile};
            end
        end
        %  xname =fullfile(pathstr,[name,'_events',ext]);
        %         figure;plot(PMI{currentsub}.data(cf).HRF.AvgC)
        if get(handles.popup_typeDC,'value')==1     %DC
            PMI{currentsub}.data(cf).HRF.AvgC = fopen_NIR(xname,NC)';
        elseif get(handles.popup_typeDC,'value')==2 %PH
            xname = [xname(1:end-4),'ph.nir'];
            PMI{currentsub}.data(cf).HRF.AvgC = fopen_NIR([xname],NC)';
        elseif get(handles.popup_typeDC,'value')==3 %AC
            xname = [xname(1:end-4),'ac.nir']
            PMI{currentsub}.data(cf).HRF.AvgC = fopen_NIR(xname,NC)';
        else
            
        end
        PMI{currentsub}.data(cf).HRF.tHRF = 1/handles.NIRS.Cf.dev.fs:1/handles.NIRS.Cf.dev.fs:size(PMI{currentsub}.data(cf).HRF.AvgC,1)*1/handles.NIRS.Cf.dev.fs;
        %Timing setting for epoch average
        if epochavg
            PMI{currentsub}.data(cf).HRF.tHRF =PMI{currentsub}.data(cf).HRF.tHRF - str2num(handles.NIRS.Dt.fir.pp(idmodule).job.choiceave.pretime)  ;
            %                 figure
            %                 xname =fullfile(pathstr,[name,'_events',ext]);
            %                 events = fopen_NIR(xname,NC)';
            %                 bar(median(events)/max(events(:)));
            %                 title('events')
        end
        
        fileid=get(handles.popupmenu_file,'value');
        
        if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); %
            if isfield(handles.NIRS.Cf.H.C,'okavg')
                PMI{currentsub}.data(cf).MeasListAct = handles.NIRS.Cf.H.C.okavg;
            else
                PMI{currentsub}.data(cf).MeasListAct = ones(size (handles.NIRS.Cf.H.C.ok,1),1);
                disp('channellist reinitialized')
            end
        else
            PMI{currentsub}.data(cf).MeasListAct = handles.NIRS.Cf.H.C.ok(:, fileid);
        end
        %PMI{currentsub}.data(cf).MeasListAct(109)
        %Find if a normalisation module in the data.
        
        
        imodulenorm = [];
        for i=1:numel(handles.NIRS.Dt.fir.pp)
            if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
                if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
                    imodulenorm  = i;        %Le dernier si plusieurs Step Detection
                end
            end
        end
        
        if ~isempty(imodulenorm)
            if isfield(handles.NIRS.Dt.fir.pp(imodulenorm),'normlist')
                NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
                if isempty(NormList)
                    set(handles.listbox_normlist,'string','ND');
                    set(handles.listbox_normlist,'value',1);
                else
                    for i = 1:size(NormList,1)
                        if NormList(i,5)==1
                            state = '_On';
                        else
                            state = '_Off';
                        end
                        nblist = sprintf('%04.0f',i);
                        chname =sprintf('%03.0f',NormList(i,1));
                        valname =sprintf('%03.0f',NormList(i,4)*100);
                        Label = [nblist,'_',valname,'Ch',chname,state];
                        NormListLabel{i}=Label;
                    end
                    set(handles.listbox_normlist,'string',NormListLabel);
                    set(handles.listbox_normlist,'value',1);
                end
            end
            set(handles.uipanel_manuelNorm,'title',['Step detection ',num2str(size(NormList,1))]);
            
            
            
        end
    end
    if isfield(handles.NIRS.Dt,'AUX') %LOAD AUX ICI AUXILIARY AND ADD IN PMI structure
        try
            iAUX = get(handles.popupauxnb,'value');
            PMI{currentsub}.data(cf).AUX = [];
            nameAUX = handles.NIRS.Dt.AUX(iAUX).pp(end).p{idfile};
            if isfield(handles.NIRS.Dt.AUX(iAUX).pp(end),'sync_timesec')
                tstart = handles.NIRS.Dt.AUX(iAUX).pp(end).sync_timesec{idfile};
            else
                tstart = 1;
            end
            if tstart <0
               %msgbox('NEED TO PAD HERE ERROR')            
                tstop = PMI{currentsub}.data(cf).HRF.tHRF(end);
                 [data,infoBV,label,ind_dur_ch] = fopen_EEG(nameAUX, tstart, tstop);
        
            else 
                tstop = tstart+PMI{currentsub}.data(cf).HRF.tHRF(end);
                  [data,infoBV,label,ind_dur_ch] = fopen_EEG(nameAUX, tstart, tstop);
            end
          
            
            
            for ich=1:numel(infoBV.name_ele)
                PMI{currentsub}.data(cf).AUX.label{ich} =[infoBV.name_ele{ich},'  on'];
                PMI{currentsub}.data(cf).AUX.data{ich} = data(:,ich);
                PMI{currentsub}.data(cf).AUX.fs{ich} =1/(infoBV.SamplingInterval/1000000); %Frequence echantillonage Hz
                PMI{currentsub}.data(cf).AUX.view{ich} = 1;
            end
            set(handles.listbox_AUXCH,'string',PMI{currentsub}.data(cf).AUX.label')
            set(handles.listbox_AUXCH,'enable','on')
        catch
            try
                disp(['AUX ',  handles.NIRS.Dt.AUX(iAUX).pp(end).p{idfile},' could not be loaded']);
            catch
                disp(['AUX bloc', sprintf('%03.0f',idfile),' do not exist, tip: concatenate then segment according to trigger']);
            end
        end
    end
    
    
    if loadnir==2 %Hold
        %LOAD DATA
        PMI{2}= PMI{1};
        currentsubhold = 2; %reserver aux data hold
        epochavg = 0;
        if numel(handles.NIRS.Dt.fir.pp(idmodulehold).pre)>15
            if strcmp(handles.NIRS.Dt.fir.pp(idmodulehold).pre(1:15), 'Epoch averaging');
                epochavg = 1;
            end
        end
        
        if epochavg
            [pathstr, name, ext] = fileparts(handles.NIRS.Dt.fir.pp(idmodulehold).p{idfile})
            if get(handles.popup_average,'value')==1 %mean
                xname =fullfile(pathstr,[name,ext]);
            elseif get(handles.popup_average,'value')==2 %mean et std
                xname =fullfile(pathstr,[name,ext]);
                xstd = fullfile(pathstr,[name,'_std',ext]);
                xnumstim = fullfile(pathstr,[name,'_events',ext]);
                AvgStd = fopen_NIR(xstd,NC)';
                numstim = fopen_NIR(xnumstim,NC)';
                PMI{currentsubhold}.data(cf).HRF.AvgStdErr = AvgStd./ numstim.^0.5;
            elseif get(handles.popup_average,'value')==3 %tval
                xname =fullfile(pathstr,[name,'_tval',ext]);
            elseif get(handles.popup_average,'value')==4 %event
                xname =fullfile(pathstr,[name,'_events',ext]);
            elseif 5 %tval & conc
                xname =fullfile(pathstr,[name,ext]);
                fopen_NIR(xname,NC)';
                PMI{currentsubhold}.data(cf).HRF = 2;
            end
        else
            xname = handles.NIRS.Dt.fir.pp(idmodulehold).p{idfile};
        end
        
        
        %  xname =fullfile(pathstr,[name,'_events',ext]);
        %         figure;plot(PMI{currentsubhold}.data(cf).HRF.AvgC)
        if get(handles.popup_typeDC,'value')==1     %DC
            PMI{currentsubhold}.data(cf).HRF.AvgC = fopen_NIR(xname,NC)';
        elseif get(handles.popup_typeDC,'value')==2 %PH
            xname = [xname(1:end-4),'ph.nir'];
            PMI{currentsubhold}.data(cf).HRF.AvgC = fopen_NIR([xname],NC)'
        elseif get(handles.popup_typeDC,'value')==3 %AC
            xname = [xname(1:end-4),'ac.nir']
            PMI{currentsubhold}.data(cf).HRF.AvgC = fopen_NIR(xname,NC)';
        else
            
        end
        PMI{currentsubhold}.data(cf).HRF.tHRF = 1/handles.NIRS.Cf.dev.fs:1/handles.NIRS.Cf.dev.fs:size(PMI{currentsubhold}.data(cf).HRF.AvgC,1)*1/handles.NIRS.Cf.dev.fs;
        %Timing setting for epoch average
        if epochavg
            PMI{currentsubhold}.data(cf).HRF.tHRF =PMI{currentsub}.data(cf).HRF.tHRF - str2num(handles.NIRS.Dt.fir.pp(idmodulehold).job.choiceave.pretime)  ;
            %                 figure
            %                 xname =fullfile(pathstr,[name,'_events',ext]);
            %                 events = fopen_NIR(xname,NC)';
            %                 bar(median(events)/max(events(:)));
            %                 title('events')
        end
        
        fileid=get(handles.popupmenu_file,'value');
        
%         if epochavg
%             PMI{currentsubhold}.data(cf).MeasListAct = handles.NIRS.Dt.fir.pp(idmodulehold).chok;
%         else
%             PMI{currentsubhold}.data(cf).MeasListAct = handles.NIRS.Cf.H.C.ok(:, fileid);
%         end
        
        %Find if a normalisation module in the data.
        
        
        imodulenorm = [];
        for i=1:numel(handles.NIRS.Dt.fir.pp)
            if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
                if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
                    imodulenorm  = i;        %Le dernier is plusieur Step Detection
                end
            end
        end
        
        if ~isempty(imodulenorm)
            if isfield(handles.NIRS.Dt.fir.pp(imodulenorm),'normlist')
                NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
                if isempty(NormList)
                    set(handles.listbox_normlist,'string','ND');
                    set(handles.listbox_normlist,'value',1);
                else
                    for i = 1:size(NormList,1)
                        if NormList(i,5)==1
                            state = '_On';
                        else
                            state = '_Off';
                        end
                        nblist = sprintf('%04.0f',i);
                        chname =sprintf('%03.0f',NormList(i,1));
                        valname =sprintf('%03.0f',NormList(i,4)*100);
                        Label = [nblist,'_',valname,'Ch',chname,state];
                        NormListLabel{i}=Label;
                    end
                    set(handles.listbox_normlist,'string',NormListLabel);
                    set(handles.listbox_normlist,'value',1);
                end
            end
            set(handles.uipanel_manuelNorm,'title',['Step detection ',num2str(size(NormList,1))]);
        end
    end
    PMI{1}.NIRSmat=handles.NIRSall;
    try
        PMI{1}.prj_name=handles.NIRS.Cf.H.prj;
    catch
    end
    idfile=get(handles.popupmenu_file,'value');
    PMI{1}.idfile=idfile ;
catch
    h=msgbox(['Data ',xname,' can''t be loaded']);
    uiwait(h)
    return
end
set(guiHOMER,'UserData',PMI);

try
    if loadvmrk
        %Load VMRK
        noise = logical(zeros(size(PMI{currentsub}.data(cf).HRF.AvgC)));
        mrk_type_arr = cellstr('bad_step');
        mrks = [];
        ind = [];
        dur = [];
        
        x = handles.NIRS.Dt.fir.pp(idmodule).p{idfile};
        vmrk_path = [x(1:end-3),'vmrk'];
        %         handles.file_vmrk = handles.NIRS.Dt.fir.pp(end).p{idfile}; %used
        %         in save noise
        [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
        if ~isempty(ind_dur_ch)
            maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
            badind = find(maxpoint>size(noise,1));
            if ~isempty(badind)
                disp(['Warning file ' vmrk_path ' marker : ' num2str(badind') ' are out of range in the data file'])
                ind_dur_ch(badind,2)=size(noise,1)- ind_dur_ch(badind,1);
            end
            for Idx = 1:size(noise,2)
                mrks = find(ind_dur_ch(:,3)==Idx);
                ind = ind_dur_ch(mrks,1);
                indf = ind + ind_dur_ch(mrks,2) - 1;
                if ~isempty(ind)
                    try
                        for i = 1:numel(ind)
                            noise(ind(i):indf(i),Idx) = 1;
                        end
                    catch
                        msgbox('Noise reading problem')
                    end
                end
            end
        end
        PMI{currentsub}.data(cf).HRF.noise  = noise;
        
    end
    set(guiHOMER,'UserData',PMI);
catch
    h=msgbox(['Data ',handles.file_vmrk,' can''t be loaded']);
    uiwait(h)
    set(guiHOMER,'UserData',PMI);
end






function edit_start_Callback(hObject, eventdata, handles)
% hObject    handle to edit_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_start as text
%        str2double(get(hObject,'String')) returns contents of edit_start as a double
updatedisplay(handles)

% --- Executes during object creation, after setting all properties.
function edit_start_CreateFcn(hObject, eventdata, ~)
% hObject    handle to edit_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_duration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_duration as text
%        str2double(get(hObject,'String')) returns contents of edit_duration as a double
updatedisplay(handles)

% --- Executes during object creation, after setting all properties.
function edit_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_plus.
function btn_plus_Callback(hObject, eventdata, handles)
% hObject    handle to btn_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

start = str2num(get(handles.edit_start,'string'));
step = str2num(get(handles.edit_step,'string'));
start = start+step;
set(handles.edit_start,'string',num2str(start));
updatedisplay(handles);

% --- Executes on button press in btn_moins.
function btn_moins_Callback(hObject, eventdata, handles)
% hObject    handle to btn_moins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
start = str2num(get(handles.edit_start,'string'));
step = str2num(get(handles.edit_step,'string'));
start = start-step;
set(handles.edit_start,'string',num2str(start));
updatedisplay(handles);


function edit_step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_step as text
%        str2double(get(hObject,'String')) returns contents of edit_step as a double


% --- Executes during object creation, after setting all properties.
function edit_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function context_restore_Callback(hObject, eventdata, handles,state)
% hObject    handle to context_restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

name=get(gco,'Displayname');
[tok,rem]=strtok(name,'ch');
izone = str2num(tok(5:end));
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
NChalf = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
plotLst = PMI{currentsub}.zone.plotLst{izone};
plotLstAll = [];
for i=1:numel(plotLst)
    chlist = plotLst(i);
    if isempty(find(chlist - NChalf < 0 ))
        plotLstAll = [plotLstAll,chlist, chlist-NChalf]; %add 830 nm channels
    else
        plotLstAll = [plotLstAll,chlist,chlist+NChalf]; %add 690 nm channels
    end
end
plotLst = plotLstAll;

posxstart = str2num(get(handles.edit_time_start,'string'));
posxstop = str2num(get(handles.edit_time_stop,'string'));
indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
PMI{currentsub}.data(cf).HRF.noise(indstart(end):indstop(end),plotLst)= state;
set(guiHOMER,'UserData',PMI);
updatedisplay(handles);

% --------------------------------------------------------------------
function context_Artefact_Callback(hObject, eventdata, handles)
% hObject    handle to context_Artefact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_savenoise.
function handles=btn_savenoise_Callback(hObject, eventdata, handles)
% hObject    handle to btn_savenoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
fileid=get(handles.popupmenu_file,'value');
x = handles.NIRS.Dt.fir.pp(end).p{fileid};
vmrk_path = [x(1:end-3),'vmrk'];
[label,ind_dur_ch] = read_vmrk_all(vmrk_path);
if get(handles.popupmenu_module,'value')~=numel(handles.NIRS.Dt.fir.pp);
    if strcmp(get(handles.context_warning_noise,'checked'),'off')
        choice = questdlg(['You are not working in ''Manual module''',...
            'This will reset your previous manual selection with the current module. Do you want to continue ?'],'Warning','Yes','No','No');
        if strcmp(choice,'No')
            return
        end
    end
end
%remove old bad steps
list_bad_step =  [];
for i=1:size(label,1)
    if strmatch(label{i,1},'bad_step')
        list_bad_step = [list_bad_step; i];
    end
end
label(list_bad_step,:)=[];
ind_dur_ch(list_bad_step,:) = [];

%%% Fill the field .noise with "Data select" info
% try
plotLst = PMI{currentsub}.plotLst;
hBrushLine = findall(gca,'tag','Brushing');
brushedData = get(hBrushLine, {'Ydata'});
a = 1;
non_ch_brush = [];
for i = 1:length(brushedData)
    if length(brushedData{i}') ~= length(PMI{currentsub}.data(cf).HRF.noise(:,1))
        non_ch_brush = [non_ch_brush, i];
    end
end
brushedData(non_ch_brush) = []; %Get rid of non-channel-related cells
for i = 1:length(brushedData)
    brushedIdx = ~isnan(brushedData{i});
    while brushedData{i}(brushedIdx) ~= (PMI{currentsub}.data(cf).HRF.AvgC(brushedIdx,plotLst(a)))'
        a = a + 1;
        if a > numel(plotLst)
            a = a - 1;
            break
        end
    end
    if sum(brushedIdx) ~= 0
        PMI{currentsub}.data(cf).HRF.noise(:,plotLst(a)) = PMI{currentsub}.data(cf).HRF.noise(:,plotLst(a)) | brushedIdx';
    end
    a = 1;
end
PMI{currentsub}.data(cf).HRF.noise  = copy_channel_noise(PMI{currentsub}.data(cf).HRF.noise);
% end


%look new old bad step
ind_dur_ch_2 = mat2d2ind_dur_ch(PMI{currentsub}.data(cf).HRF.noise');
label_2 = cell(size(ind_dur_ch_2,1),2);
label_2(:,1)={'bad_step'};
label_2(:,2)={'manual'};
ind_dur_ch = [ind_dur_ch;ind_dur_ch_2];
label = [label;label_2];

write_vmrk_all(vmrk_path,ind_dur_ch,label);
Nc = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
MeasListAct = [double(sum(reshape(PMI{currentsub}.data(cf).MeasListAct,Nc,2),2)>1);double(sum(reshape(PMI{currentsub}.data(cf).MeasListAct,Nc,2),2)>1)];
PMI{currentsub}.data(cf).MeasListAct = MeasListAct;
  if ~isempty( strfind(handles.NIRS.Dt.fir.pp(end).pre, 'Epoch averaging')); %
    handles.NIRS.Cf.H.C.okavg(:, fileid)=MeasListAct;
  else      
    handles.NIRS.Cf.H.C.ok(:, fileid)=MeasListAct;
  end
NIRS = handles.NIRS;
save(handles.NIRSpath{1},'NIRS');
pause(0.5)
set(guiHOMER,'UserData',PMI);
%guidata(hObject,handles);
%updatedisplay(handles);




% --------------------------------------------------------------------
function context_restoreall_Callback(hObject, eventdata, handles,state)
% hObject    handle to context_restoreall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
NChalf = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
plotLst = PMI{currentsub}.plotLst;
plotLstAll = [];
for i=1:numel(plotLst)
    chlist = plotLst(i);
    if isempty(find(chlist - NChalf < 1 ))
        plotLstAll = [plotLstAll,chlist, chlist-NChalf]; %add 830 nm channels
    else
        plotLstAll = [plotLstAll,chlist,chlist+NChalf]; %add 690 nm channels
    end
end
plotLst = plotLstAll;
posxstart = str2num(get(handles.edit_time_start,'string'));
posxstop = str2num(get(handles.edit_time_stop,'string'));
indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
if isempty(indstart)
    indstart = 1;
end
indstart = indstart(end);
indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
if isempty(indstop)
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
end
indstop = indstop(end);
PMI{currentsub}.data(cf).HRF.noise(indstart:indstop,plotLst)= state;
set(guiHOMER,'UserData',PMI);
updatedisplay(handles);


% --- Executes on button press in btn_PCA.
function btn_PCA_Callback(hObject, eventdata, handles)
% hObject    handle to btn_PCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.radio_PCA,'value',1);


guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart = 1;
    end
end

if isempty(time_stop);
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
end
gui_PCAio(indstart(end),indstop(end));

% if isempty(time_start)|isempty(time_stop)
%     set(handles.text_Advice,'string','Please enter time start and time stop to define time window')
%     return
% end
%
% tstart=find(PMI{currentsub}.data(cf).HRF.tHRF<time_start);
% if isempty(tstart);tstart=1;end
% tstop=find(PMI{currentsub}.data(cf).HRF.tHRF<time_stop);
% indt = tstart(end):tstop(end);
% intensnorm = d(indt,:);
%
% %FILTER
% AdvOptions.FiltOptions.FilterType=1;  %Default Butterworth
% AdvOptions.FiltOptions.FilterOrder=3;  %Default 3rd order
% lpf =str2num(get(handles.edit_LPF,'string'));
% fs = handles.NIRS.Cf.dev.fs;
% [fb,fa] = MakeFilter(AdvOptions.FiltOptions.FilterType,AdvOptions.FiltOptions.FilterOrder,fs,lpf,'low');
% intensnorm = filtfilt(fb,fa,intensnorm);
%
% MeasList = PMI{currentsub}.data(cf).MeasList;
% MeasListAct = PMI{currentsub}.data(cf).MeasListAct;
% coefcorr = str2num(get(handles.edit_correlation,'string'));
% usewv = get(handles.popup_corrwavelenght,'value');
% zone = correlation_channel_gen(intensnorm,MeasListAct,MeasList,coefcorr,usewv);
%
% t = PMI{currentsub}.data(cf).HRF.tHRF;
% [h,dataPCA,measlistact]  = Peak_Marking_PCA(intensnorm,zone,t(indt),MeasListAct );    %d1ok segment of normalise intensity
% % figure
% % plot(dataPCA)
% dnew = d;
% dnew(indt,:) = d(indt,:) - dataPCA;
% list = find(sum(dataPCA)~=0)
%
% %FIT END remove step with the rest of the data at the begining and the end
% %of the data
% nb = numel(t)
% if indt(end)<nb
%     delta = d(indt(end)+1,list)-dnew(indt(end),list);
%     for i=1:numel(list)
%         dnew(indt(end):nb,list(i)) = d(indt(end):nb,list(i))-delta(i);
%     end
% end
% if indt(1)>2
%     delta = d(indt(1),list)-dnew(indt(1)-1,list);
%     for i=1:numel(list)
%         dnew(1:indt(1),list(i)) = d(1:indt(1),list(i))-delta(i);
%     end
% end
%
% figure;
% hold on;
% step=str2num(get(handles.edit_stepvalue,'string'));
% offset = 0;
% if (indt(1)-10) < 1
%     indt(1)=1;
% else
%     indt(end)=indt(end)-10;
% end
% if (indt(end)+10) > numel(t)
%     indt(end)=numel(t);
% else
%     indt(end)=indt(end)+10;
% end
% indtlarge = indt(1):indt(end);
% set(gca,'fontsize',14);
% for ilist = 1:numel(list)
%     h= plot(t(indtlarge),d(indtlarge,list(ilist))+offset,'b');
% 	h= plot(t(indtlarge),dnew(indtlarge,list(ilist))+offset,'r');
%     set(h,'displayname',['CH',num2str( list(ilist))]);
%     offset = offset-step;
% end
% xlim([t(indtlarge(1)),t(indtlarge(end))]);
% title('Modif suggested in red','fontsize',14)
% set(gcf,'unit','normalized','position',[0,0,0.3,1])
% PMI{currentsub}.data(cf).HRF.AvgC = dnew;
% handles.newlist=1; %update helmet channel
% PMI{currentsub}.plotLst =  list;
% PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(list,1),PMI{currentsub}.data(cf).MeasList(list,2)];
% set(guiHOMER,'UserData',PMI);
updatedisplay(handles);

% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

trig_ind = get(handles.popupmenu8, 'value');
trig_name = get(handles.popupmenu8, 'string');
if isfield(handles,'trigger')
    handles.selected_trig = handles.triggers(ismember(handles.triggers(:,1), str2num(char(trig_name{trig_ind}))),:);
else
    
end
guidata(hObject,handles);
updatedisplay(handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_brushingselectdata.
function radio_brushingselectdata_Callback(hObject, eventdata, handles)
% hObject    handle to radio_brushingselectdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
brush_on = get(hObject,'Value');
axes(handles.display_axes);
if brush_on
    brush on
else
    brush off
end
% Pour aller chercher les variables sélectionnées avec le brush. Pourrait
% être utile de relier cette fonction à celles qui nécessitent des time
% start/stop.
% handles.edit_time_start & handles.edit_time_stop???

% Hint: get(hObject,'Value') returns toggle state of radio_brushingselectdata


% --- Executes on button press in radio_showrejectedchannels.
function radio_showrejectedchannels_Callback(hObject, eventdata, handles)
% hObject    handle to radio_showrejectedchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updatedisplay(handles);

%Remove rejected channels from this plotting list
%Create a second plotting list with only bad channels present in the first
%   -> Plot


% Hint: get(hObject,'Value') returns toggle state of radio_showrejectedchannels


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

choice = questdlg(['Are you sure you want to modify the data ?'],'Be carefull', ...
    'Yes', ...
    'No',...
    'No');
if strcmp(choice,'No')
    return
end
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
plotLst = PMI{currentsub}.plotLst;
d1 = PMI{currentsub}.data(cf).HRF.AvgC;
%     Save new d
if 0 %enlever vieille version Hubert à vérifier l'interpollation se fera avec le bruit dans la toolbox
    h = handles.display_axes;
    lineH = findobj(h,'type','line','-and','-not','Marker','*'); %Get only channels data, not markers
    newData = get(lineH, {'xdata','ydata'});
    plotlist_tag = get(lineH, 'tag');
    for i = 1:numel(plotlist_tag)
        if iscell(plotlist_tag)
            plot_ch = str2double(plotlist_tag{i});
        else
            plot_ch = str2double(plotlist_tag(i));
        end
        if size(newData{i,1}) == size(PMI{currentsub}.data(cf).HRF.tHRF) %Look if points were removed
            d1(:,plot_ch) = newData{i,2};
        else %Interpolate
            time_to_interp = ~ismember(PMI{currentsub}.data(cf).HRF.tHRF, newData{i,1});
            newData{i,2}(~time_to_interp) = newData{i,2}; %Make room for interpolation
            newData{i,2}(time_to_interp) = NaN;
            %         ind_interp = find(time_to_interp);
            %         diff_interp = diff(ind_interp);
            %         for j = 1:numel(diff_interp)
            %
            %         end
            %         x1 = ind_interp(1)-1;
            %         x2 = ind_interp(end)+1;
            %         a = (d1(x2,i)-d1(x1,i))/(x2-x1);
            %         b = d1(x1,i) - a*x1;
            %         newData{i,2}(time_to_interp) = a.*(x1+1:x2-1) + b;
            %Save data in d1
            d1(:,plot_ch) = newData{i,2};
        end
    end
end

%Save data in PMI structure
PMI{currentsub}.data(cf).HRF.AvgC  = d1;

idfile = get(handles.popupmenu_file,'value');
outfile = handles.NIRS.Dt.fir.pp(end).p{idfile};
fwrite_NIR(outfile,d1');

set(guiHOMER,'UserData',PMI);
updatedisplay(handles);


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8


% --------------------------------------------------------------------
function EditMenu_Callback(hObject, eventdata, handles)
% hObject    handle to EditMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function modNIRS_Callback(hObject, eventdata, handles)
% hObject    handle to modNIRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.NIRS] = NIRS_GUI(handles.NIRSpath{handles.subjectnb,1});
[handles] = set_popupmenus(handles);
guidata(hObject, handles);


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
fileid=get(handles.popupmenu_file,'value');
x = handles.NIRS.Dt.fir.pp(end).p{fileid};
vmrk_path = [x(1:end-3),'vmrk'];
[label,ind_dur_ch] = read_vmrk_all(vmrk_path);
%remove old bad steps
list_bad_step =  [];
for i=1:size(label,1)
    if strmatch(label{i,1},'bad_step')
        list_bad_step = [list_bad_step; i];
    end
end
label(list_bad_step,:)=[];
ind_dur_ch(list_bad_step,:) = [];

%%% Fill the field .noise with "Data select" info
% try
plotLst = PMI{currentsub}.plotLst;
hBrushLine = findall(gca,'tag','Brushing');
brushedData = get(hBrushLine, {'Ydata'});
a = 1;
non_ch_brush = [];
for i = 1:length(brushedData)
    if length(brushedData{i}') ~= length(PMI{currentsub}.data(cf).HRF.noise(:,1))
        non_ch_brush = [non_ch_brush, i];
    end
end
brushedData(non_ch_brush) = []; %Get rid of non-channel-related cells
for i = 1:length(brushedData)
    brushedIdx = ~isnan(brushedData{i});
    while brushedData{i}(brushedIdx) ~= (PMI{currentsub}.data(cf).HRF.AvgC(brushedIdx,plotLst(a)))'
        a = a + 1;
        if a>numel(plotLst)
            a = a-1;
            break
        end
    end
    if sum(brushedIdx) ~= 0
        ch_noise = PMI{currentsub}.data(cf).HRF.noise(:,plotLst(a));
        ch_noise(brushedIdx) = 0;
        PMI{currentsub}.data(cf).HRF.noise(:,plotLst(a)) = ch_noise;
    end
    a = 1;
end
%Copy to the co-channel
noise = PMI{currentsub}.data(cf).HRF.noise;
NC = size(noise,2);
noise(:,1:NC/2) = noise(:,1:NC/2) & noise(:,NC/2+1:end);
noise(:,NC/2+1:end) = noise(:,1:NC/2);
PMI{currentsub}.data(cf).HRF.noise = noise;
% end

%look new old bad step
ind_dur_ch_2 = mat2d2ind_dur_ch(PMI{currentsub}.data(cf).HRF.noise');
label_2 = cell(size(ind_dur_ch_2,1),2);
label_2(:,1)={'bad_step'};
label_2(:,2)={'manual'};
ind_dur_ch = [ind_dur_ch;ind_dur_ch_2];
label = [label;label_2];

write_vmrk_all(vmrk_path,ind_dur_ch,label)

set(guiHOMER,'UserData',PMI);
updatedisplay(handles);


% --------------------------------------------------------------------
function remove_entire_ch_Callback(hObject, eventdata, handles)
% hObject    handle to remove_entire_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
fileid=get(handles.popupmenu_file,'value');
x = handles.NIRS.Dt.fir.pp(end).p{fileid};

%Get the co-channel too
ch1 = str2double(get(gco, 'Tag'));
if ch1 <= handles.NIRS.Cf.H.C.N/2
    ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
else
    ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
end
%ici
idmodule = get(handles.popupmenu_module, 'value');

epochavg = 0;
if numel(handles.NIRS.Dt.fir.pp(idmodule).pre)>15
    if strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre(1:15), 'Epoch averaging');
        epochavg = 1;
    end
end


if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
    handles.NIRS.Cf.H.C.okavg([ch1 ch2],:) =0;
else
    handles.NIRS.Cf.H.C.ok([ch1 ch2], fileid) = 0;
end
NIRS=handles.NIRS;
save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');

PMI{currentsub}.data(cf).MeasListAct([ch1 ch2]) = 0;
set(guiHOMER,'UserData',PMI)
handles.newlist=1; %update helmet channel
guidata(hObject,handles);
updatedisplay(handles);



% --------------------------------------------------------------------
function remove_bad_ch_Callback(hObject, eventdata, handles)
% hObject    handle to remove_bad_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
namech = get(gco,'displayname');
set(handles.menu_selectch,'label',namech);




% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function restore_entire_ch_Callback(hObject, eventdata, handles)
% hObject    handle to restore_entire_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
fileid=get(handles.popupmenu_file,'value');
x = handles.NIRS.Dt.fir.pp(end).p{fileid};

%Get the co-channel too
ch1 = str2double(get(gco, 'Tag'));
if ch1 <= handles.NIRS.Cf.H.C.N/2
    ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
else
    ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
end
idmodule = get(handles.popupmenu_module, 'value');

epochavg = 0;
if numel(handles.NIRS.Dt.fir.pp(idmodule).pre)>15
    if strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre(1:15), 'Epoch averaging');
        epochavg = 1;
    end
end

if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
    handles.NIRS.Cf.H.C.okavg([ch1 ch2],:) =1;
else
    handles.NIRS.Cf.H.C.ok([ch1 ch2], fileid) = 1;
end
NIRS=handles.NIRS;
save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');

PMI{currentsub}.data(cf).MeasListAct([ch1 ch2])=1;
set(guiHOMER,'UserData',PMI);
handles.newlist=1; %update helmet channel
guidata(hObject,handles);
updatedisplay(handles);


% --- Executes on selection change in popup_typeDC.
function popup_typeDC_Callback(hObject, eventdata, handles)
% hObject    handle to popup_typeDC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_typeDC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_typeDC
updatedata(handles,1,1)
updatedisplay(handles)

% --- Executes during object creation, after setting all properties.
function popup_typeDC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_typeDC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_reject_ch.
function btn_reject_ch_Callback(hObject, eventdata, handles)
% hObject    handle to btn_reject_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
badch = (sum(PMI{currentsub}.data(cf).HRF.noise)./size(PMI{currentsub}.data(cf).HRF.noise,1)>0.80)
chbad = find(badch(1:end/2)| badch(end/2+1:end));
nbch = numel(badch(1:end/2));
fileid=get(handles.popupmenu_file,'value');

if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
    handles.NIRS.Cf.H.C.okavg([ch1 ch2],:) =0;
else
    handles.NIRS.Cf.H.C.ok([chbad chbad+nbch], fileid) = 0;
end
NIRS=handles.NIRS;
save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');

PMI{currentsub}.data(cf).MeasListAct([chbad chbad+nbch])=0;
set(guiHOMER,'UserData',PMI);


function edit_zonelist_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zonelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zonelist as text
%        str2double(get(hObject,'String')) returns contents of edit_zonelist as a double
handles.newlist=1; %Update 3dhelmet channel display
updatedisplay(handles);

% --- Executes during object creation, after setting all properties.
function edit_zonelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zonelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_TrialSelection.
function btn_TrialSelection_Callback(hObject, eventdata, handles)
% hObject    handle to btn_TrialSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


idmodule = get(handles.popupmenu_module, 'value');
idfile = get(handles.popupmenu_file,'value');
xname = handles.NIRS.Dt.fir.pp(idmodule).p{end};
try
    load([xname(1:end-3),'mat'])
catch %VIEW MULTIMODDULE
    pathout = uigetdir();
    if isnumeric(pathout)
        return
    end
    pathout = [pathout,'/'];
    namemodule = get(handles.popupmenu_module, 'string');
    moduleselected = [];
    for i=1:numel(namemodule)
        a = namemodule{i}
        if strcmp(a(end-2:end),' on')
            moduleselected=[ moduleselected,i];
        end
    end
    if isempty(moduleselected)
        return
    end
    htitle = figure
    set(gca,'visible','off')
    for imodule=1:numel(moduleselected)
        texttoadd = handles.NIRS.Dt.fir.pp(moduleselected(imodule)).pre
        cla
        text(0,0,[texttoadd],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
        saveas(htitle ,'temp','tif');
        imgall = imread('temp','tif');
        indy = find(sum(sum(imgall,3)~=765,1));
        indx = find(sum(sum(imgall,3)~=765,2));
        indx2=indx(1):indx(end);
        indy2=indy(1):indy(end);
        titre{imodule} = imgall(indx2,indy2,:);
    end
    close(htitle)
    
    for ibloc = 1:numel(get(handles.popupmenu_file,'string'))
        gcaall = [];
        for imodule = 1:numel(moduleselected)
            %%%POPUP FILE
            set(handles.popupmenu_file,'value',ibloc)
            set(handles.popupmenu_module,'value',moduleselected(imodule))
            lst = length(handles.NIRS.Dt.fir.pp); %get the number of processing step
            %Get all the processing steps (modules)
            handles.module_list = [];
            handles.module_path = [];
            
            trig_ind = get(handles.popupmenu8, 'value');
            idfile = get(handles.popupmenu_file,'value');
            if isfield(handles.NIRS.Dt.fir,'aux5')
                if numel(handles.NIRS.Dt.fir.aux5)>=idfile
                    handles.triggers = handles.NIRS.Dt.fir.aux5{idfile};
                    unique_trig = unique(handles.triggers(:,1));
                    trig_list = cellstr(num2str(unique_trig));
                    set(handles.popupmenu8, 'String', trig_list);
                    if numel(trig_list)< trig_ind
                        set(handles.popupmenu8, 'value', numel(trig_list));
                    end
                else
                    handles.triggers = [1,1];
                end
            else
                handles.triggers = [];
            end
            trig_ind = get(handles.popupmenu8, 'value');
            trig_name = get(handles.popupmenu8, 'string');
            try
                for itrig = 1:numel(trig_ind)
                    handles.selected_trig(itrig,:) = handles.triggers(ismember(handles.triggers(:,1), str2num(trig_name{trig_ind(itrig)})),:);
                end
            catch
            end
            updatedata(handles,1,1)
            handles.newlist = 1; %Reactualiser la liste de canaux de 3DHelmet
            updatedisplay(handles)
            handles.oldfile = get(handles.popupmenu_file, 'Value');
            guidata(hObject, handles);
            %FIN POPUP FILE
            pos = get(gca,'position');
            posfig =  get(gcf,'position');
            cap.x = posfig(1)+ pos(1)*posfig(3)-20;
            cap.y = posfig(2)+ pos(2)*posfig(4)-20;
            cap.widthSS =  pos(3)*posfig(3)+30;
            cap.heightSS = pos(4)*posfig(4)+30;
            figure(handles.figure1);
            gcaview = ScreenShotRGB(cap.x ,cap.y,cap.widthSS, cap.heightSS);
            
            gcaview(1:size(titre{imodule},1),1:size(titre{imodule},2),:) = titre{imodule};
            gcaall = [gcaall;gcaview];
            
            %         filename = ['module',sprintf('%03.0f',moduleselected(imodule))]
            %         saveas(gcf,filename,'jpg')
        end %FIN MODULE
        %filename = ['module',sprintf('%03.0f',moduleselected(imodule)),'Bloc', sprintf('%03.0f',ibloc),'.jpg']
        filename = [pathout,'Blocks', sprintf('%03.0f',ibloc),'.jpg'];
        imwrite(gcaall,filename,'jpg');
    end %FIN BLOC
    
    return
end
%open figure
hfigure(1) = NIRS_AvgTrial_Noise;%figure;
hold on
cla
title('Trial of the average')

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
idfile = get(handles.popupmenu_file,'value');
cf = PMI{currentsub}.currentFile;
plotLst = PMI{currentsub}.plotLst;




%Initialisation du NOISE lire dans les vmrk correspondants au fichier
%du InfoTrig selectionner dans manual selection les valeurs actuels du
%bruit pour chaque trial.
for itrial = 1:numel(InfoTrig)
    k = strfind(handles.NIRS.Dt.fir.pp(end).p(InfoTrig(itrial).filenb),  InfoTrig(itrial).file)
    if isempty(k)
        msgbox('The file does''nt match with the manual file please be sure you have selected the data and not the average for this step')
        return
    end
    filenblist = InfoTrig(itrial).filenb
end
for ifilenb =1:filenblist
    xname = handles.NIRS.Dt.fir.pp(end).p{ifilenb};
    vmrk_path = [xname(1:end-3),'vmrk'];
    mrk_type_arr = cellstr('bad_step');
    mrks = [];
    ind = [];
    dur = [];
    %Initialisation du Noise pour chaque fichier
    noise = logical(zeros(InfoTrig(itrial).nbpoint,size(A,1)));
    [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
    if ~isempty(ind_dur_ch)
        for Idx = 1:size(noise,2)
            mrks = find(ind_dur_ch(:,3)==Idx);
            ind = ind_dur_ch(mrks,1);
            indf = ind + ind_dur_ch(mrks,2) - 1;
            if ~isempty(ind)
                try
                    for i = 1:numel(ind)
                        noise(ind(i):indf(i),Idx) = 1;
                    end
                catch
                    msgbox('Noise reading problem')
                end
            end
        end
    end
    filenoise{ifilenb}.noise=noise;
    filenoise{ifilenb}.name=vmrk_path;
end




for itrial = 1:numel(InfoTrig)
    tbloc =  1/handles.NIRS.Cf.dev.fs:1/handles.NIRS.Cf.dev.fs:size(A,2)*1/handles.NIRS.Cf.dev.fs;
    %nbavg=(InfoTrig(itrial).stoppoint-InfoTrig(itrial).startpoint)
    for ind_ch = 1:numel(plotLst) %HbO %les canaux
        ich = plotLst(ind_ch);
        figure(hfigure(1));
        fig_handles = guihandles(hfigure(1));
        isnoise = find(noise(InfoTrig(itrial).startpoint:InfoTrig(itrial).stoppoint,ich));
        h = plot(tbloc,A(ich,:,itrial),'linewidth',2);
        if handles.NIRS.Cf.H.C.wl(ich)==2
            set(h,'linestyle','--');
        end
        set(h,'uicontextmenu',fig_handles.Noise)
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList,ich, numel(PMI{currentsub}.color)/3);
        set(h,'color',PMI{currentsub}.color(idxc,:));
        noise = filenoise{InfoTrig(itrial).filenb}.noise;
        set(h,'Displayname',['file',num2str(InfoTrig(itrial).filenb),'ch_',num2str(ich),'ind_',num2str(InfoTrig(itrial).startpoint),':',num2str(InfoTrig(itrial).stoppoint)]);
        plot(tbloc(isnoise),A(ich,isnoise,itrial),'xy');
    end
end
PMI{currentsub}.filenoise  = filenoise;
plot([tbloc(1),tbloc(end)],[0,0],'k')

set(guiHOMER,'UserData',PMI)



function edit_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min as text
%        str2double(get(hObject,'String')) returns contents of edit_min as a double
updatedisplay(handles)


% --- Executes during object creation, after setting all properties.
function edit_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max as text
%        str2double(get(hObject,'String')) returns contents of edit_max as a double
updatedisplay(handles)


% --- Executes during object creation, after setting all properties.
function edit_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_autoaxis.
function radio_autoaxis_Callback(hObject, eventdata, handles)
% hObject    handle to radio_autoaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_autoaxis
updatedisplay(handles)


% --------------------------------------------------------------------
function menu_selectch_Callback(hObject, eventdata, handles)
% hObject    handle to menu_selectch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
plotLst = str2double(get(gco, 'Tag'));
PMI{currentsub}.plotLst=plotLst;
PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
set(handles.edit_plotLst,'string',num2str(plotLst));
set(guiHOMER,'UserData',PMI);
handles.newlist=1;
guidata(hObject, handles);
updatedisplay(handles)


% --- Executes on button press in radio_GLMresult.
function radio_GLMresult_Callback(hObject, eventdata, handles)
% hObject    handle to radio_GLMresult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_GLMresult
load([handles.NIRS.SPM{1},'/SPM.mat'],'-mat')
plotLst = str2num(get(handles.edit_plotLst,'string'));
for indbeta = 1:1%size(SPM.xXn{1}.beta,1)
    Ball = [];
    Tall = [];
    for ibloc = 1:numel(SPM.xXn)
        Tall = [Tall;SPM.xXn{ibloc}.t(indbeta ,plotLst)];
        Ball = [Ball;SPM.xXn{ibloc}.beta(indbeta,plotLst)];
    end
end
figure
subplot(2,1,1)
hold on
bar(mean(Ball))
errorbar(mean(Ball),std(Ball),'x')
xlabel('channel')
ylabel('Beta val')

subplot(2,1,2)
hold on
bar(mean(Tall))
errorbar(mean(Tall),std(Tall),'x')
xlabel('channel')
ylabel('Tval')

figure
subplot(2,1,1)
imagesc(Ball)
xlabel('Blocks')
ylabel('Channel')
title('Beta estimation')
subplot(2,1,2)
imagesc(Tall)
xlabel('Blocks')
ylabel('Channel')
title('T val estimation')
colorbar

updatedisplay(handles)

% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stepvalue_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stepvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stepvalue as text
%        str2double(get(hObject,'String')) returns contents of edit_stepvalue as a double
updatedisplay(handles)

% --- Executes during object creation, after setting all properties.
function edit_stepvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stepvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_moinslst.
function push_moinslst_Callback(hObject, eventdata, handles)
% hObject    handle to push_moinslst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.popup_view,'value')~=5 %Autre mode que zone list
    list = str2num(get(handles.edit_plotLst,'string'));
    step = str2num(get(handles.edit_zonelist,'string'));
    if isempty(step)
        msgbox('Please ajust the increment value one will be place as default')
        set(handles.edit_zonelist,'string','1')
        step = 1;
    end
    list = list-step(1);
    global currentsub;
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    PMI = get(guiHOMER,'UserData');
    cf  = 1;
    if min(list) < 1 % Check if minimum channel is reached and return if next channel cannot be selected.
        % errordlg('Minimum channel reached.','Operation error') % Inform the user that the operation cannot be done.
        return
    end
    set(handles.edit_plotLst,'string',num2str(list));
    PMI{currentsub}.plotLst = list;
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(list,1),PMI{currentsub}.data(cf).MeasList(list,2)];
    handles.newlist=1;
    set(guiHOMER,'UserData',PMI)
    guidata(hObject, handles);
    updatedisplay(handles)
elseif get(handles.popup_view,'value')==5 %Mode zone
    list = str2num(get(handles.edit_zonelist,'string'));
    list = list - 1;
    minlist = min(list);
    if minlist >= 1
        set(handles.edit_zonelist,'string',num2str(list))
    end
end
updatedisplay(handles)

% --- Executes on button press in push_pluslst.
function push_pluslst_Callback(hObject, eventdata, handles)
% hObject    handle to push_pluslst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.popup_view,'value')~=5 %Autre mode que zone list
    list = str2num(get(handles.edit_plotLst,'string'));
    step = str2num(get(handles.edit_zonelist,'string'));
    if isempty(step)
        msgbox('Please ajust the increment value one will be place as default')
        set(handles.edit_zonelist,'string','1')
        step = 1;
    end
    list = list+step(1);
    global currentsub;
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    PMI = get(guiHOMER,'UserData');
    cf  = 1;
    if max(list) > length(PMI{currentsub}.data(cf).MeasList) % Check if maximum channel is reached and return if next channel cannot be selected.
        % errordlg('Maximum channel reached.','Operation error') % Inform the user that the operation cannot be done.
        return
    end
    set(handles.edit_plotLst,'string',num2str(list));
    PMI{currentsub}.plotLst = list;
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(list,1),PMI{currentsub}.data(cf).MeasList(list,2)];
    handles.newlist=1;
    set(guiHOMER,'UserData',PMI)
    guidata(hObject, handles);
    updatedisplay(handles)
elseif get(handles.popup_view,'value')==5 %Mode zone
    list = str2num(get(handles.edit_zonelist,'string'));
    list = list + 1;
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    maxlist = max(list);
    if maxlist <=numel(PMI{currentsub}.zone.plotLst)
        set(handles.edit_zonelist,'string',num2str(list))
    end
end
updatedisplay(handles)


% --- Executes on button press in btn_getzone.
function btn_getzone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_getzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
plotLst_measlist = find(DOT{currentsub}.data(cf).MeasListAct(DOT{currentsub}.plotLst));
plot_measlist =  DOT{currentsub}.plot(plotLst_measlist,:);
if ~isempty(plotLst_measlist) % Make sure that at least one channel is selected to avoid rendering errors later on.
    label = Zonename;
    if ~isempty(label) % Make sure that a zone label is provided.
        if ~isfield(DOT{currentsub},'zone')
            DOT{currentsub}.zone.plot = [{plot_measlist }];
            DOT{currentsub}.zone.plotLst = [{DOT{currentsub}.plotLst(plotLst_measlist) }];
            DOT{currentsub}.zone.label =  [{label} ];
            DOT{currentsub}.zone.color = colorcube(1);
        else
            if ~isfield(DOT{currentsub}.zone,'color')
                num_zone = numel(DOT{currentsub}.zone.label);
                DOT{currentsub}.zone.color = colorcube(num_zone);
            end
            DOT{currentsub}.zone.plot = [DOT{currentsub}.zone.plot {plot_measlist }];
            DOT{currentsub}.zone.plotLst = [DOT{currentsub}.zone.plotLst {DOT{currentsub}.plotLst(plotLst_measlist) }];
            DOT{currentsub}.zone.label =  [DOT{currentsub}.zone.label {label} ];
            DOT{currentsub}.zone.color = [DOT{currentsub}.zone.color;colorcube(1)];
        end

        set(guiHOMER,'UserData',DOT);
        set(handles.popupmenu_zone,'string',DOT{currentsub}.zone.label);
        guidata(handles.figure1,handles);
    else
        errordlg('A name is required to create a zone.', 'Operation error') % Inform the user that a label is required to create a zone.
    end
else
    errordlg('Select at least one channel to create a zone.','Operation error') % Inform the user that the operation wasn't succesful.
end


% --- Executes on button press in btn_deletezone.
function btn_deletezone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_deletezone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
if isfield(DOT{currentsub}, 'zone') % Check if there exists a zone.
    if  ~isempty(DOT{currentsub}.zone.label)
        num = get(handles.popupmenu_zone,'value');
        DOT{currentsub}.zone.plot(num) = [];
        DOT{currentsub}.zone.plotLst(num)= [];
        DOT{currentsub}.zone.label(num) = [];
        DOT{currentsub}.zone.color(num,:)=[];
        if  ~isempty(DOT{currentsub}.zone.label)
            set(guiHOMER,'UserData',DOT);
            set(handles.popupmenu_zone,'value',1);
            set(handles.popupmenu_zone,'string',DOT{currentsub}.zone.label);
            guidata(handles.figure1,handles);
        else
            set(handles.popupmenu_zone,'value',1);
            set(handles.popupmenu_zone,'string',' ');
            if isfield(DOT{currentsub}, 'zone')
                DOT{currentsub} = rmfield(DOT{currentsub}, 'zone');
            end
            set(guiHOMER,'UserData',DOT);
            guidata(handles.figure1,handles);
        end
    else
        set(handles.popupmenu_zone,'value',1);
        set(handles.popupmenu_zone,'string',' ');
        guidata(handles.figure1,handles);
    end
else
    % errordlg('A zone must exist to be deleted.','Operation error')
    % Inform the user (isn't necessary since they already achieved their goal.)
end

% --- Executes during object creation, after setting all properties.
function btn_deletezone_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btn_deletezone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in btn_modifyzone.
function btn_modifyzone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_modifyzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
num = get(handles.popupmenu_zone,'value');
button = questdlg('Are you sure to modify the current zone',DOT{currentsub}.zone.label{num},'Yes','No','Yes');
if strcmp(button,'Yes')
    DOT{currentsub}.zone.plotLst{num} = DOT{currentsub}.plotLst;
    DOT{currentsub}.zone.plot{num} = DOT{currentsub}.plot;
    set(guiHOMER,'UserData',DOT);
end


% --- Executes on button press in btn_savezone.
function btn_savezone_Callback(hObject, eventdata, handles,type)
% hObject    handle to btn_savezone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global currentsub;

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
if type == 1
    %  button = questdlg('Do you want to save zone and reuse it for a different measurement nirs file?','Save option','Yes','No','No')
    %     if strcmp(button,'Yes')
    %         [name,path]=uigetfile('.nirs');
    %         B = load([path,name],'-mat');
    %         ML_new = B.ml;
    %         ML_old = DOT{currentsub}.data(cf).MeasList;
    %         zone = DOT{currentsub}.zone;
    %         for i = 1:numel(DOT{currentsub}.zone.plotLst)
    %             plotLst =  DOT{currentsub}.zone.plotLst{i};
    %             plotLstnew = [];
    %             plotold = DOT{currentsub}.zone.plot{i};
    %             plotnew = [];
    %             for indplot = 1:numel(plotLst)
    %                 a = plotLst(indplot);
    %                 lista = ML_old(a,:);
    %                 newid = find(ML_new(:,1) == lista(1) & ML_new(:,2) == lista(2) & ML_new(:,3) == lista(3) & ML_new(:,4) == lista(4));
    %                 if ~isempty(newid)
    %                     plotLstnew=[plotLstnew;newid];
    %                     plotnew = [plotnew;plotold(indplot,:)];
    %                 end
    %             end
    %             zone.plotLst{i} = plotLstnew;
    %             zone.plot{i}=plotnew;
    %         end
    %         [name, path]= uiputfile([name,'.zone']);
    %         save( [path, name],'zone');
    DOT{currentsub}.zone.ml = DOT{currentsub}.data(cf).MeasList;
    SD.SrcPos = handles.NIRS.Cf.H.S.r.o.mm.p';
    SD.DetPos = handles.NIRS.Cf.H.D.r.o.mm.p';
    SD.Lambda = handles.NIRS.Cf.dev.wl;
    SD.nSrcs  = numel(SD.SrcPos)/3;
    SD.nDets  = numel(SD.DetPos)/3;
    SD.SrcAmp = ones(SD.nSrcs ,2);
    SD.DetAmp = ones(SD.nDets,2);
    %DOT{currentsub}.zone.SD = DOT{currentsub}.SD;
    for id =1:max(size(DOT{currentsub}.data(cf).MeasList))
        idsrc=DOT{currentsub}.zone.ml(id,1);
        iddet=DOT{currentsub}.zone.ml(id,2);
        SrcPos=SD.SrcPos(idsrc,:);
        DetPos=SD.DetPos(iddet,:);
        pos(id,1) = (DetPos(1)-SrcPos(1))/2+SrcPos(1);
        pos(id,2) = (DetPos(2)-SrcPos(2))/2+SrcPos(2);
        pos(id,3) = (DetPos(3)-SrcPos(3))/2+SrcPos(3);
        pos(id,4) = sqrt( (DetPos(1)-SrcPos(1))^2 + (DetPos(2)-SrcPos(2))^2+ (DetPos(3)-SrcPos(3))^2);           %'distance entre les canaux
    end
    DOT{currentsub}.zone.pos = pos;
    DOT{currentsub}.zone.SD = SD;
    [name, path]= uiputfile('.zone');
    zone = DOT{currentsub}.zone;
    save( [path, name],'zone');
    
    
else
    [name, path]= uigetfile('.zone') ;
    try
    load([path, name],'-mat');
    catch
       return
    end 
    if isfield(zone,'ml')
        if numel(zone.ml) ~=  numel(DOT{currentsub}.data(cf).MeasList)
            msgbox('error in zone concordance');return
        end
        if find(zone.ml ~=  DOT{currentsub}.data(cf).MeasList)
            msgbox('error in zone concordance')
        end
    else
        msgbox('Montage concordance could no be check')
    end
    DOT{currentsub}.data(cf).MeasList;
    if isfield(DOT{currentsub},'zone')
        if zone.plot{1}== 0
            for izone = 1:numel(zone.plot)
                zone.plot{izone} = [DOT{currentsub}.data(cf).MeasList(zone.plotLst{izone},1),...
                    DOT{currentsub}.data(cf).MeasList(zone.plotLst{izone},2)];
            end
        end
        DOT{currentsub}.zone.plot = [DOT{currentsub}.zone.plot,zone.plot];
        %Check the plotLst with the measure list
        DOT{currentsub}.zone.plotLst = [DOT{currentsub}.zone.plotLst, zone.plotLst];
        DOT{currentsub}.zone.label = [DOT{currentsub}.zone.label, zone.label];
        DOT{currentsub}.zone.color =  [DOT{currentsub}.zone.color; zone.color];
        set(handles.popupmenu_zone,'string',DOT{currentsub}.zone.label);
        
        
        for i = 1:numel(DOT{currentsub}.zone.plotLst) %Modification pour chaque zone du plotLst en fonction des canaux de plot
            ML_new = DOT{1}.data.MeasList;
            plotLst =  DOT{currentsub}.zone.plotLst{i};
            plotLstnew = [];
            plotold = DOT{currentsub}.zone.plot{i};
            plotnew = []; 
            for indplot = 1:numel(plotLst)
                a = plotLst(indplot);
                newid = find(ML_new(:,1) == plotold(indplot,1) & ML_new(:,2) == plotold(indplot,2) & ML_new(:,3) == 1 & ML_new(:,4) == 1);
                if ~isempty(newid)
                    plotLstnew=[plotLstnew;newid];
                    plotnew = [plotnew;plotold(indplot,:)];
                end
            end
            zone.plotLst{i} = plotLstnew;
            zone.plot{i}=plotnew;
        end
        DOT{currentsub}.zone.plot = zone.plot;
        DOT{currentsub}.zone.plotLst = zone.plotLst;
        set(handles.popupmenu_zone,'string',DOT{currentsub}.zone.label);
    else
        if zone.plot{1}== 0
            for izone = 1:numel(zone.plot)
                zone.plot{izone} = [DOT{currentsub}.data(cf).MeasList(zone.plotLst{izone},1),...
                    DOT{currentsub}.data(cf).MeasList(zone.plotLst{izone},2)];
            end
        end
        DOT{currentsub}.zone.plot = zone.plot;
        DOT{currentsub}.zone.plotLst = zone.plotLst;
        DOT{currentsub}.zone.label = zone.label;
        DOT{currentsub}.zone.color =  zone.color;
        set(handles.popupmenu_zone,'string',DOT{currentsub}.zone.label);
        set(handles.popupmenu_zone,'value',numel(DOT{currentsub}.zone.label));
    end
end

set(guiHOMER,'UserData',DOT);


% --- Executes on button press in btn_openzone.
function btn_openzone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_openzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_CleanZone.
function btn_CleanZone_Callback(hObject, eventdata, handles)
% hObject    handle to btn_CleanZone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData')
cf = DOT{currentsub}.currentFile;
num = get(handles.popupmenu_zone,'value');
DOT{currentsub}.zone.plot = {};
DOT{currentsub}.zone.plotLst = {};
DOT{currentsub}.zone.label = {};
DOT{currentsub}.zone.color = [];
set(handles.popupmenu_zone,'value',1);
set(handles.popupmenu_zone,'string',' ');
set(guiHOMER,'UserData',DOT);



function edit_noisemarker_Callback(hObject, eventdata, handles)
% hObject    handle to edit_noisemarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_noisemarker as text
%        str2double(get(hObject,'String')) returns contents of edit_noisemarker as a double

updatedisplay(handles)

% --- Executes during object creation, after setting all properties.
function edit_noisemarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_noisemarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in listbox_normlist.
function listbox_normlist_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_normlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_normlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_normlist
imodulenorm = [];
for i=1:numel(handles.NIRS.Dt.fir.pp)
    if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
        if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
            imodulenorm  = i;
        end
    end
end
idmodule = get(handles.popupmenu_module, 'value');
epochavg = 0;
if numel(handles.NIRS.Dt.fir.pp(idmodule).pre)>15
    if strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre(1:15), 'Epoch averaging')|strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre(1:15), 'Manual Gui Epoc');
        epochavg = 1;
    end
end

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
fileid=get(handles.popupmenu_file,'value');
PMI = get(guiHOMER,'UserData');

list = [];

if isfield(PMI{currentsub},'zone')
    if ~isempty(PMI{currentsub}.zone.label)
        num = get(handles.popupmenu_zone,'value');
        chzone = PMI{currentsub}.zone.plotLst{num};
    end
else
    set(handles.menu_navigatezone,'checked','off')
end
if idmodule == imodulenorm
    set(handles.popup_view,'value',2)
end
idfile = get(handles.popupmenu_file,'value');
if ~isempty(imodulenorm)
    NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
    if strcmp('off',get(handles.menu_navigatezone,'checked')) %all
        actNormList = ones(size(NormList,1),1); %Active normlist
        if ~get(handles.radio_showrejectedchannels,'value')
            actNormlistchok = zeros(size(NormList,1),1);
            chok_zone = find(handles.NIRS.Cf.H.C.ok(:,fileid));
            for i=1:numel(chok_zone)
                chok = find(NormList(:,1)==chok_zone(i));
                actNormlistchok(chok) = 1;
            end
            actNormList  = actNormList & actNormlistchok;
        end
        
    elseif strcmp('on',get(handles.menu_navigatezone,'checked'))%only zone selected channel
        actNormList = zeros(size(NormList,1),1);
        for i = 1:numel(chzone)
            chok = find(NormList(:,1)==chzone(i));
            actNormList(chok) = 1;
        end
        if ~get(handles.radio_showrejectedchannels,'value')
            actNormlistchok = zeros(size(NormList,1),1);
            chok_zone = find(handles.NIRS.Cf.H.C.ok(:,fileid));
            for i=1:numel(chok_zone)
                chok = find(NormList(:,1)==chok_zone(i));
                actNormlistchok(chok) = 1;
            end
            actNormList  = actNormList & actNormlistchok;
        end
        
    end
    
    NormList_ind = get(handles.listbox_normlist, 'value');
    if actNormList(NormList_ind)%ok display this index value
        
    else %Show the next ok
        indlist =  find(actNormList);
        nextindlist = find(indlist>NormList_ind);
        if ~isempty(nextindlist)
            NormList_ind= indlist(nextindlist(1));
        else
            NormList_ind = indlist(1);
            idfile = idfile + 1;
            listfile = get(handles.popupmenu_file,'string');
            if idfile <= numel(listfile)
                set(handles.popupmenu_file,'value',idfile)
                popupmenu_file_Callback(hObject, eventdata, handles)
            end
            
        end
        set(handles.listbox_normlist, 'value',NormList_ind);
    end
    
    LabelList = get(handles.listbox_normlist, 'string');
    if NormList(NormList_ind,2)==1
        set(handles.radio_on_norm,'enable','off')
        set(handles.radio_off_norm,'enable','off')
    else
        set(handles.radio_on_norm,'enable','on')
        set(handles.radio_off_norm,'enable','on')
    end
    if strfind(LabelList{NormList_ind(end)},'Off')
        set(handles.radio_on_norm,'value',0);
        set(handles.radio_off_norm,'value',1);
    else
        set(handles.radio_on_norm,'value',1);
        set(handles.radio_off_norm,'value',0);
    end
    list = NormList(NormList_ind,1);
end

if ~isempty(list) % Only proceed if a list value was found.
    set(handles.edit_plotLst,'String',num2str(list));
    currrentsub = 1;
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    PMI = get(guiHOMER,'UserData');
    cf  = 1;
    ch1 = list(1);
    if ch1 <= handles.NIRS.Cf.H.C.N/2
        ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
    else
        ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
    end
    list = sort([ch1,ch2]);


    PMI{currentsub}.plotLst = list;
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(list,1),PMI{currentsub}.data(cf).MeasList(list,2)];
    handles.newlist=1;

    if epochavg==0;
        set(handles.edit_time_start,'string',num2str(PMI{currentsub}.data(cf).HRF.tHRF(NormList(NormList_ind(end),2))));
        set(handles.edit_time_stop,'string',num2str(PMI{currentsub}.data(cf).HRF.tHRF(NormList(NormList_ind(end),3))));
    end
    set(guiHOMER,'UserData',PMI);
    guidata(hObject, handles);
    updatedisplay(handles);
else
    errordlg('Artifact review requires an artifact analysis batch module.','Operation error')
end



% --- Executes during object creation, after setting all properties.
function listbox_normlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_normlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Addnorm.
function btn_Addnorm_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Addnorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function uipanel_Norm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_Norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in radio_on_norm.
function radio_on_norm_Callback(hObject, eventdata, handles)
% hObject    handle to radio_on_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_on_norm
if get(handles.radio_on_norm,'value')
    set(handles.radio_off_norm,'value',0);
else
    set(handles.radio_off_norm,'value',1);
end

imodulenorm = [];
for i=1:numel(handles.NIRS.Dt.fir.pp)
    if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
        if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
            imodulenorm  = i;        %Le dernier is plusieur Step Detection
        end
    end
end
idfile = get(handles.popupmenu_file,'value');

NormList_Label = [];
if ~isempty(imodulenorm)
    NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
    NormList_ind = get(handles.listbox_normlist, 'value');
    NormList_Label = get(handles.listbox_normlist, 'string');
    for i=1:numel(NormList_ind)
        Label = NormList_Label{NormList_ind(i)};
        if get(handles.radio_on_norm,'value') %Mettre à ON
            if strfind(Label,'Off')
                NormList_Label{NormList_ind(i)} = [Label(1:end-3),'On'];
                NormList(NormList_ind,5)=1;
            end
        else %Mettre à Off
            if strfind(Label,'On')
                NormList_Label{NormList_ind(i)} = [Label(1:end-2),'Off'];
                NormList(NormList_ind,5)=0;
            end
        end
    end
    handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1} = NormList;
end

if ~isempty(NormList_Label) % Verify if the operation had been succesful.
    set(handles.listbox_normlist, 'string', NormList_Label);
    NIRS = handles.NIRS;
    NIRSpath = get(handles.edit_nirsmat,'string');
    subjectnb = get(handles.edit_nirsmat,'value');
    save(NIRSpath{subjectnb,1},'NIRS','-mat');
    updatedisplay(handles);
    guidata(handles.figure1, handles);
else
    errordlg('Artifact review requires an artifact analysis batch module.','Operation error')
    set(handles.radio_on_norm,'value',0); % Set back the buttons to default.
    set(handles.radio_off_norm,'value',1);
end

% --- Executes on button press in radio_off_norm.
function radio_off_norm_Callback(hObject, eventdata, handles)
% hObject    handle to radio_off_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_off_norm

if get(handles.radio_off_norm,'value')
    set(handles.radio_on_norm,'value',0);
else
    set(handles.radio_on_norm,'value',1);
end
imodulenorm = [];
for i=1:numel(handles.NIRS.Dt.fir.pp)
    if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
        if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
            imodulenorm  = i;        %Le dernier is plusieur Step Detection
        end
    end
end
idfile = get(handles.popupmenu_file,'value');

NormList_Label = [];
if ~isempty(imodulenorm)
    NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
    NormList_ind = get(handles.listbox_normlist, 'value');
    NormList_Label = get(handles.listbox_normlist, 'string');
    for i=1:numel(NormList_ind)
        Label = NormList_Label{NormList_ind(i)};
        if ~get(handles.radio_off_norm,'value') %Mettre à ON
            if strfind(Label,'Off')
                NormList_Label{NormList_ind(i)} = [Label(1:end-3),'On'];
                NormList(NormList_ind,5)=1;
            end
        else %Mettre à Off
            if strfind(Label,'On')
                NormList_Label{NormList_ind(i)} = [Label(1:end-2),'Off'];
                NormList(NormList_ind,5)=0;
            end
        end
    end
    handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1} = NormList;
end

if ~isempty(NormList_Label) % Verify if the operation has been succesful.
    set(handles.listbox_normlist, 'string', NormList_Label);
    NIRS = handles.NIRS;
    NIRSpath = get(handles.edit_nirsmat,'string');
    subjectnb = get(handles.edit_nirsmat,'value');
    save(NIRSpath{subjectnb,1},'NIRS','-mat');

    updatedisplay(handles);
    guidata(handles.figure1, handles);
else
    errordlg('Artifact review requires an artifact analysis batch module.','Operation error')
    set(handles.radio_on_norm,'value',0); % Set back the buttons to default.
    set(handles.radio_off_norm,'value',1);
end


% --- Executes on key press with focus on listbox_normlist and none of its controls.
function listbox_normlist_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox_normlist (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

artefact_keypressfunction(eventdata, handles)


% --- Executes on button press in btn_findbadch.
function btn_findbadch_Callback(hObject, eventdata, handles)
% hObject    handle to btn_findbadch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty( indstart);
        indstart = 1;
    end
    indstart= indstart(end);
end

if isempty(time_stop);
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop);
        indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
    indstop = indstop(end);
end


nch = size(PMI{currentsub}.data(cf).HRF.AvgC,2);
plotLst= find( reshape(sum(PMI{currentsub}.data(cf).HRF.noise(indstart:indstop,:))~=0,nch,1) & reshape(PMI{currentsub}.data(cf).MeasListAct,nch,1))
if ~isempty(plotLst)
    PMI{currentsub}.plotLst = plotLst;
    set(handles.edit_plotLst,'string',mat2str(plotLst))
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
end
handles.newlist=1;
set(guiHOMER,'UserData',PMI)
guidata(hObject, handles);
updatedisplay(handles)
if isempty(plotLst)
    msgbox('Any artefact are found in this period of time')
end

% --- Executes on button press in btn_forcetozero.
function btn_forcetozero_Callback(hObject, eventdata, handles)
% hObject    handle to btn_forcetozero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = 1;
posxstart = str2num(get(handles.edit_time_start,'string'));
indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
indstart = indstart(end);
offset = ones(size(PMI{currentsub}.data(cf).HRF.AvgC,1),1)*PMI{currentsub}.data(cf).HRF.AvgC(indstart,:);
PMI{currentsub}.data(cf).HRF.AvgC = PMI{currentsub}.data(cf).HRF.AvgC-offset;
set(guiHOMER,'UserData',PMI)
updatedisplay(handles)


% --- Executes on button press in radio_colorzone.
function radio_colorzone_Callback(hObject, eventdata, handles)
% hObject    handle to radio_colorzone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_colorzone

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = 1;
if isfield(PMI{currentsub}, 'zone') % Check if a zone exists.
    if isfield(PMI{currentsub}.zone, 'plotLst')
        izonelist = 1:numel(PMI{currentsub}.zone.plotLst);
        if get(handles.radio_colorzone,'value')==1
            color = zeros(size(handles.NIRS.Cf.H.C.id,2),3);
            for izone=izonelist;
                plotLst= PMI{currentsub}.zone.plotLst{izone};
                for izoneplotLst=1:numel(plotLst);
                    color(plotLst(izoneplotLst),:) = PMI{currentsub}.zone.color(izone,:);
                end
            end

            PMI{1}.color = color;
        else
            PMI{1}.color=lines(size(handles.NIRS.Cf.H.C.id,2)); %affichage de couleur des canaux = couleur des zones
        end
        set(guiHOMER,'UserData',PMI);
        handles.newlist=1; %Update 3dhelmet channel display
        updatedisplay(handles);
    else
        set(handles.radio_colorzone, 'value', 0) % Unselect the radio button.
        errordlg('A zone must exist and be selected to enable color.', 'Operation error') % Notify user and proceed as normal.
    end
else
    set(handles.radio_colorzone, 'value', 0) % Unselect the radio button.
    errordlg('A zone must exist and be selected to enable color.', 'Operation error') % Notify user and proceed as normal.
end


% --------------------------------------------------------------------
function remove_ch_interval_Callback(hObject, eventdata, handles)
% hObject    handle to remove_ch_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
fileid=get(handles.popupmenu_file,'value');
x = handles.NIRS.Dt.fir.pp(end).p{fileid};
%Get the co-channel too
ch1 = str2double(get(gco, 'Tag'));
if ch1 <= handles.NIRS.Cf.H.C.N/2
    ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
else
    ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
end
%ici
idmodule = get(handles.popupmenu_module, 'value');
posxstart = str2num(get(handles.edit_time_start,'string'));
posxstop = str2num(get(handles.edit_time_stop,'string'));
indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
PMI{currentsub}.data(cf).HRF.noise(indstart(end):indstop(end),[ch1 ch2])= 1;
set(guiHOMER,'UserData',PMI)
guidata(hObject,handles);
updatedisplay(handles);


% --------------------------------------------------------------------
function restore_ch_interval_Callback(hObject, eventdata, handles)
% hObject    handle to restore_ch_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
fileid=get(handles.popupmenu_file,'value');
x = handles.NIRS.Dt.fir.pp(end).p{fileid};
%Get the co-channel too
ch1 = str2double(get(gco, 'Tag'));
if ch1 <= handles.NIRS.Cf.H.C.N/2
    ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
else
    ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
end
%ici
idmodule = get(handles.popupmenu_module, 'value');
posxstart = str2num(get(handles.edit_time_start,'string'));
posxstop = str2num(get(handles.edit_time_stop,'string'));
indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
PMI{currentsub}.data(cf).HRF.noise(indstart(end):indstop(end),[ch1 ch2])= 0;
set(guiHOMER,'UserData',PMI)
guidata(hObject,handles);
updatedisplay(handles);


% --- Executes on button press in btn_findgoodch.
function btn_findgoodch_Callback(hObject, eventdata, handles)
% hObject    handle to btn_findgoodch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart = 1;
    end
    indstart =indstart(end);
end

if isempty(time_stop);
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
    indstop = indstop(end);
end



nch = size(PMI{currentsub}.data(cf).HRF.AvgC,2);
plotLst = find(reshape(sum(PMI{currentsub}.data(cf).HRF.noise(indstart:indstop,:))==0,nch,1)& reshape(PMI{currentsub}.data(cf).MeasListAct',nch,1)& ~isnan(sum(PMI{currentsub}.data(cf).HRF.AvgC(indstart:indstop,:)))');

if ~isempty(plotLst)
    PMI{currentsub}.plotLst = plotLst;
    set(handles.edit_plotLst,'string',mat2str(plotLst))
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
end
handles.newlist=1;
set(guiHOMER,'UserData',PMI)
guidata(hObject, handles);
updatedisplay(handles)
if isempty(plotLst)
    msgbox('Any artefact are found in this period of time')
end


% --- Executes on button press in btn_yaxis_up.
function btn_yaxis_up_Callback(hObject, eventdata, handles)
% hObject    handle to btn_yaxis_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
min=str2num(get(handles.edit_min,'string'));
max=str2num(get(handles.edit_max,'string'));
offset = str2num(get(handles.edit_yaxis_offset,'string'));
newmin = min+abs(offset);
newmax = max+abs(offset);
set(handles.edit_min,'string',num2str(newmin));
set(handles.edit_max,'string',num2str(newmax));
updatedisplay(handles)


% --- Executes on button press in btn_yaxis_down.
function btn_yaxis_down_Callback(hObject, eventdata, handles)
% hObject    handle to btn_yaxis_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
min=str2num(get(handles.edit_min,'string'));
max=str2num(get(handles.edit_max,'string'));
offset = str2num(get(handles.edit_yaxis_offset,'string'));
newmin = min-abs(offset);
newmax = max-abs(offset);
set(handles.edit_min,'string',num2str(newmin));
set(handles.edit_max,'string',num2str(newmax));
updatedisplay(handles)

function edit_yaxis_offset_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yaxis_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yaxis_offset as text
%        str2double(get(hObject,'String')) returns contents of edit_yaxis_offset as a double


% --- Executes during object creation, after setting all properties.
function edit_yaxis_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yaxis_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_createzone_fromvhdr.
function btn_createzone_fromvhdr_Callback(hObject, eventdata, handles)
% hObject    handle to btn_createzone_fromvhdr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% file = handles.NIRS.Dt.fir.pp(1).p{1};
% fileVHDR = [file(1:end-3),'vhdr'];
% zone = CreateZoneFromVHDR(fileVHDR);
%view
if get(handles.popupmenu_ROIEEG,'value')
    global currentsub
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    PMI = get(guiHOMER,'UserData');
    cf  = 1;
    newmfile = mfilename('fullpath');
[path,name,ext1] = fileparts(newmfile);
load(fullfile(path,'DoubleBananaBrainstormSansT5O1'))
SD.Lambda = handles.NIRS.Cf.dev.wl;
SD.SrcPos = handles.NIRS.Cf.H.S.r.o.mm.p'/100; %en cm
SD.DetPos = handles.NIRS.Cf.H.D.r.o.mm.p'/100;
SD.nSrc = size(SD.SrcPos,1);
SD.nDet = size(SD.DetPos,1);
SD.SrcAmp = ones(SD.nSrc,2);
SD.DetAmp = ones(SD.nDet,2);
ml = PMI{currentsub}.data.MeasList;

for id = 1:size(ml,1)
    idsrc=ml(id,1);
    iddet=ml(id,2);
    SrcPos=SD.SrcPos(idsrc,:);
    DetPos=SD.DetPos(iddet,:);
    %position du centre du canal xyz et  distance en 4
    
    pos(id,1) = (DetPos(1)-SrcPos(1))/2+SrcPos(1);
    pos(id,2) = (DetPos(2)-SrcPos(2))/2+SrcPos(2);
    pos(id,3) = (DetPos(3)-SrcPos(3))/2+SrcPos(3);
    pos(id,4) = sqrt( (DetPos(1)-SrcPos(1))^2 + (DetPos(2)-SrcPos(2))^2+ (DetPos(3)-SrcPos(3))^2);           %'distance entre les canaux
end


%matrice distance ROI EEG BANANA AND POS CENTER CANAL
for idch = 1:size(ml,1)/2
    for idROI = 1:numel(Channel)
        Posch= [pos(idch,1), pos(idch,2), pos(idch,3)];
        Posroi = Channel(1,idROI).Loc;
        if pos(idch ,4)<0.05
            DistROI(idch,idROI) = sqrt( (Posroi(1)-Posch(1))^2 + (Posroi(2)-Posch(2))^2+ (Posroi(3)-Posch(3))^2);
        else
            DistROI(idch,idROI)= 0;
        end
    end
end


%attribuer les canaux a la zone EEG double banana la plus proche.
[C,I] = min(DistROI');
val = find(C==0); % channel trop long a éliminer
if ~isempty(val)
    I(val) = 0;
end
valid = find(C>0.03); %channel trop loin d'une region a eliminer

if ~isempty(valid)
    I(valid) = 0;
end

izone= 0;
for idROI = 1:numel(Channel)
    find(I==1);
    plotLst = find(I==idROI);
    
    if ~isempty(plotLst)
        izone = izone + 1;
        zone.plotLst{izone}= plotLst;
        zone.label{izone} = Channel(1,idROI).Name{1};
        zone.plot{izone} = [ml(plotLst,1), ml(plotLst,2)];
        zone.color(izone,:) = Channel(1,idROI).color;
    end
end


zone.pos = pos;
zone.ml = ml;
zone.SD = SD;
PMI{currentsub}.zone = zone;
[name, path]= uiputfile(['.zone']);
if ~ischar([name, path]) % Check if a file exists and was selected.
    errordlg('Aborting: A .zone file is required to save the zone.','Input error') % Inform the user of a failure.
    return % End the function.
end
save( [path, name],'zone');
set(guiHOMER,'UserData',PMI);
set(handles.popupmenu_zone,'string',zone.label);
set(handles.popupmenu_zone,'value',1);
handles.newlist=1; %Update 3dhelmet channel display
updatedisplay(handles);
end

% --- Executes on selection change in popup_corrwavelenght.
function popup_corrwavelenght_Callback(hObject, eventdata, handles)
% hObject    handle to popup_corrwavelenght (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_corrwavelenght contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_corrwavelenght


% --- Executes during object creation, after setting all properties.
function popup_corrwavelenght_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_corrwavelenght (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LPF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LPF as text
%        str2double(get(hObject,'String')) returns contents of edit_LPF as a double


% --- Executes during object creation, after setting all properties.
function edit_LPF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_pruning.
function btn_pruning_Callback(hObject, eventdata, handles)
% hObject    handle to btn_pruning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt={'Type min or max',...
    'Nb of channel'};
name='Prunning higher channel';
numlines=1;
defaultanswer={'max','5'};
answer=inputdlg(prompt,name,numlines,defaultanswer);

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = 1;
posxstart = str2num(get(handles.edit_time_start,'string'));
posxstop = str2num(get(handles.edit_time_stop,'string'));
if isempty(posxstart)|isempty(posxstop)
    msgbox('Please select a time interval before')
    return
end

indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
if isempty(indstart)
    indstart = 1;
end
indstart = indstart(end);
indstop = indstop(end);
lstok = find(PMI{currentsub}.data(cf).MeasListAct(PMI{currentsub}.plotLst));
plotLst = PMI{currentsub}.plotLst(lstok);

nch = size(PMI{currentsub}.data(cf).HRF.AvgC,2);
[val,ind]=sort(sum(PMI{currentsub}.data(cf).HRF.AvgC(indstart:indstop,plotLst)));
nbselect = str2num(answer{2});
if (numel(ind)>=nbselect)
    if strcmp(answer{1},'max')
        plotLst = plotLst(ind(end-nbselect+1:end));
    else
        plotLst = plotLst(ind(1:nbselect));
    end
else
    msgbox(['Only ', num2str(numel(ind)),' valid  channel are available select mode to find the ', num2str(nbselect), ' ' ,answer{1}]);
    return
end
if ~isempty(plotLst)
    PMI{currentsub}.plotLst = plotLst;
    set(handles.edit_plotLst,'string',mat2str(plotLst))
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
end
handles.newlist=1;
set(guiHOMER,'UserData',PMI)
guidata(hObject, handles);
updatedisplay(handles)
if isempty(plotLst)
    msgbox('Any artefact are found in this period of time')
end


% --- Executes during object creation, after setting all properties.
function btn_TrialSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btn_TrialSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function menu_trialview_select_Callback(hObject, eventdata, handles)
% hObject    handle to menu_trialview_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
module = get(handles.popupmenu_module,'string')
imodule = get(handles.popupmenu_module,'value')
a = module{imodule}
module{imodule} = [a(1:end-4),'  on'];
set(handles.popupmenu_module,'string',module)

% --------------------------------------------------------------------
function menu_trialview_unselect_Callback(hObject, eventdata, handles)
% hObject    handle to menu_trialview_unselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
module = get(handles.popupmenu_module,'string')
imodule = get(handles.popupmenu_module,'value')
a = module{imodule}
module{imodule} = [a(1:end-4),' off'];
set(handles.popupmenu_module,'string',module)

% --------------------------------------------------------------------
function menu_Module_trialview_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Module_trialview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_Manual_Norm_new_Callback(hObject, eventdata, handles)
% hObject    handle to context_Manual_Norm_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.listbox_normlist,'value')
idfile = get(handles.popupmenu_file,'value');
imodulenorm = [];
for i=1:numel(handles.NIRS.Dt.fir.pp)
    if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
        if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
            imodulenorm  = i;        %Le dernier is plusieur Step Detection
        end
    end
end
if ~isempty(imodulenorm)
    NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
end
%ind tstart tstop
time_start = str2num(get(handles.edit_time_start,'string'));
time_stop = str2num(get(handles.edit_time_stop,'string'));
guiHOMER = getappdata(0,'gui_SPMnirsHSJ')
PMI = get(guiHOMER,'UserData');
cf  = 1;
tstart=find(PMI{1}.data(cf).HRF.tHRF<=time_start);
if isempty(tstart);tstart=1;end
tstop=find(PMI{1}.data(cf).HRF.tHRF<=time_stop);
tstart=tstart(end);
tstop=tstop(end);

NormList_ind = get(handles.listbox_normlist, 'value');
ch = NormList(val,1);
allintch=find(NormList(:,1)==ch);
%Determiner Listmoins et Listplus (autre canaux pour recrer la
%liste complete
if allintch(1)>1
    chmoins = allintch(1)-1;
    Listmoins = NormList(1:chmoins,:);
else
    Listmoins = [];
end
if allintch(end)<size(NormList,1)
    chplus = allintch(end)+1
    Listplus = NormList(chplus:end,:);
else
    Listplus = [];
end

NormListpre =NormList(allintch,:);  %debug
startint =NormList(allintch,2);
stopint = NormList(allintch,3);

allint = [startint,stopint];
allint = sort(allint(:));
idinter = find(allint>=tstart & allint<=tstop);
%CAS interval >
%-> combiner les intervals
if ~isempty(idinter)
    if idinter(1)>1
        idmoins = idinter(1)-1;
        allmoins = allint(1:idmoins); %list pre
    else
        allmoins = [];
    end
    if idinter(end)<numel(allint)
        listplus  = idinter(end)+1
        allplus = allint(listplus:end);
    else
        allplus = [];
    end
    
    if mod(numel(allmoins),2)  %Impair fini par un start
        allnew = [allmoins; tstart-1 ; tstart];
    else %Pair fini par un stop
        allnew = [allmoins;tstart;];
    end
    if mod(numel(allplus),2) %Impair commence par un stop
        allnew = [allnew; tstop;tstop+1; allplus];
    else
        allnew = [allnew; tstop; allplus];
    end
    
    Normlistch = [];
    for i = 1:floor(numel(allnew)/2)
        new = [ch, allnew(i*2-1) , allnew(i*2), 0,0,0]
        Normlistch = [Normlistch;new];
    end
    1
else %CAS interval <
    idpre = find(allint<=tstart);
    idpre = idpre(end);
    allmoins = allint(1:idpre);
    allplus =  allint(idpre+1:end);
    allnew = [allmoins;tstart-1;tstart; tstop; tstop+1;allplus];
    Normlistch = [];
    for i = 1:floor(numel(allnew)/2)
        new = [ch, allnew(i*2-1) , allnew(i*2), 0,0,0]
        Normlistch = [Normlistch;new];
    end
end
idnull = find(Normlistch(:,2)==Normlistch(:,3));
if ~isempty(idnull)
    Normlistch(idnull,:) = [];
end

NormList = [Listmoins;Normlistch;Listplus];

%Sauvegarde new list
handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1} = NormList;
NIRS = handles.NIRS;
NIRSpath = get(handles.edit_nirsmat,'string');
subjectnb = get(handles.edit_nirsmat,'value');
save(NIRSpath{subjectnb,1},'NIRS','-mat');

%View new list
for i = 1:size(NormList,1)
    if NormList(i,5)==1
        state = '_On';
    else
        state = '_Off';
    end
    nblist = sprintf('%04.0f',i);
    chname =sprintf('%03.0f',NormList(i,1));
    valname =sprintf('%03.0f',NormList(i,4)*100);
    Label = [nblist,'_',valname,'Ch',chname,state];
    NormListLabel{i}=Label;
end
set(handles.listbox_normlist,'string',NormListLabel);
guidata(handles.figure1, handles);










% --------------------------------------------------------------------
function context_Manual_Norm_Callback(hObject, eventdata, handles)
% hObject    handle to context_Manual_Norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_createinterval_Callback(hObject, eventdata, handles)
% hObject    handle to context_createinterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%loader la liste
guiHOMER = getappdata(0,'gui_SPMnirsHSJ')
PMI = get(guiHOMER,'UserData');

val = get(handles.listbox_normlist,'value');
idfile = get(handles.popupmenu_file,'value');
imodulenorm = [];
for i=1:numel(handles.NIRS.Dt.fir.pp)
    if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
        if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
            imodulenorm  = i;        %Le dernier is plusieur Step Detection
        end
    end
end
if ~isempty(imodulenorm)
    NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
end
NormList_ind = get(handles.listbox_normlist, 'value');
ch = NormList(val,1);
allintch=find(NormList(:,1)==ch);

%Determiner Listmoins et Listplus (autre canaux pour recrer la
%liste complete
if allintch(1)>1
    chmoins = allintch(1)-1;
    Listmoins = NormList(1:chmoins,:);
else
    Listmoins = [];
end
if allintch(end)<size(NormList,1)
    chplus = allintch(end)+1
    Listplus = NormList(chplus:end,:);
else
    Listplus = [];
end

stepmat = PMI{1}.data(1).HRF.noise';
d = PMI{1}.data(1).HRF.AvgC';
Idx = ch;
dmeanintervallist  = [];
dmeaninterval = zeros(size(d));
indbad = find(stepmat(Idx,:)==1);
if numel(indbad)>3
    indbad(end)=[]; %pour éviter les dépassements d'indice
    indbad(1)=[];
    inddiff = diff(indbad);
    inddiff(end+1) = 0;
    step_dur = 0;
    for i = 1:numel(indbad)
        if inddiff(i) == 1
            step_dur = step_dur+1;
        else
            temp_ind = indbad(i)-step_dur;
            temp_dur = step_dur+1;
            step_dur = 0; %reinitialisation
            interval = [temp_ind temp_ind+temp_dur]; % interval = [ind_start ind_end]
            dmeaninterval(Idx, interval(1):interval(2)) = mean(d(Idx, interval(1):interval(2))');
            dmeanintervallist = [dmeanintervallist;Idx,interval(1),interval(2),mean(d(Idx, interval(1):interval(2))'),0,0];
        end
    end
end
indgood = find(stepmat(Idx,:)==0);
if numel(indgood)>3
    %indgood(end)=[];  %pour éviter les dépassements d'indice
    %indgood(1)=[];
    inddiff = diff(indgood);
    inddiff(end+1) = 0; %on force une fin
    step_dur = 0;
    
    for i = 1:numel(indgood)
        if inddiff(i) == 1
            step_dur = step_dur+1;
        else
            temp_ind = indgood(i)-step_dur;
            temp_dur = step_dur+1;
            step_dur = 0; %reinitialisation
            interval = [temp_ind temp_ind+temp_dur]; % interval = [ind_start ind_end]
            %Eviter les dépassements
            if interval(1)<1
                interval(1)=1;
            end
            if interval(2) > size(dmeaninterval,2)
                interval(2)=size(dmeaninterval,2);
            end
            dmeaninterval(Idx, interval(1):interval(2)) = mean(d(Idx, interval(1):interval(2))');
            dmeanintervallist = [dmeanintervallist;Idx,interval(1),interval(2),mean(d(Idx, interval(1):interval(2))'),0,1];
            %champ 1 Le numéro du canal
            %champ 2 Le temps debut de l'interval
            %champ 3 Le temps de fin de l'interval
            %champ 4 La moyenne pour cet interval
            %champ 5 cet interval doit être traité ou
            %champ 6 si c'est un bon interval 1 si
            %artefacté 0
            %nom pour la renormalisation
        end
    end
end

indgood = find(dmeanintervallist(:,6)==1);
dmeanintervallist_good = dmeanintervallist(indgood,:);

NormList =  [Listmoins;dmeanintervallist_good ;Listplus];
%View new list
for i = 1:size(NormList,1)
    if NormList(i,5)==1
        state = '_On';
    else
        state = '_Off';
    end
    nblist = sprintf('%04.0f',i);
    chname =sprintf('%03.0f',NormList(i,1));
    valname =sprintf('%03.0f',NormList(i,4)*100);
    Label = [nblist,'_',valname,'Ch',chname,state];
    NormListLabel{i}=Label;
end

%Sauvegarde new list
handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1} = NormList;
NIRS = handles.NIRS;

NIRSpath = get(handles.edit_nirsmat,'string');
subjectnb = get(handles.edit_nirsmat,'value');
save(NIRSpath{subjectnb,1},'NIRS','-mat');
set(handles.listbox_normlist,'string',NormListLabel);
guidata(handles.figure1, handles);


% --------------------------------------------------------------------
function menu_navigatezone_Callback(hObject, eventdata, handles)
% hObject    handle to menu_navigatezone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.menu_navigatezone,'checked');
if strcmp(val,'on')
    set(handles.menu_navigatezone,'checked','off');
else
    set(handles.menu_navigatezone,'checked','on');
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
artefact_keypressfunction(eventdata, handles)

function artefact_keypressfunction(eventdata, handles)
hObject= handles.listbox_normlist;
if strcmp(get(handles.radio_on_norm,'enable'),'on')
    if strcmp(eventdata.Key,'leftarrow')
        set(handles.radio_on_norm,'value',1);
        radio_on_norm_Callback(hObject, eventdata, handles);
    end
    if strcmp(eventdata.Key,'rightarrow')
        set(handles.radio_off_norm,'value',1);
        radio_off_norm_Callback(hObject, eventdata, handles);
    end
    
end
idmodule = get(handles.popupmenu_module, 'value');
if strcmp(eventdata.Key,'d')   %deleted ch for this event
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    fileid=get(handles.popupmenu_file,'value');
    x = handles.NIRS.Dt.fir.pp(end).p{fileid};
    if size(PMI{currentsub}.plotLst) > 1 % Verify if the channel removal can be done.
        %Get the other wavelenght too
        
    nbhalf=numel(PMI{currentsub}.data(cf).MeasListAct)/2;
    ch1 = PMI{currentsub}.plotLst;%(1);
    idfirst = find(ch1<=nbhalf);
    idlast = find(ch1>nbhalf);
    halflist = zeros(nbhalf,2);
    halflist(ch1(idfirst),:)=1;
    halflist(ch1(idlast)-nbhalf,:)=1;  
   
    if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
         handles.NIRS.Cf.H.C.okavg(find(halflist(:)))=0;
    else
         handles.NIRS.Cf.H.C.ok(find(halflist(:)), fileid) = 0;
    end
  
% old
%         if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
%             handles.NIRS.Cf.H.C.okavg([ch1,ch2])=0;
%         else
%             handles.NIRS.Cf.H.C.ok([ch1 ch2], fileid) = 0;
%         end

        NIRS=handles.NIRS;
        save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');
        PMI{currentsub}.data(cf).MeasListAct([ch1 ch2])=0;
        set(guiHOMER,'UserData',PMI);
        handles.newlist=1; %update helmet channel
        guidata(hObject,handles);
        updatedisplay(handles);
    else
        errordlg('Cannot delete last channel.','Operation error') % Otherwise indicate why the operation cannot be completed.
    end
end

if strcmp(eventdata.Key,'r')   %restore ch for this event
    %RESTORE CH
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    fileid=get(handles.popupmenu_file,'value');
    x = handles.NIRS.Dt.fir.pp(end).p{fileid};
    
    %Get the co-channel too
    nbhalf=numel(PMI{currentsub}.data(cf).MeasListAct)/2;
    ch1 = PMI{currentsub}.plotLst;%(1);
    idfirst = find(ch1<=nbhalf);
    idlast = find(ch1>nbhalf);
    halflist = zeros(nbhalf,2);
    halflist(ch1(idfirst),:)=1;
    halflist(ch1(idlast)-nbhalf,:)=1;
    idmodule = get(handles.popupmenu_module, 'value');
    

  
    if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
         handles.NIRS.Cf.H.C.okavg(find(halflist(:)))=1;
    else
         handles.NIRS.Cf.H.C.ok(find(halflist(:)), fileid) = 1;
    end
    NIRS=handles.NIRS;
    save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');
    PMI{currentsub}.data(cf).MeasListAct(find(halflist(:)))=1;
    set(guiHOMER,'UserData',PMI);
    handles.newlist=1; %update helmet channel
    guidata(hObject,handles);
    updatedisplay(handles);
end
if  strcmp(eventdata.Key,'a')
    button = questdlg('Are you sure you want to remove this channel for all trial ?','Confirm','Yes','No','No');
    if strcmp(button,'Yes')
        guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
        currentsub=1;
        PMI = get(guiHOMER,'UserData');
        cf = PMI{currentsub}.currentFile;   
        measlistact = ones(size(PMI{currentsub}.data.MeasList,1),1)
          tmp = reshape( measlistact,numel(measlistact)/2, 2)
         tmp( PMI{currentsub}.plotLst)=0;
      tmptoremove =  find( [sum(tmp,2)<2; sum(tmp,2)<2])
        if numel(  tmptoremove) ==2 %% numel(PMI{currentsub}.plotLst) == 2 % Check if the channel can be removed.
            %Get the co-channel too
%             ch1 = PMI{currentsub}.plotLst(1);
%             ch2 = PMI{currentsub}.plotLst(2);
            idmodule = get(handles.popupmenu_module, 'value');
            if size(handles.NIRS.Cf.H.C.ok,1)~=size(PMI{currentsub}.data.MeasListAct,1)
                msgbox('DEBUG MEASLIST HIGHER reinitialisation')
                handles.NIRS.Cf.H.C.ok = ones(size(PMI{currentsub}.data.MeasList,1),size(handles.NIRS.Cf.H.C.ok,2))
                PMI{currentsub}.data.MeasListAct = PMI{currentsub}.data.MeasListAct(1:size(PMI{currentsub}.data.MeasList,1))
            end
                
            
            if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
                 handles.NIRS.Cf.H.C.okavg(tmptoremove,:)=0;
            else
                 handles.NIRS.Cf.H.C.ok(tmptoremove,:) = 0;
            end
            
          %  handles.NIRS.Cf.H.C.ok(ch2,:) = 0;
            NIRS=handles.NIRS;
            save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');
            PMI{currentsub}.data(cf).MeasListAct( tmptoremove)=0;
            set(guiHOMER,'UserData',PMI);
            handles.newlist=1; %update helmet channel
            guidata(hObject,handles);
            updatedisplay(handles);
        else
            errordlg('Cannot delete last channel.Select both wavelenght', 'Operation error') % Inform the user why they cannot remove the last channel.
        end
    end
end

if strcmp(eventdata.Key,'t')
    button = questdlg('Are you sure you want to remove this trial for all channel ?','Confirm','Yes','No','No');
    if strcmp(button,'Yes')
        guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
        currentsub=1;
        PMI = get(guiHOMER,'UserData');
        cf = PMI{currentsub}.currentFile;
        fileid=get(handles.popupmenu_file,'value');
        
        
         if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
                 handles.NIRS.Cf.H.C.okavg(:, fileid)=0;
         else
                 handles.NIRS.Cf.H.C.ok(:, fileid) = 0;
         end
        
        NIRS=handles.NIRS;
        save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');
        PMI{currentsub}.data(cf).MeasListAct(:)=0;
        set(guiHOMER,'UserData',PMI);
        handles.newlist=1; %update helmet channel
        guidata(hObject,handles);
        updatedisplay(handles);
    end
end
if strcmp(eventdata.Key,'y')
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    NChalf = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
    plotLst = PMI{currentsub}.plotLst;
    plotLstAll = [];
    for i=1:numel(plotLst)
        chlist = plotLst(i);
        if isempty(find(chlist - NChalf < 1 ))
            plotLstAll = [plotLstAll,chlist, chlist-NChalf]; %add 830 nm channels
        else
            plotLstAll = [plotLstAll,chlist,chlist+NChalf]; %add 690 nm channels
        end
    end
    plotLst = plotLstAll;
    posxstart = str2num(get(handles.edit_time_start,'string'));
    posxstop = str2num(get(handles.edit_time_stop,'string'));
    indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart=1;
    end
    indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
    PMI{currentsub}.data(cf).HRF.noise(indstart(end):indstop(end),plotLst)= 1;
    set(guiHOMER,'UserData',PMI);
    updatedisplay(handles);
end
if strcmp(eventdata.Key,'u')
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    NChalf = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
    plotLst = PMI{currentsub}.plotLst;
    plotLstAll = [];
    for i=1:numel(plotLst)
        chlist = plotLst(i);
        if isempty(find(chlist - NChalf < 1 ))
            plotLstAll = [plotLstAll,chlist, chlist-NChalf]; %add 830 nm channels
        else
            plotLstAll = [plotLstAll,chlist,chlist+NChalf]; %add 690 nm channels
        end
    end
    plotLst = plotLstAll;
    posxstart = str2num(get(handles.edit_time_start,'string'));
    posxstop = str2num(get(handles.edit_time_stop,'string'));
    indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart = 1;
    end
    indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
    
    PMI{currentsub}.data(cf).HRF.noise(indstart(end):indstop(end),plotLst)= 0;
    set(guiHOMER,'UserData',PMI);
    updatedisplay(handles);
end
if strcmp(eventdata.Key,'c')
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    NChalf = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
    plotLst = PMI{currentsub}.plotLst;
    plotLstAll = [];
    for i=1:numel(plotLst)
        chlist = plotLst(i);
        if isempty(find(chlist - NChalf < 1 ))
            plotLstAll = [plotLstAll,chlist, chlist-NChalf]; %add 830 nm channels
        else
            plotLstAll = [plotLstAll,chlist,chlist+NChalf]; %add 690 nm channels
        end
    end
    plotLst = plotLstAll;
    posxstart = str2num(get(handles.edit_time_start,'string'));
    posxstop = str2num(get(handles.edit_time_stop,'string'));
    if isempty(posxstart)
        indstart=1;
    else
        indstart= find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
    end
    if isempty(posxstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    else
        indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
    end
    PMI{currentsub}.data(cf).HRF.noise(:,plotLst)= 0;
    set(guiHOMER,'UserData',PMI);
    updatedisplay(handles);
end


if strcmp(eventdata.Key,'w') %Whole (Whole channels restore for this trial)
    button = questdlg('Are you sure you want to restore whole channels for this trial ?','Confirm','Yes','No','No');
    if strcmp(button,'Yes')
        guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
        currentsub=1;
        PMI = get(guiHOMER,'UserData');
        cf = PMI{currentsub}.currentFile;
        fileid=get(handles.popupmenu_file,'value');
         if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
                 handles.NIRS.Cf.H.C.okavg(:, fileid)=1;
         else
                 handles.NIRS.Cf.H.C.ok(:, fileid) = 1;
         end
 
        NIRS=handles.NIRS;
        save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');
        PMI{currentsub}.data(cf).MeasListAct(:)=1;
        set(guiHOMER,'UserData',PMI);
        handles.newlist=1; %update helmet channel
        guidata(hObject,handles);
        updatedisplay(handles);
    end
end

if strcmp(eventdata.Key,'i')
    context_createinterval_Callback(hObject, eventdata, handles)
end
if strcmp(eventdata.Key,'uparrow')
    valuel = get(handles.listbox_normlist,'value');
    if valuel>1
        valuel = valuel-1;
    else
        button = questdlg('Do you want to process the previous bloc','Previous bloc','yes','no','yes');
        if strcmp(button,'yes')
            listfile = get(handles.popupmenu_file,'string');
            idfile = get(handles.popupmenu_file,'value');
            if idfile>1
                idfile=idfile-1;
                set(handles.popupmenu_file,'value',idfile);
            end
            popupmenu_file_Callback(hObject, eventdata, handles);
            list = get(handles.listbox_normlist,'string');
            valuel=numel(list);
            handles.oldfile = idfile; %get(handles.popupmenu_file, 'Value');
            guidata(hObject, handles);
        end
    end
    set(handles.listbox_normlist,'value',valuel);
    listbox_normlist_Callback(hObject, eventdata, handles);
    
end
if strcmp(eventdata.Key,'downarrow')
    valuel = get(handles.listbox_normlist,'value');
    list = get(handles.listbox_normlist,'string');
    if valuel<numel(list)
        valuel = valuel+1;
    else
        listfile = get(handles.popupmenu_file,'string');
        idfile = get(handles.popupmenu_file,'value');
        if idfile==numel(listfile)
            msgbox('End of the verification','DONE')
            return
        end
        button = questdlg('Do you want to process the next bloc','Next bloc','yes','no','yes');
        if strcmp(button,'yes')
            listfile = get(handles.popupmenu_file,'string');
            idfile = get(handles.popupmenu_file,'value');
            if idfile<numel(listfile)
                idfile=idfile+1;
                set(handles.popupmenu_file,'value',idfile);
            end
            popupmenu_file_Callback(hObject, eventdata, handles);
            handles.oldfile = idfile; %get(handles.popupmenu_file, 'Value');
            guidata(hObject, handles);
            valuel=1;
        end
    end
    set(handles.listbox_normlist,'value',valuel);
    listbox_normlist_Callback(hObject, eventdata, handles);
end



% --------------------------------------------------------------------
function context_warning_noise_Callback(hObject, eventdata, handles)
% hObject    handle to context_warning_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.context_warning_noise,'checked'),'on')
    set(handles.context_warning_noise,'checked','off')
else
    set(handles.context_warning_noise,'checked','on')
end


% --------------------------------------------------------------------
function context_disablewarning_Callback(hObject, eventdata, handles)
% hObject    handle to context_disablewarning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menucontext_CreateIntervalOn_Callback(hObject, eventdata, handles)
% hObject    handle to menucontext_CreateIntervalOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
val = get(handles.listbox_normlist,'value');
idfile = get(handles.popupmenu_file,'value');
imodulenorm = [];

for i=1:numel(handles.NIRS.Dt.fir.pp)
    if numel(handles.NIRS.Dt.fir.pp(i).pre)>=14
        if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
            imodulenorm  = i;        %Le dernier is plusieur Step Detection
        end
    end
end
if ~isempty(imodulenorm)
    NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
end
NormList_ind = get(handles.listbox_normlist, 'value');
plotLst = PMI{1}.plotLst;

for ichplotLst =  1:numel(plotLst)
    
    ch = plotLst(ichplotLst); %NormList(val,1);
    if ch < size(PMI{1}.data.MeasList,1)/2 %Juste les hBO
        allintch=find(NormList(:,1)==ch);
        
        %Determiner Listmoins et Listplus (autre canaux pour recrer la
        %liste complete
        if ~isempty(allintch)
            if allintch(1)>1
                chmoins = allintch(1)-1;
                Listmoins = NormList(1:chmoins,:);
            else
                Listmoins = [];
            end
            if allintch(end)<size(NormList,1)
                chplus = allintch(end)+1;
                Listplus = NormList(chplus:end,:);
            else
                Listplus = [];
            end
            
            stepmat = PMI{1}.data(1).HRF.noise';
            d = PMI{1}.data(1).HRF.AvgC';
            Idx = ch;
            dmeanintervallist  = [];
            dmeaninterval = zeros(size(d));
            indbad = find(stepmat(Idx,:)==1);
            if numel(indbad)>3
                indbad(end)=[]; %pour éviter les dépassements d'indice
                indbad(1)=[];
                inddiff = diff(indbad);
                inddiff(end+1) = 0;
                step_dur = 0;
                for i = 1:numel(indbad)
                    if inddiff(i) == 1
                        step_dur = step_dur+1;
                    else
                        temp_ind = indbad(i)-step_dur;
                        temp_dur = step_dur+1;
                        step_dur = 0; %reinitialisation
                        interval = [temp_ind temp_ind+temp_dur]; % interval = [ind_start ind_end]
                        dmeaninterval(Idx, interval(1):interval(2)) = mean(d(Idx, interval(1):interval(2))');
                        dmeanintervallist = [dmeanintervallist;Idx,interval(1),interval(2),mean(d(Idx, interval(1):interval(2))'),0,0];
                    end
                end
            end
            indgood = find(stepmat(Idx,:)==0);
            if numel(indgood)>3
                %indgood(end)=[];  %pour éviter les dépassements d'indice
                %indgood(1)=[];
                inddiff = diff(indgood);
                inddiff(end+1) = 0; %on force une fin
                step_dur = 0;
                
                for i = 1:numel(indgood)
                    if inddiff(i) == 1
                        step_dur = step_dur+1;
                    else
                        temp_ind = indgood(i)-step_dur;
                        temp_dur = step_dur+1;
                        step_dur = 0; %reinitialisation
                        interval = [temp_ind temp_ind+temp_dur]; % interval = [ind_start ind_end]
                        %Eviter les dépassements
                        if interval(1)<1
                            interval(1)=1;
                        end
                        if interval(2) > size(dmeaninterval,2)
                            interval(2)=size(dmeaninterval,2);
                        end
                        dmeaninterval(Idx, interval(1):interval(2)) = mean(d(Idx, interval(1):interval(2))');
                        dmeanintervallist = [dmeanintervallist;Idx,interval(1),interval(2),mean(d(Idx, interval(1):interval(2))'),0,1];
                        %champ 1 Le numéro du canal
                        %champ 2 Le temps debut de l'interval
                        %champ 3 Le temps de fin de l'interval
                        %champ 4 La moyenne pour cet interval
                        %champ 5 cet interval doit être traité ou
                        %champ 6 si c'est un bon interval 1 si
                        %artefacté 0
                        %nom pour la renormalisation
                    end
                end
            end
        end
        %
        indgood = find(dmeanintervallist(:,6)==1);
        dmeanintervallist_good = dmeanintervallist(indgood,:);
        dmeanintervallist_good(2:end,5)=1;
        NormList =  [Listmoins;dmeanintervallist_good ;Listplus];
    end
end

%View new list
for i = 1:size(NormList,1)
    if NormList(i,5)==1
        state = '_On';
    else
        state = '_Off';
    end
    nblist = sprintf('%04.0f',i);
    chname =sprintf('%03.0f',NormList(i,1));
    valname =sprintf('%03.0f',NormList(i,4)*100);
    Label = [nblist,'_',valname,'Ch',chname,state];
    NormListLabel{i}=Label;
end

%Sauvegarde new list
handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1} = NormList;
NIRS = handles.NIRS;

NIRSpath = get(handles.edit_nirsmat,'string');
subjectnb = get(handles.edit_nirsmat,'value');
save(NIRSpath{subjectnb,1},'NIRS','-mat');
set(handles.listbox_normlist,'string',NormListLabel);
guidata(handles.figure1, handles);


% --------------------------------------------------------------------
function menu_exportdata_img_Callback(hObject, eventdata, handles)
% hObject    handle to menu_exportdata_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ')
PMI = get(guiHOMER,'UserData');
pos = get(handles.display_axes,'CurrentPoint');
posx =pos(1,1);
ptime = find(PMI{1}.data.HRF.tHRF<posx);
d1 = PMI{1}.data.HRF.AvgC(ptime(end),:);

NIRS = handles.NIRS;

if isfield(NIRS.Cf.dev,'mask')
    [Mask.hdr Mask.filetype, Mask.fileprefix, Mask.machine] = load_nii_hdr([NIRS.Cf.dev.mask]);
    [Mask.img,hdr] = load_nii_img(Mask.hdr, Mask.filetype, Mask.fileprefix, Mask.machine);
    listvoxmask=find(Mask.img>0);
    img = single(zeros(Mask.hdr.dime.dim(2),Mask.hdr.dime.dim(3),Mask.hdr.dime.dim(4)));
    img(listvoxmask)=single(d1);%uint16(d1);
    % load_nii_img(Mask.hdr, Mask.filetype, Mask.fileprefix, Mask.machine);
end
file = 'C:\data\Antonin\TESTMASK_EXPORT\GLM\Regresseur_Cerebelum\beta_0001.img'
[Mask.hdr Mask.filetype, Mask.fileprefix, Mask.machine] = load_nii_hdr([file]);
Mask.img = img;
% nii.img = int32(imgmask);
% nii.hdr = hdr;
fileout = ['NEW.img'];
save_nii(Mask, fileout)



% --- Executes on selection change in popupmenu_module_hold.
function popupmenu_module_hold_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_module_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_module_hold contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_module_hold

idmodule = get(handles.popupmenu_module,'value');
%avg option
set(handles.popup_average,'enable','off');
set(handles.popup_average,'value',1);
if numel(handles.NIRS.Dt.fir.pp(idmodule).pre)>15
    if strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre(1:15), 'Epoch averaging');
        set(handles.popup_average,'enable','on');
    end
end

updatedata(handles,2,2)
updatedisplay(handles)

handles.oldmodule = get(handles.popupmenu_module,'value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_module_hold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_module_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_hold.
function radiobutton_hold_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_hold
idmodule = get(handles.popupmenu_module,'value');
%avg option
set(handles.popup_average,'enable','off');
set(handles.popup_average,'value',1);
if numel(handles.NIRS.Dt.fir.pp(idmodule).pre)>15
    if strcmp(handles.NIRS.Dt.fir.pp(idmodule).pre(1:15), 'Epoch averaging');
        set(handles.popup_average,'enable','on');
    end
end



updatedata(handles,2,2)
updatedisplay(handles)

handles.oldmodule = get(handles.popupmenu_module,'value');
guidata(hObject, handles);

if get(handles.popup_view,'val')>1
    msgbox('Please use ''Normal'' view to visualised an overlay')
end

% --- Executes on button press in btn_PARAFAC.
function btn_PARAFAC_Callback(hObject, eventdata, handles)
% hObject    handle to btn_PARAFAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set(handles.radio_Parafac,'value',1)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart = 1;
    end
end

if isempty(time_stop);
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
    
end
indstart = indstart(end);
indstop = indstop(end);
% if isempty(time_start)|isempty(time_stop)
%     set(handles.text_Advice,'string','Please enter time start and time stop to define time window')
%     return
% end

gui_ParafacIO(indstart,indstop);




% --------------------------------------------------------------------
function remove_copy_Callback(hObject, eventdata, handles)
% hObject    handle to remove_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataX = get(gco,'Xdata');
dataY = get(gco,'Ydata');
label = get(gco,'displayname');
tmp = sprintf('%s\t%s\n','Time',label);
for i=1:numel(dataY);
    str = sprintf('%f\t%f\n',dataX(i),dataY(i));
    tmp = [tmp,str];
end
clipboard('copy', tmp);


% --- Executes on button press in btn_substractPARAFAC.
function btn_substractPARAFAC_Callback(hObject, eventdata, handles)
% Substract any component
% hObject    handle to btn_substractPARAFAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
d = PMI{currentsub}.data(cf).HRF.AvgC;
idval = get(handles.popupmethodselected,'val');
listmethod = get(handles.popupmethodselected,'string');
idmodule = get(handles.popupmenu_module,'value');
if idmodule < numel(handles.NIRS.Dt.fir.pp)
    msgbox('modification could be done only on the last manual step')
    return
end
%what to substract!
[pathstr, name, ext] = fileparts(handles.NIRSpath{1});

if idval ~= 5
    msgbox('To subtract a previously identified component, you must choose the option ''component'' on the menu.')
    return
end

if strcmp(listmethod{idval},'Parafac')  %substractPARAFAC
    indt = [PMI{currentsub}.tmpPARAFAC.indt(1):PMI{currentsub}.tmpPARAFAC.indt(2)];%Time indice
    intensnorm = d(indt,:);
    %Detrent DATA segment for centrering
    X = 1:1:size(intensnorm,1);
    Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
    Mb2 =  intensnorm(1,:)'; %offset
    A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
    spar = intensnorm - A;
    A = PMI{currentsub}.tmpPARAFAC.Factors{1};
    B = PMI{currentsub}.tmpPARAFAC.Factors{2};
    C = PMI{currentsub}.tmpPARAFAC.Factors{3};
    listgood = PMI{currentsub}.tmpPARAFAC.listgood;
    ComponentToKeep = PMI{1}.tmpPARAFAC.selected;
    Ac = A(:,ComponentToKeep); Bc = B(:,ComponentToKeep); Cc = C(:,ComponentToKeep);
    % figure;plot(reshape(Xm,[numel(indt),size(d,2)]))
    [Xm]=nmodel(({Ac,Bc,Cc}));
    data = cat(3,d(indt,1:end/2),d(indt,end/2+1:end));
    data(:,listgood,:) = data(:,listgood,:)-Xm;    
    %     if numel(data(:,listgood,:))==numel(Xm)    %check sum protection
    %         data(:,listgood,:) = data(:,listgood,:)-Xm;
    %     else
    %         msgbox('Please update the decomposition, this is not the one for the current data')
    %         return
    %     end
    %
    %save in a structure all apply correction
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    try
        load(fullfile(pathstr,'CorrectionApply.mat'));
        newfile = 0;
    catch
        %donot exist create the stucture
        PARCORR.file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR.filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR.module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR.listgood =  listgood;
        PARCORR.indt = indt;
        PARCORR.data = data(:,listgood,:);
        PARCORR.Xm = Xm;
        PARCORR.FacA = PMI{currentsub}.tmpPARAFAC.Factors{1};
        PARCORR.FacB = PMI{currentsub}.tmpPARAFAC.Factors{2};
        PARCORR.FacC = PMI{currentsub}.tmpPARAFAC.Factors{3};
        PARCORR.ComponentToKeep = ComponentToKeep;
        PARCORR.label= ['PARAFAC' sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR.type = 'PARAFAC';
        FacSpatial = PMI{currentsub}.tmpPARAFAC.Factors{2};
        selected = PMI{currentsub}.tmpPARAFAC.selected;
        PARCORR.topo = FacSpatial(:,selected);
        newfile = 1;
    end
    if newfile == 0
        id = numel(PARCORR);
        PARCORR(id+1).file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR(id+1).filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR(id+1).module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR(id+1).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR(id+1).listgood =  listgood;
        PARCORR(id+1).indt = indt; %indice de temps.
        PARCORR(id+1).data = data(:,listgood,:);
        PARCORR(id+1).Xm = Xm;
        PARCORR(id+1).FacA = PMI{currentsub}.tmpPARAFAC.Factors{1};
        PARCORR(id+1).FacB = PMI{currentsub}.tmpPARAFAC.Factors{2};
        PARCORR(id+1).FacC = PMI{currentsub}.tmpPARAFAC.Factors{3};
        PARCORR(id+1).ComponentToKeep = ComponentToKeep;
        PARCORR(id+1).label= ['PARAFAC' sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR(id+1).type = 'PARAFAC';
        FacSpatial = PMI{currentsub}.tmpPARAFAC.Factors{2};
        selected = PMI{currentsub}.tmpPARAFAC.selected;
        PARCORR(id+1).topo = FacSpatial(:,selected);
    end
    
    for i=1:numel(PARCORR)
        CORRlist{i} = PARCORR(i).label;
    end
    set(handles.listbox_CorrectionDecomposition,'string',CORRlist);
    PMI{currentsub}.data(cf).HRF.AvgC(indt,:) = reshape(data,[numel(indt),size(d,2)]);
    
elseif strcmp(listmethod{idval},'PCA') %substract PCA
    %Detrent DATA segment for centrering
    indt = [PMI{currentsub}.tmpPCA.indt(1):PMI{currentsub}.tmpPCA.indt(2)];%Time indice
    intensnorm = d(indt,:);
    X = 1:1:size(intensnorm,1);
    Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
    Mb2 =  intensnorm(1,:)'; %offset
    A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
    spar = intensnorm - A;
    listgood = PMI{currentsub}.tmpPCA.listgood;
    lstSV = PMI{1}.tmpPCA.selected;
    u =PMI{currentsub}.tmpPCA.u ;
    s =PMI{currentsub}.tmpPCA.s ;
    v =PMI{currentsub}.tmpPCA.v ;
    Xm = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
    data = d(indt,:);
    data(:,listgood) = data(:,listgood)- Xm;
    %     if sum( data(:) )==PMI{currentsub}.tmpPCA.checksumd    %check sum protection
    %         data(:,listgood) = data(:,listgood)- temp;
    %     else
    %         msgbox('Please update the decomposition, this is not the one for the current data')
    %         return
    %     end
    %
    %save in a structure all apply correction    
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    try
        load(fullfile(pathstr,'CorrectionApply.mat'));
        newfile = 0;
    catch
        %donot exist create the stucture
        PARCORR.file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR.filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR.module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR.listgood = listgood;
        PARCORR.indt = indt; %indice de temps.
        PARCORR.data = data(:,listgood,:);
        PARCORR.Xm = Xm;
        PARCORR.u =  PMI{currentsub}.tmpPCA.u;
        PARCORR.s = PMI{currentsub}.tmpPCA.s;
        PARCORR.v = PMI{currentsub}.tmpPCA.v;
        PARCORR.Xm = PARCORR.u*PARCORR.s*PARCORR.v';
        PARCORR.ComponentToKeep = PMI{1}.tmpPCA.selected;
        PARCORR.label= ['PCA' sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR.type = 'PCA';
        PARCORR.topo = s(lstSV,lstSV)*v(:,lstSV)';
        newfile = 1;
    end
    if newfile == 0
        id = numel(PARCORR);
        PARCORR(id+1).file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR(id+1).filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR(id+1).module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR(id+1).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR(id+1).listgood = listgood;
        PARCORR(id+1).indt = indt; %indice de temps.
        PARCORR(id+1).data = data(:,listgood,:);
        PARCORR(id+1).u =  PMI{currentsub}.tmpPCA.u;
        PARCORR(id+1).s = PMI{currentsub}.tmpPCA.s;
        PARCORR(id+1).v = PMI{currentsub}.tmpPCA.v;
        PARCORR(id+1).Xm = PARCORR(id+1).u*PARCORR(id+1).s*PARCORR(id+1).v';
        PARCORR(id+1).ComponentToKeep = PMI{1}.tmpPCA.selected;
        PARCORR(id+1).label= ['PCA' sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR(id+1).type = 'PCA';
        PARCORR(id+1).topo = s(lstSV,lstSV)*v(:,lstSV)';
        newfile = 1;
    end
    PMI{currentsub}.data(cf).HRF.AvgC(indt,:) = data;
elseif strcmp(listmethod{idval},'ICA')    %ICA HERE    
   indt = [PMI{currentsub}.tmpICA.indt(1):PMI{currentsub}.tmpICA.indt(2)];%Time indice
   listgood = PMI{currentsub}.tmpICA.listgood;
   selected = PMI{currentsub}.tmpICA.selected;
   A =  PMI{currentsub}.tmpICA.Factors{1};
   C = PMI{currentsub}.tmpICA.Factors{3};
   Xm = (C(:,PMI{currentsub}.tmpICA.selected)*A(PMI{currentsub}.tmpICA.selected,:))';    
         try
            load(fullfile(pathstr,'CorrectionApply.mat'))
            newfile = 0
        catch
        PARCORR.file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR.filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR.module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR.listgood =  listgood;
        PARCORR.indt =  PMI{currentsub}.tmpICA.indt; %indice de temps.
        PARCORR.data =  PMI{currentsub}.tmpICA.spar;
        PARCORR.Factors =  PMI{currentsub}.tmpICA.Factors ;        
        PARCORR.ComponentToKeep = PMI{1}.tmpICA.selected;
        PARCORR.Xm = Xm;
        PARCORR.label= ['ICA' sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR.type = 'ICA';
        PARCORR.topo = C(:,PMI{currentsub}.tmpICA.selected);
         newfile = 1;
        end
         if newfile == 0
            id = numel(PARCORR);
                    PARCORR(id+1).file= get(handles.popupmenu_file,'value');
                    fileall = get(handles.popupmenu_file,'string');
                    PARCORR(id+1).filestr =  fileall{get(handles.popupmenu_file,'value')};
                    PARCORR(id+1).module  =get(handles.popupmenu_module, 'value');
                    moduleall = get(handles.popupmenu_module,'string');
                    PARCORR(id+1).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
                    PARCORR(id+1).listgood =  listgood;
                    PARCORR(id+1).indt =  PMI{currentsub}.tmpICA.indt; %indice de temps.
                    PARCORR(id+1).data =  PMI{currentsub}.tmpICA.spar;
                    PARCORR(id+1).Factors =  PMI{currentsub}.tmpICA.Factors ;        
                    PARCORR(id+1).ComponentToKeep = PMI{1}.tmpICA.selected;
                    PARCORR(id+1).Xm = Xm;
                    PARCORR(id+1).label= ['ICA' sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
                    PARCORR(id+1).type = 'ICA';
                    PARCORR(id+1).topo = C(:,PMI{currentsub}.tmpICA.selected);
                     newfile = 1;
         end
         PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood)-Xm;  
elseif strcmp(listmethod{idval},'GLM')
    if ~isfield(PMI{currentsub}, 'tmpGLM') % TODO: Explain what this field is for and the consequences of not having it.
        return
    end
    indt = [PMI{currentsub}.tmpGLM.indt(1):PMI{currentsub}.tmpGLM.indt(2)];%Time indice
    listgood = PMI{currentsub}.tmpGLM.listgood;
    lstSV = PMI{1}.tmpGLM.selected;
    beta = PMI{currentsub}.tmpGLM.beta ;
    selected = PMI{currentsub}.tmpGLM.selected;
    idreg = PMI{currentsub}.tmpGLM.idreg;
    Xm = PMI{1}.tmpGLM.AUX.data{idreg(selected)}*beta(selected,:);
    label =  PMI{1}.tmpGLM.AUX.label{idreg(selected)};
    try
        load(fullfile(pathstr,'CorrectionApply.mat'));
        newfile = 0;
    catch
        %donot exist create the stucture
        PARCORR.file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR.filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR.module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR.listgood =  listgood;
        PARCORR.beta =  PMI{currentsub}.tmpGLM.beta;
        PARCORR.std =  PMI{currentsub}.tmpGLM.std;
        PARCORR.AUX = PMI{currentsub}.tmpGLM.AUX;
        PARCORR.indt = indt; %indice de temps.
        PARCOMP.data = d(indt,listgood) ;
        PARCORR.Xm = Xm;
        PARCORR.ComponentToKeep = selected;
        PARCORR.idreg = PMI{1}.tmpGLM.idreg;
        PARCORR.label= ['GLM', label, sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR.type = 'GLM';
        PARCORR.topo =  beta(selected,:);
        newfile = 1;
    end
    if newfile == 0
        id = numel(PARCORR);
        PARCORR(id+1).file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR(id+1).filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR(id+1).module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR(id+1).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR(id+1).listgood =  listgood;
        PARCORR(id+1).beta =  PMI{currentsub}.tmpGLM.beta;
        PARCORR(id+1).std =  PMI{currentsub}.tmpGLM.std;
        PARCORR(id+1).AUX = PMI{currentsub}.tmpGLM.AUX;
        PARCORR(id+1).indt = indt; %indice de temps.
        PARCORR(id+1).data =  d(indt,listgood);
        PARCORR(id+1).Xm = Xm;
        PARCORR(id+1).ComponentToKeep = PMI{1}.tmpGLM.selected;
        PARCORR(id+1).idreg = PMI{1}.tmpGLM.idreg;
        PARCORR(id+1).label= ['GLM',label sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR(id+1).type = 'GLM';
        PARCORR(id+1).topo =  beta(selected,:);
        
    end
    PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) = d(indt,listgood)-Xm;
    
    
elseif strcmp(listmethod{idval},'Offset Adjustment')
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    d = PMI{currentsub}.data(cf).HRF.AvgC;
    time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
    time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
    if isempty(time_start); indstart = 1;else;
        indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
        if isempty(indstart)
            indstart = 1;
        end
    end
    
    if isempty(time_stop);
        indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
    else
        indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
        if isempty(indstop)
            indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
        end
        
    end
    indstart = indstart(end);
    indstop = indstop(end);
    
    %offset
    %pre no change
    
    %PMI{currentsub}.data(cf).HRF.AvgC(1:indstart,:)
    %between detrend
    indt =  indstart:indstop;
    
    %forcer le changement au 2 longueur d'onde
    MeasListActplotLst = zeros(size(d,2),1);
    MeasListActplotLst(PMI{1}.plotLst,:)=1;
    temp = reshape( MeasListActplotLst, size(d,2)/2,2)
    listgood = find([sum(temp'),sum(temp')]);
    
    intensnorm = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood);
    X = 1:1:size(intensnorm,1);
    Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
    Mb2 =  intensnorm(1,:)'; %offset
    A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
    spar = intensnorm - A;
    Xm = A-1;
    PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) = spar + ones(numel(indt),1)*PMI{currentsub}.data(cf).HRF.AvgC(indstart-1,listgood) ;
    %post substract
    PMI{currentsub}.data(cf).HRF.AvgC(indstop+1:end,listgood) =  PMI{currentsub}.data(cf).HRF.AvgC(indstop+1:end,listgood)-...
        ones(size(PMI{currentsub}.data(cf).HRF.AvgC(indstop+1:end,1)))* (PMI{currentsub}.data(cf).HRF.AvgC(indstop+1,listgood) - PMI{currentsub}.data(cf).HRF.AvgC(indstart-1,listgood) );
    set(guiHOMER,'UserData',PMI);
    updatedisplay(handles);
    
    %figure;plot( PMI{currentsub}.data(cf).HRF.AvgC)
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    try
        load(fullfile(pathstr,'CorrectionApply.mat'))
        newfile = 0;
    catch
        PARCORR.file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR.filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR.module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR.listgood =listgood;
        PARCORR.indt = indstart:indstop;%indice de temps.
        PARCORR.Xm = Xm; 
        PARCORR.Offset =  (PMI{currentsub}.data(cf).HRF.AvgC(indstop+1,:) - PMI{currentsub}.data(cf).HRF.AvgC(indstart-1,:) );
        newfile = 1;
        PARCORR.label= ['Offset Adjustment' sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR.type = 'Offset Adjustment';
    end
    if newfile == 0
        id = numel(PARCORR);
        PARCORR(1+id).file= get(handles.popupmenu_file,'value');
        fileall = get(handles.popupmenu_file,'string');
        PARCORR(1+id).filestr =  fileall{get(handles.popupmenu_file,'value')};
        PARCORR(1+id).module  =get(handles.popupmenu_module, 'value');
        moduleall = get(handles.popupmenu_module,'string');
        PARCORR(1+id).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
        PARCORR(1+id).listgood =   listgood;
        PARCORR(1+id).indt = indstart:indstop;
        PARCORR(1+id).Xm = Xm; %indice de temps.
        PARCORR(1+id).Offset =  (PMI{currentsub}.data(cf).HRF.AvgC(indstop+1,:) - PMI{currentsub}.data(cf).HRF.AvgC(indstart-1,:) );
        newfile = 1;
        PARCORR(1+id).label= ['Offset Adjustment' sprintf('%03.0f',size(PARCORR,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
        PARCORR(1+id).type = 'Offset Adjustment';
    end
    
elseif strcmp(listmethod{idval},'Component') %substract selected componenet PCA or PARAFAC
    idcomp = get(handles.listbox_Component,'value');
    try
        load(fullfile(pathstr,'SelectedFactors.mat'))
    catch
        msgbox('No selected factors') 
        return
    end
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    idcomp = get(handles.listbox_Component,'value');
    %substraction from the list update list CorrectionApply
    if strcmp(PARCOMP(idcomp).type,'PARAFAC') %in substract component
        try
            load(fullfile(pathstr,'CorrectionApply.mat'))
            newfile = 0;
        catch
            %donot exist create the stucture
            PARCORR.file= PARCOMP(idcomp).file;
            PARCORR.filestr = PARCOMP(idcomp).filestr;
            PARCORR.module = PARCOMP(idcomp).module;
            PARCORR.modulestr = PARCOMP(idcomp).modulestr;
            PARCORR.listgood = PARCOMP(idcomp).listgood ;
            PARCORR.indt = PARCOMP(idcomp).indt;
            PARCORR.data = PARCOMP(idcomp).data;
            PARCORR.Xm = PARCOMP(idcomp).Xm;
            PARCORR.FacA = PARCOMP(idcomp).FacA;
            PARCORR.FacB = PARCOMP(idcomp).FacB;
            PARCORR.FacC = PARCOMP(idcomp).FacC;
            PARCORR.topo = PARCOMP(idcomp).topo;
            PARCORR.ComponentToKeep = PARCOMP(idcomp).ComponentToKeep;
            PARCORR.label = PARCOMP(idcomp).label;
            PARCORR.type = PARCOMP(idcomp).type;
            newfile = 1
        end
        if newfile == 0
            id = numel(PARCORR);
            PARCORR(id+1).file= PARCOMP(idcomp).file;
            PARCORR(id+1).filestr = PARCOMP(idcomp).filestr;
            PARCORR(id+1).module = PARCOMP(idcomp).module;
            PARCORR(id+1).modulestr = PARCOMP(idcomp).modulestr;
            PARCORR(id+1).listgood = PARCOMP(idcomp).listgood ;
            PARCORR(id+1).indt = PARCOMP(idcomp).indt;
            PARCORR(id+1).data = PARCOMP(idcomp).data;
            PARCORR(id+1).Xm = PARCOMP(idcomp).Xm;
            PARCORR(id+1).FacA = PARCOMP(idcomp).FacA;
            PARCORR(id+1).FacB = PARCOMP(idcomp).FacB;
            PARCORR(id+1).FacC = PARCOMP(idcomp).FacC;
            PARCORR(id+1).topo = PARCOMP(idcomp).topo;
            PARCORR(id+1).ComponentToKeep = PARCOMP(idcomp).ComponentToKeep;
            PARCORR(id+1).label = PARCOMP(idcomp).label;
            PARCORR(id+1).type = PARCOMP(idcomp).type;
        end
        indt = PARCOMP(idcomp).indt(1):PARCOMP(idcomp).indt(end);
        listgood = PARCOMP(idcomp).listgood;
        Xm = PARCOMP(idcomp).Xm;
        d = PMI{currentsub}.data(cf).HRF.AvgC;
        data = cat(3,d(indt,1:end/2),d(indt,end/2+1:end));
        data(:,listgood,:) = data(:,listgood,:)-Xm;
        PMI{currentsub}.data(cf).HRF.AvgC(indt,:) = reshape(data,[numel(indt),size(d,2)]);
        
    elseif strcmp(PARCOMP(idcomp).type,'PCA')%in substract component
        try
            load(fullfile(pathstr,'CorrectionApply.mat'))
            newfile = 0
        catch
            %donot exist create the stucture
            PARCORR.file= PARCOMP(idcomp).file;
            PARCORR.filestr = PARCOMP(idcomp).filestr;
            PARCORR.module = PARCOMP(idcomp).module;
            PARCORR.modulestr = PARCOMP(idcomp).modulestr;
            PARCORR.listgood = PARCOMP(idcomp).listgood ;
            PARCORR.indt = PARCOMP(idcomp).indt;
            PARCORR.data = PARCOMP(idcomp).data;
            PARCORR.u = PARCOMP(idcomp).u;
            PARCORR.s = PARCOMP(idcomp).s;
            PARCORR.v = PARCOMP(idcomp).v;
            PARCORR.topo = PARCOMP(idcomp).topo;
            PARCORR.ComponentToKeep = PARCOMP(idcomp).ComponentToKeep;
            PARCORR.Xm = PARCOMP(idcomp).Xm;
            PARCORR.label = PARCOMP(idcomp).label;
            PARCORR.type = PARCOMP(idcomp).type;
            newfile = 1
        end
        if newfile == 0
            id = numel(PARCORR);
            PARCORR(id+1).file= PARCOMP(idcomp).file;
            PARCORR(id+1).filestr = PARCOMP(idcomp).filestr;
            PARCORR(id+1).module = PARCOMP(idcomp).module;
            PARCORR(id+1).modulestr = PARCOMP(idcomp).modulestr;
            PARCORR(id+1).listgood = PARCOMP(idcomp).listgood ;
            PARCORR(id+1).indt = PARCOMP(idcomp).indt;
            PARCORR(id+1).data = PARCOMP(idcomp).data;
            PARCORR(id+1).u = PARCOMP(idcomp).u;
            PARCORR(id+1).s = PARCOMP(idcomp).s;
            PARCORR(id+1).v = PARCOMP(idcomp).v;
            PARCORR(id+1).topo = PARCOMP(idcomp).topo;
            PARCORR(id+1).ComponentToKeep = PARCOMP(idcomp).ComponentToKeep;
            PARCORR(id+1).Xm = PARCOMP(idcomp).Xm;
            PARCORR(id+1).label = PARCOMP(idcomp).label;
            PARCORR(id+1).type = PARCOMP(idcomp).type;
        end
         indt = PARCOMP(idcomp).indt(1):PARCOMP(idcomp).indt(end);
         listgood = PARCOMP(idcomp).listgood;
         Xm = PARCOMP(idcomp).Xm; 
        PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood)-Xm;
    elseif strcmp(PARCOMP(idcomp).type,'ICA')
        
         try
            load(fullfile(pathstr,'CorrectionApply.mat'))
            newfile = 0
        catch
            %donot exist create the stucture
            PARCORR.file= PARCOMP(idcomp).file;
            PARCORR.filestr = PARCOMP(idcomp).filestr;
            PARCORR.module = PARCOMP(idcomp).module;
            PARCORR.modulestr = PARCOMP(idcomp).modulestr;
            PARCORR.listgood = PARCOMP(idcomp).listgood ;
            PARCORR.indt = PARCOMP(idcomp).indt;
            PARCORR.data = PARCOMP(idcomp).data;
            PARCORR.Factors = PARCOMP(idcomp).Factors ;
            PARCORR.topo = PARCOMP(idcomp).topo;
            PARCORR.ComponentToKeep = PARCOMP(idcomp).ComponentToKeep;
            PARCORR.Xm = PARCOMP(idcomp).Xm;
            PARCORR.label = PARCOMP(idcomp).label;
            PARCORR.type = PARCOMP(idcomp).type;
            newfile = 1
        end
         if newfile == 0
            id = numel(PARCORR);
            PARCORR(id+1).file= PARCOMP(idcomp).file;
            PARCORR(id+1).filestr = PARCOMP(idcomp).filestr;
            PARCORR(id+1).module = PARCOMP(idcomp).module;
            PARCORR(id+1).modulestr = PARCOMP(idcomp).modulestr;
            PARCORR(id+1).listgood = PARCOMP(idcomp).listgood ;
            PARCORR(id+1).indt = PARCOMP(idcomp).indt;
            PARCORR(id+1).data = PARCOMP(idcomp).data;
            PARCORR(id+1).Factors = PARCOMP(idcomp).Factors ;
            PARCORR(id+1).topo = PARCOMP(idcomp).topo;
            PARCORR(id+1).ComponentToKeep = PARCOMP(idcomp).ComponentToKeep;
            PARCORR(id+1).Xm = PARCOMP(idcomp).Xm;
            PARCORR(id+1).label = PARCOMP(idcomp).label;
            PARCORR(id+1).type = PARCOMP(idcomp).type;
         end
          indt = PARCOMP(idcomp).indt(1):PARCOMP(idcomp).indt(end);
         listgood = PARCOMP(idcomp).listgood;
         Xm = PARCOMP(idcomp).Xm; 
        PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood)-Xm;
        %figure;plot( PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) )
    elseif strcmp(PARCOMP(idcomp).type,'GLM')%in substract component
        try
            load(fullfile(pathstr,'CorrectionApply.mat'));
            newfile = 0;
        catch
            %donot exist create the stucture
            PARCORR.file= PARCOMP(idcomp).file;;
            PARCORR.filestr= PARCOMP(idcomp).filestr;
            PARCORR.module  = PARCOMP(idcomp).module;
            PARCORR.modulestr = PARCOMP(idcomp).modulestr;
            PARCORR.listgood  = PARCOMP(idcomp).listgood;
            PARCORR.beta =  PARCOMP(idcomp).beta;
            PARCORR.std =  PARCOMP(idcomp).std;
            PARCORR.AUX = PARCOMP(idcomp).AUX;
            PARCORR.indt = PARCOMP(idcomp).indt; %indice de temps.
            PARCORR.data = PARCOMP(idcomp).data  ;
            PARCORR.Xm = PARCOMP(idcomp).Xm;
            PARCORR.ComponentToKeep = PARCOMP(idcomp).ComponentToKeep ;
            PARCORR.topo = PARCOMP(idcomp).topo;
            PARCORR.idreg = PARCOMP(idcomp).idreg;
            PARCORR.label= PARCOMP(idcomp).label;
            PARCORR.type = 'GLM';
            PARCORR.topo = PARCOMP(idcomp).topo;
            newfile = 1;
        end
        if newfile == 0
            id = numel(PARCORR);
            PARCORR(id+1).file= PARCOMP(idcomp).file;
            PARCORR(id+1).filestr= PARCOMP(idcomp).filestr;
            PARCORR(id+1).module  = PARCOMP(idcomp).module;
            PARCORR(id+1).modulestr = PARCOMP(idcomp).modulestr;
            PARCORR(id+1).listgood  = PARCOMP(idcomp).listgood;
            PARCORR(id+1).beta =  PARCOMP(idcomp).beta;
            PARCORR(id+1).std =  PARCOMP(idcomp).std;
            PARCORR(id+1).AUX = PARCOMP(idcomp).AUX;
            PARCORR(id+1).indt = PARCOMP(idcomp).indt; %indice de temps.
            PARCORR(id+1).data = PARCOMP(idcomp).data  ;
            PARCORR(id+1).Xm = PARCOMP(idcomp).Xm;
            PARCORR(id+1).ComponentToKeep = PARCOMP(idcomp).ComponentToKeep ;
            PARCORR(id+1).topo = PARCOMP(idcomp).topo;
            PARCORR(id+1).idreg = PARCOMP(idcomp).idreg;
            PARCORR(id+1).label= PARCOMP(idcomp).label;
            PARCORR(id+1).type = 'GLM';
            PARCORR(id+1).topo = PARCOMP(idcomp).topo;
            
        end
         indt = PARCOMP(idcomp).indt(1):PARCOMP(idcomp).indt(end);
         listgood = PARCOMP(idcomp).listgood;
         Xm = PARCOMP(idcomp).Xm; 
        PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood)-Xm;
    end
    
    %delete this component from the list
    if numel(PARCOMP)>1
        PARCOMP(idcomp) = [];
        save(fullfile(pathstr,'SelectedFactors.mat'),'PARCOMP');
        for i=1:numel(PARCOMP)
            if iscell(PARCOMP(i).label);
                COMPlist{i} = cat(2,PARCOMP(i).label(:));
            else
                COMPlist{i} = PARCOMP(i).label
            end
        end
        set(handles.listbox_Component,'string',COMPlist);
        if (idcomp-1)>1
            set(handles.listbox_Component,'value',idcomp-1);
        else
            set(handles.listbox_Component,'value',idcomp);
        end
    else
        delete(fullfile(pathstr,'SelectedFactors.mat'));
        set(handles.listbox_Component,'string',' ');
    end
    listbox_Component_Callback(hObject, eventdata, handles)
end
%sauvegarde list correction
for i=1:numel(PARCORR)
    if iscell(PARCORR(i).label)
        CORRlist{i} = cat(2,PARCORR(i).label(:));
    else
        CORRlist{i} = PARCORR(i).label;
    end
end
set(handles.listbox_CorrectionDecomposition,'string',CORRlist);
if numel(PARCORR)>1
    set(handles.listbox_CorrectionDecomposition,'value',numel(PARCORR));
end
save(fullfile(pathstr,'CorrectionApply.mat'),'PARCORR');
d1 = PMI{currentsub}.data(cf).HRF.AvgC;
idfile = get(handles.popupmenu_file,'value');
outfile = handles.NIRS.Dt.fir.pp(end).p{idfile};
fwrite_NIR(outfile,d1');


set(guiHOMER,'UserData',PMI);
updatedisplay(handles);

% --- Executes on selection change in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu13


% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Restore.
function btn_Restore_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set(handles.btn_Restore,'enable','off')
val = get(handles.listbox_CorrectionDecomposition,'value');
label = get(handles.listbox_CorrectionDecomposition,'string');
listmethod = get(handles.popupmethodselected,'string');
valmethod = get(handles.popupmethodselected,'value');
if ~strcmp(listmethod{valmethod},'Component') %Only additive on for component
    msgbox('Restore noise by selecting component')
elseif strcmp(listmethod{valmethod},'Component')
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    try
    load(fullfile(pathstr,'CorrectionApply.mat'));
    catch
        msgbox('To be restore you must have correction apply and selected')
        return
    end
    idcomp = get(handles.listbox_CorrectionDecomposition,'value');
    if strcmp( PARCORR(1,idcomp).type,'Offset Adjustment')
         msgbox('Offset ajustment could not be added')
         return
    else       
        choice = questdlg(['Add', label{val}], ...
            'Yes', ...
            'No');
        switch choice
            case 'No'
                return
            case 'Cancel'
                return
        end
    end
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    %time = PMI{1,1}.data.HRF.tHRF
    %ADD PARAFAC        
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    try
        load(fullfile(pathstr,'CorrectionApply.mat'));
    catch
        
    end
    idcomp = get(handles.listbox_CorrectionDecomposition,'value');
    PARCORR(1,idcomp).file;
    PARCORR(1,idcomp).module;
    listgood = PARCORR(1,idcomp).listgood;
    Xm = PARCORR(1,idcomp).Xm;
    indt = PARCORR(1,idcomp).indt(1):PARCORR(1,idcomp).indt(end);
    d = PMI{currentsub}.data(cf).HRF.AvgC;
    
    switch PARCORR(1,idcomp).type
        case 'PARAFAC'
            data = cat(3,d(indt,1:end/2),d(indt,end/2+1:end));
            data(:,listgood,:) = data(:,listgood,:)+Xm;
            PMI{currentsub}.data(cf).HRF.AvgC(indt,:) = reshape(data,[numel(indt),size(d,2)]);
            try
                load(fullfile(pathstr,'SelectedFactors.mat'))
                newfile = 0;
            catch
                %donot exist create the stucture
                PARCOMP.file= PARCORR(idcomp).file;
                PARCOMP.filestr = PARCORR(idcomp).filestr;
                PARCOMP.module = PARCORR(idcomp).module;
                PARCOMP.modulestr = PARCORR(idcomp).modulestr;
                PARCOMP.listgood = PARCORR(idcomp).listgood ;
                PARCOMP.indt = PARCORR(idcomp).indt;
                PARCOMP.data = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood);
                PARCOMP.Xm = PARCORR(idcomp).Xm;
                PARCOMP.FacA = PARCORR(idcomp).FacA;
                PARCOMP.FacB = PARCORR(idcomp).FacB;
                PARCOMP.FacC = PARCORR(idcomp).FacC;
                PARCOMP.topo = PARCORR(idcomp).topo;
                PARCOMP.ComponentToKeep = PARCORR(idcomp).ComponentToKeep;
                PARCOMP.label = PARCORR(idcomp).label;
                PARCOMP.type = PARCORR(idcomp).type;
                newfile = 1;
            end
            if newfile == 0
                id = numel(PARCOMP);
                PARCOMP(id+1).file= PARCORR(idcomp).file;
                PARCOMP(id+1).filestr = PARCORR(idcomp).filestr;
                PARCOMP(id+1).module = PARCORR(idcomp).module;
                PARCOMP(id+1).modulestr = PARCORR(idcomp).modulestr;
                PARCOMP(id+1).listgood = PARCORR(idcomp).listgood ;
                PARCOMP(id+1).indt = PARCORR(idcomp).indt;
                PARCOMP(id+1).data = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood);
                PARCOMP(id+1).Xm = PARCORR(idcomp).Xm;
                PARCOMP(id+1).FacA = PARCORR(idcomp).FacA;
                PARCOMP(id+1).FacB = PARCORR(idcomp).FacB;
                PARCOMP(id+1).FacC = PARCORR(idcomp).FacC;
                if isfield(PARCORR(idcomp),'topo')
                PARCOMP(id+1).topo = PARCORR(idcomp).topo;
                else
                end
                PARCOMP(id+1).ComponentToKeep = PARCORR(idcomp).ComponentToKeep;
                PARCOMP(id+1).label = PARCORR(idcomp).label;
                PARCOMP(id+1).type = PARCORR(idcomp).type;
            end
            
            
        case 'PCA'
            PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) =  d(indt,listgood) + Xm;
            try
                load(fullfile(pathstr,'SelectedFactors.mat'))
                newfile = 0;
            catch
                %donot exist create the stucture
                PARCOMP.file= PARCORR(idcomp).file;
                PARCOMP.filestr = PARCORR(idcomp).filestr;
                PARCOMP.module = PARCORR(idcomp).module;
                PARCOMP.modulestr = PARCORR(idcomp).modulestr;
                PARCOMP.listgood = PARCORR(idcomp).listgood ;
                PARCOMP.indt = PARCORR(idcomp).indt;
                PARCOMP.data = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood);
                PARCOMP.u = PARCORR(idcomp).u;
                PARCOMP.s = PARCORR(idcomp).s;
                PARCOMP.v = PARCORR(idcomp).v;
                PARCOMP.ComponentToKeep = PARCORR(idcomp).ComponentToKeep;
                PARCOMP.topo = PARCORR(idcomp).topo;
                PARCOMP.Xm = PARCORR(idcomp).Xm;
                PARCOMP.label = PARCORR(idcomp).label;
                PARCOMP.type = PARCORR(idcomp).type;
                newfile = 1;
            end
            if newfile == 0
                id = numel(PARCOMP);
                PARCOMP(id+1).file= PARCORR(idcomp).file;
                PARCOMP(id+1).filestr = PARCORR(idcomp).filestr;
                PARCOMP(id+1).module = PARCORR(idcomp).module;
                PARCOMP(id+1).modulestr = PARCORR(idcomp).modulestr;
                PARCOMP(id+1).listgood = PARCORR(idcomp).listgood ;
                PARCOMP(id+1).indt = PARCORR(idcomp).indt;
                PARCOMP(id+1).data = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood);
                PARCOMP(id+1).u = PARCORR(idcomp).u;
                PARCOMP(id+1).s = PARCORR(idcomp).s;
                PARCOMP(id+1).v = PARCORR(idcomp).v;
                PARCOMP(id+1).ComponentToKeep = PARCORR(idcomp).ComponentToKeep;
                PARCOMP(id+1).topo = PARCORR(idcomp).topo;
                PARCOMP(id+1).Xm = PARCORR(idcomp).Xm;
                PARCOMP(id+1).label = PARCORR(idcomp).label;
                PARCOMP(id+1).type = PARCORR(idcomp).type;
            end
         case 'ICA'
            PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) =  d(indt,listgood) + Xm;
            
            try
                load(fullfile(pathstr,'SelectedFactors.mat'));
                newfile = 0;
            catch
                %donot exist create the stucture
                PARCOMP.file= PARCORR(idcomp).file;
                PARCOMP.filestr= PARCORR(idcomp).filestr;
                PARCOMP.module  = PARCORR(idcomp).module;
                PARCOMP.modulestr = PARCORR(idcomp).modulestr;
                PARCOMP.listgood  = PARCORR(idcomp).listgood;
                PARCOMP.indt = PARCORR(idcomp).indt; %indice de temps.
                PARCOMP.data = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood)  ;
                PARCOMP.Factors = PARCORR(idcomp).Factors ;
                PARCOMP.Xm = PARCORR(idcomp).Xm;
                PARCOMP.ComponentToKeep = PARCORR(idcomp).ComponentToKeep ;
                PARCOMP.label= PARCORR(idcomp).label;
                PARCOMP.type = 'ICA';
                PARCOMP.topo = PARCORR(idcomp).topo;
                newfile = 1;
            end
            if newfile == 0
                id = numel(PARCOMP);
                PARCOMP(id+1).file= PARCORR(idcomp).file;
                PARCOMP(id+1).filestr= PARCORR(idcomp).filestr;
                PARCOMP(id+1).module  = PARCORR(idcomp).module;
                PARCOMP(id+1).modulestr = PARCORR(idcomp).modulestr;
                PARCOMP(id+1).listgood  = PARCORR(idcomp).listgood;    
                PARCOMP(id+1).indt = PARCORR(idcomp).indt; %indice de temps.
                PARCOMP(id+1).data = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood)  ;
                PARCOMP(id+1).Xm = PARCORR(idcomp).Xm;
                PARCOMP(id+1).ComponentToKeep = PARCORR(idcomp).ComponentToKeep ;
                PARCOMP(id+1).Factors = PARCORR(idcomp).Factors ;
                PARCOMP(id+1).label= PARCORR(idcomp).label;
                PARCOMP(id+1).type = 'ICA';
                PARCOMP(id+1).topo = PARCORR(idcomp).topo;
                
            end
            
        case 'GLM'
            PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood) =  d(indt,listgood) + Xm;
            
            try
                load(fullfile(pathstr,'SelectedFactors.mat'));
                newfile = 0;
            catch
                %donot exist create the stucture
                PARCOMP.file= PARCORR(idcomp).file;
                PARCOMP.filestr= PARCORR(idcomp).filestr;
                PARCOMP.module  = PARCORR(idcomp).module;
                PARCOMP.modulestr = PARCORR(idcomp).modulestr;
                PARCOMP.listgood  = PARCORR(idcomp).listgood;
                PARCOMP.beta =  PARCORR(idcomp).beta;
                PARCOMP.std =  PARCORR(idcomp).std;
                PARCOMP.AUX = PARCORR(idcomp).AUX;
                PARCOMP.indt = PARCORR(idcomp).indt; %indice de temps.
                PARCOMP.data = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood)  ;
                PARCOMP.Xm = PARCORR(idcomp).Xm;
                PARCOMP.ComponentToKeep = PARCORR(idcomp).ComponentToKeep ;
                PARCOMP.idreg = PARCORR(idcomp).idreg;
                PARCOMP.label= PARCORR(idcomp).label;
                PARCOMP.type = 'GLM';
                PARCOMP.topo = PARCORR(idcomp).topo;
                newfile = 1;
            end
            if newfile == 0
                id = numel(PARCOMP);
                PARCOMP(id+1).file= PARCORR(idcomp).file;
                PARCOMP(id+1).filestr= PARCORR(idcomp).filestr;
                PARCOMP(id+1).module  = PARCORR(idcomp).module;
                PARCOMP(id+1).modulestr = PARCORR(idcomp).modulestr;
                PARCOMP(id+1).listgood  = PARCORR(idcomp).listgood;
                PARCOMP(id+1).beta = PARCORR(idcomp).beta;
                PARCOMP(id+1).std =  PARCORR(idcomp).std;
                PARCOMP(id+1).AUX = PARCORR(idcomp).AUX;
                PARCOMP(id+1).indt = PARCORR(idcomp).indt; %indice de temps.
                PARCOMP(id+1).data = PMI{currentsub}.data(cf).HRF.AvgC(indt,listgood)  ;
                PARCOMP(id+1).Xm = PARCORR(idcomp).Xm;
                PARCOMP(id+1).ComponentToKeep = PARCORR(idcomp).ComponentToKeep ;
                PARCOMP(id+1).idreg = PARCORR(idcomp).idreg;
                PARCOMP(id+1).label= PARCORR(idcomp).label;
                PARCOMP(id+1).type = 'GLM';
                PARCOMP(id+1).topo = PARCORR(idcomp).topo;
                
            end
            
        case 'Offset Adjustment'
            msgbox('Offset ajustment could not be add')
            return
    end
    
    set(guiHOMER,'UserData',PMI);
    d1 = PMI{currentsub}.data(cf).HRF.AvgC;
    %update list CorrectionApply
    if numel(PARCORR)>1
        PARCORR(idcomp) = [];
        save(fullfile(pathstr,'CorrectionApply.mat'),'PARCORR');
        for i=1:numel(PARCORR)
            CORRlist{i} = PARCORR(i).label;
        end
        set(handles.listbox_CorrectionDecomposition,'string',CORRlist);
        set(handles.listbox_CorrectionDecomposition,'value',1);
    else
        delete(fullfile(pathstr,'CorrectionApply.mat'));
        set(handles.listbox_CorrectionDecomposition,'string',' ');
    end
    
    save(fullfile(pathstr,'SelectedFactors.mat'),'PARCOMP');
    for i=1:numel(PARCOMP)
        COMPlist{i} = PARCOMP(i).label;
    end
    set(handles.listbox_Component,'string',COMPlist);
    set(handles.listbox_Component,'value',1);
    %save new data
    idfile = get(handles.popupmenu_file,'value');
    outfile = handles.NIRS.Dt.fir.pp(end).p{idfile};
    fwrite_NIR(outfile,d1');
    updatedisplay(handles);
end
%set(handles.listbox_CorrectionDecomposition,'value',1)

% --- Executes on button press in radiobutton16.
function radiobutton16_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton16


% --- Executes on button press in btn_ICA.
function btn_ICA_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
currentsub=1;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart = 1;
    end
end

if isempty(time_stop);
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
end
indstop = indstop(end);
indstart = indstart(end);
set(handles.radio_ICA,'value',1)

gui_icaIO(indstart,indstop);

% --- Executes on button press in gui_EEG.
function gui_EEG_Callback(hObject, eventdata, handles)
% hObject    handle to gui_EEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%
time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
currentsub=1;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

ifile = get(handles.popupmenu_file,'value') ;
if isfield(handles.NIRS.Dt,'EEG')
    filelist = handles.NIRS.Dt.EEG.pp(end).p{ifile};
    if isfield(handles.NIRS.Dt.EEG.pp(end),'sync_timesec')
        tstart = handles.NIRS.Dt.EEG.pp(end).sync_timesec{ifile};
        tstop = tstart+PMI{currentsub}.data(cf).HRF.tHRF(end);
        offset = tstart + time_start
            try              
                 set(handles.text_Advice,'string',['Time ',sprintf('%02.0f',offset),' ', filelist ]);
         catch
         %   set(handles.text_Advice,'string',['Time ', num2str(round(offset+time_start)),'s']);
            end
    else
        msgbox('No synchronised by trig')
        return
    end
else
    [file,path] = uigetfile('');
    filelist = [path,file];
    handles.NIRS.Dt.EEG.pp.p{ifile}=filelist;
end

if isempty(time_start)|isempty(time_stop)
    set(handles.text_Advice,'string','Please enter time start and time stop to define time window')
    return
end

timelist = 0;


[pathstr, name, ext] = fileparts(filelist);
if strcmp(ext,'.dat')| strcmp(ext,'.eeg')%Generic Data export analyser
    %read the trig !
    [label_all,ind_dur_ch_all] = read_vmrk_all(fullfile(pathstr,[name,'.vmrk']));
    aux5 = handles.NIRS.Dt.fir(end).aux5{ifile};%trig  indice time
    fs = handles.NIRS.Cf.dev.fs;
    %infoBV = read_vhdr_BV(fullfile(pathstr,[name,'.vhdr']));
    aux5t=aux5(:,2).*1/fs
    idtrig = find(aux5t(:,1)>time_start & aux5t(:,1)<time_stop);
    if ~isempty(idtrig)
        trig=  [aux5(idtrig,1), aux5t(idtrig,1) - time_start]
    else
        trig = [255,1]
    end
    % dall = fopen_NIR(filelist,infoBV.NumberOfChannels);
    [dall,info,label,ind_dur_ch] = fopen_EEG(filelist, tstart, tstop);
    timeEEG = ind_dur_ch_all(:,1)*info.SamplingInterval/1e6;
    tall = info.SamplingInterval/1e6:info.SamplingInterval/1e6:size(dall,1)*info.SamplingInterval/1e6;
    idstart = find(tall<=time_start);
    idstop =  find(tall<=time_stop);
    
    PMI{currentsub}.data(cf).EEG.data= dall(idstart(end):idstop(end),:);
    PMI{currentsub}.data(cf).EEG.time = tall(idstart(end):idstop(end));
    PMI{currentsub}.data(cf).EEG.SamplingRateHz  =1/(info.SamplingInterval/1e6);
    PMI{currentsub}.data(cf).EEG.label = info.name_ele; %infoBV.label;
    PMI{currentsub}.data(cf).EEG.trig = trig;
    set(guiHOMER,'UserData',PMI);
    if 1 %BASIC VIEW EEGLAB DIRECT
        guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
        currentsub=1;
        PMI = get(guiHOMER,'UserData');
        cf = PMI{currentsub}.currentFile;
        var{1,1} = 'srate';
        var{1,2} = PMI{currentsub}.data(cf).EEG.SamplingRateHz;
        var{1,3} = 'title';
        var{1,4} = 'EEG';
        var{1,5} = 'limits';
        var{1,6} = [ 0  PMI{currentsub}.data(cf).EEG.time(end)]; %time in ms
        var{1,6} = [PMI{currentsub}.data(cf).EEG.time(1) PMI{currentsub}.data(cf).EEG.time(end)];
        var{1,7} = 'command';
        var{1,8} ='temp';
        var{1,9}= 'events';
        %     var{1,10}.type =[];
        %     var{1,10}.latency = 1;
        %     var{1,10}.urevent = 1;
        PMI{currentsub}.data(cf).EEG.trig;
        for it=1:size(   PMI{currentsub}.data(cf).EEG.trig,1)
            events(it).type = PMI{currentsub}.data(cf).EEG.trig(it,1);
            events(it).latency = PMI{currentsub}.data(cf).EEG.trig(it,2)  *PMI{currentsub}.data(cf).EEG.SamplingRateHz+PMI{currentsub}.data(cf).EEG.time(1)*PMI{currentsub}.data(cf).EEG.SamplingRateHz;
            events(it).urevent  = PMI{currentsub}.data(cf).EEG.trig(it,2)  *PMI{currentsub}.data(cf).EEG.SamplingRateHz+PMI{currentsub}.data(cf).EEG.time(1)*PMI{currentsub}.data(cf).EEG.SamplingRateHz;
        end
        var{1,10} =  events;
        var{1,11} = 'eloc_file';
        for i=1:numel(PMI{currentsub}.data(cf).EEG.label)
            A(i).labels = PMI{currentsub}.data(cf).EEG.label{i};
            A(i).ref    = ' '
            A(i).theta = [];
            A(i).radius = [];
            A(i).X   = [];
            A(i).Y = [];
            A(i).Z = [];
            A(i).ref    = ' '
            A(i).sph_theta = [];
            A(i).sph_phi = [];
            A(i).sph_rad = [];
            type = '';
            A(i).urchan = [];
        end
        var{1,12}= A;
        var{1,13}='winlength';
        var{1,14}=size(PMI{currentsub}.data(cf).EEG.data,1)*1/PMI{currentsub}.data(cf).EEG.SamplingRateHz; %time lenght
        eegplot(-single(PMI{currentsub}.data(cf).EEG.data(:,:)'),var{1},var{2},var{3},var{4},var{5},var{6},var{7},var{8},var{9},var{10},var{11},var{12},var{13},var{14});
        t = get(gca,'xtick')
        
    else %ADVANCE VIEW TO BE DEVELOP
        gui_EEGIO;
    end
    
end


%
%
% guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
% currentsub=1;
% PMI = get(guiHOMER,'UserData');
% cf = PMI{currentsub}.currentFile;
% d = PMI{currentsub}.data(cf).HRF.AvgC;
% tstart=find(PMI{currentsub}.data(cf).HRF.tHRF<time_start);
% tstop=find(PMI{currentsub}.data(cf).HRF.tHRF<time_stop);
% indt = tstart(end):tstop(end);
% plotLst = PMI{1}.plotLst
% intensnorm = d(indt,:);
% figure
% plot(PMI{currentsub}.data(cf).HRF.tHRF(indt)-PMI{currentsub}.data(cf).HRF.tHRF(indt(1)),...
%     intensnorm(:,plotLst))
% %Detrent DATA segment for centrering
% figure; plot(  timeEEG ,EEG.data')
% dtEEG = timeEEG(1): infoBV.SamplingInterval/1e6 : timeEEG(2)
% timeEEG(1,:)
%
% try
%     eegplot(single(EEG.data'),var{1},var{2},var{3},var{4},var{5},var{6},var{7},var{8},var{9},var{10},var{11},var{12})
% catch
%     msgbox('eegplot is not found in the path')
% end
% eegplot() - Scroll (horizontally and/or vertically) through multichannel data.
%             Allows vertical scrolling through channels and manual marking
%             and unmarking of data stretches or epochs for rejection.
% Usage:
%           >> eegplot(data, 'key1', value1 ...); % use interface buttons, etc.
%      else
%           >> eegplot('noui', data, 'key1', value1 ...); % no user interface;
%                                                         % use for plotting
% Menu items:
%    "Figure > print" - [menu] Print figure in portrait or landscape.
%    "Figure > Edit figure" - [menu] Remove menus and buttons and call up the standard
%                  Matlab figure menu. Select "Tools > Edit" to format the figure
%                  for publication. Command line equivalent: 'noui'
%    "Figure > Accept and Close" - [menu] Same as the bottom-right "Reject" button.
%    "Figure > Cancel and Close" - [menu] Cancel all editing, same as the "Cancel" button.
%    "Display > Marking color" > [Hide|Show] marks" - [menu] Show or hide patches of
%                  background color behind the data. Mark stretches of *continuous*
%                  data (e.g., for rejection) by dragging the mouse horizontally
%                  over the activity. With *epoched* data, click on the selected epochs.
%                  Clicked on a marked region to unmark it. Called from the
%                  command line, marked data stretches or epochs are returned in
%                  the TMPREJ variable in the global workspace *if/when* the "Reject"
%                  button is pressed (see Outputs); called from pop_eegplot() or
%                  eeglab(), the marked data portions are removed from the current
%                  dataset, and the dataset is automatically updated.
%     "Display > Marking color > Choose color" - [menu] Change the background marking
%                  color. The marking color(s) of previously marked trials are preserved.
%                  Called from command line, subsequent functions eegplot2event() or
%                  eegplot2trials() allow processing trials marked with different colors
%                  in the TMPREJ output variable. Command line equivalent: 'wincolor'.
%     "Display > Grid > ..." - [menu] Toggle (on or off) time and/or channel axis grids
%                  in the activity plot. Submenus allow modifications to grid aspects.
%                  Command line equivalents: 'xgrid' / 'ygrid'
%     "Display > Show scale" - [menu] Show (or hide if shown) the scale on the bottom
%                  right corner of the activity window. Command line equivalent: 'scale'
%     "Display > Title" - [menu] Change the title of the figure. Command line equivalent:
%                  'title'
%     "Settings > Time range to display"  - [menu] For continuous EEG data, this item
%                  pops up a query window for entering the number of seconds to display
%                  in the activity window. For epoched data, the query window asks
%                  for the number of epochs to display (this can be fractional).
%                  Command line equivalent: 'winlength'
%     "Settings > Number of channels to display" - [menu] Number of channels to display
%                  in the activity window.  If not all channels are displayed, the
%                  user may scroll through channels using the slider on the left
%                  of the activity plot. Command line equivalent: 'dispchans'
%     "Settings > Channel labels > ..."  - [menu] Use numbers as channel labels or load
%                  a channel location file from disk. If called from the eeglab() menu or
%                  pop_eegplot(), the channel labels of the dataset will be used.
%                  Command line equivalent: 'eloc_file'
%     "Settings > Zoom on/off" - [menu] Toggle Matlab figure zoom on or off for time and
%                  electrode axes. left-click to zoom (x2); right-click to reverse-zoom.
%                  Else, draw a rectange in the activity window to zoom the display into
%                  that region. NOTE: When zoom is on, data cannot be marked for rejection.
%     "Settings > Events" - [menu] Toggle event on or off (assuming events have been
%                  given as input). Press "legend" to pop up a legend window for events.
%
% Display window interface:
%    "Activity plot" - [main window] This axis displays the channel activities.  For
%                  continuous data, the time axis shows time in seconds. For epoched
%                  data, the axis label indicate time within each epoch.
%    "Cancel" - [button] Closes the window and cancels any data rejection marks.
%    "Event types" - [button] pop up a legend window for events.
%    "<<" - [button] Scroll backwards though time or epochs by one window length.
%    "<"  - [button] Scroll backwards though time or epochs by 0.2 window length.
%    "Navigation edit box" - [edit box] Enter a starting time or epoch to jump to.
%    ">"  - [button] Scroll forward though time or epochs by 0.2 window length.
%    ">>" - [button] Scroll forward though time or epochs by one window length.
%    "Chan/Time/Value" - [text] If the mouse is within the activity window, indicates
%                  which channel, time, and activity value the cursor is closest to.
%    "Scale edit box" - [edit box] Scales the displayed amplitude in activity units.
%                  Command line equivalent: 'spacing'
%    "+ / -" - [buttons] Use these buttons to +/- the amplitude scale by 10%.
%    "Reject" - [button] When pressed, save rejection marks and close the figure.
%                  Optional input parameter 'command' is evaluated at that time.
%                  NOTE: This button's label can be redefined from the command line
%                  (see 'butlabel' below). If no processing command is specified
%                  for the 'command' parameter (below), this button does not appear.
%    "Stack/Spread" - [button] "Stack" collapses all channels/activations onto the
%                  middle axis of the plot. "Spread" undoes the operation.
%    "Norm/Denorm" - [button] "Norm" normalizes each channel separately such that all
%                  channels have the same standard deviation without changing original
%                  data/activations under EEG structure. "Denorm" undoes the operation.
%
% Required command line input:
%    data        - Input data matrix, either continuous 2-D (channels,timepoints) or
%                  epoched 3-D (channels,timepoints,epochs). If the data is preceded
%                  by keyword 'noui', GUI control elements are omitted (useful for
%                  plotting data for presentation). A set of power spectra at
%                  each channel may also be plotted (see 'freqlimits' below).
% Optional command line keywords:
%    'srate'      - Sampling rate in Hz {default|0: 256 Hz}
%    'spacing'    - Display range per channel (default|0: max(whole_data)-min(whole_data))
%    'eloc_file'  - Electrode filename (as in  >> topoplot example) to read
%                    ascii channel labels. Else,
%                   [vector of integers] -> Show specified channel numbers. Else,
%                   [] -> Do not show channel labels {default|0 -> Show [1:nchans]}
%    'limits'     - [start end] Time limits for data epochs in ms (for labeling
%                   purposes only).
%    'freqs'      - Vector of frequencies (If data contain  spectral values).
%                   size(data, 2) must be equal to size(freqs,2).
%                   *** This option must be used ALWAYS with 'freqlimits' ***
%    'freqlimits' - [freq_start freq_end] If plotting epoch spectra instead of data, frequency
%                   limits to display spectrum. (Data should contain spectral values).
%                   *** This option must be used ALWAYS with 'freqs' ***
%    'winlength'  - [value] Seconds (or epochs) of data to display in window {default: 5}
%    'dispchans'  - [integer] Number of channels to display in the activity window
%                   {default: from data}.  If < total number of channels, a vertical
%                   slider on the left side of the figure allows vertical data scrolling.
%    'title'      - Figure title {default: none}
%    'plottitle'  - Plot title {default: none}
%    'xgrid'      - ['on'|'off'] Toggle display of the x-axis grid {default: 'off'}
%    'ygrid'      - ['on'|'off'] Toggle display of the y-axis grid {default: 'off'}
%    'ploteventdur' - ['on'|'off'] Toggle display of event duration { default: 'off' }
%    'data2'      - [float array] identical size to the original data and
%                   plotted on top of it.
%
% Additional keywords:
%    'command'    - ['string'] Matlab command to evaluate when the 'REJECT' button is
%                   clicked. The 'REJECT' button is visible only if this parameter is
%                   not empty. As explained in the "Output" section below, the variable
%                   'TMPREJ' contains the rejected windows (see the functions
%                   eegplot2event() and eegplot2trial()).
%    'butlabel'   - Reject button label. {default: 'REJECT'}
%    'winrej'     - [start end R G B e1 e2 e3 ...] Matrix giving data periods to mark
%                    for rejection, each row indicating a different period
%                      [start end] = period limits (in frames from beginning of data);
%                      [R G B] = specifies the marking color;
%                      [e1 e2 e3 ...] = a (1,nchans) logical [0|1] vector giving
%                         channels (1) to mark and (0) not mark for rejection.
%    'color'      - ['on'|'off'|cell array] Plot channels with different colors.
%                   If an RGB cell array {'r' 'b' 'g'}, channels will be plotted
%                   using the cell-array color elements in cyclic order {default:'off'}.
%    'wincolor'   - [color] Color to use to mark data stretches or epochs {default:
%                   [ 0.7 1 0.9] is the default marking color}
%    'events'     - [struct] EEGLAB event structure (EEG.event) to use to show events.
%    'submean'    - ['on'|'off'] Remove channel means in each window {default: 'on'}
%    'position'   - [lowleft_x lowleft_y width height] Position of the figure in pixels.
%    'tag'        - [string] Matlab object tag to identify this eegplot() window (allows
%                    keeping track of several simultaneous eegplot() windows).
%    'children'   - [integer] Figure handle of a *dependent* eegplot() window. Scrolling
%                    horizontally in the master window will produce the same scroll in
%                    the dependent window. Allows comparison of two concurrent datasets,
%                    or of channel and component data from the same dataset.
%    'scale'      - ['on'|'off'] Display the amplitude scale {default: 'on'}.
%    'mocap'      - ['on'|'off'] Display motion capture data in a separate figure.
%                     To run, select an EEG data period in the scolling display using
%                     the mouse. Motion capture (mocap) data should be
%                     under EEG.moredata.mocap.markerPosition in xs, ys and zs fields which are
%                     (number of markers, number of time points) arrays.
%                    {default: 'off'}.
%    'selectcommand' - [cell array] list of 3 commands (strings) to run when the mouse
%                      button is down, when it is moving and when the mouse button is up.
%    'ctrlselectcommand' - [cell array] same as above in conjunction with pressing the
%                      CTRL key.
% Outputs:
%    TMPREJ       -  Matrix (same format as 'winrej' above) placed as a variable in
%                    the global workspace (only) when the REJECT button is clicked.
%                    The command specified in the 'command' keyword argument can use
%                    this variable. (See eegplot2trial() and eegplot2event()).
%
% Author: Arnaud Delorme & Colin Humphries, CNL/Salk Institute, SCCN/INC/UCSD, 1998-2001
%
% See also: eeg_multieegplot(), eegplot2event(), eegplot2trial(), eeglab()

% deprecated
%    'colmodif'   - nested cell array of window colors that may be marked/unmarked. Default
%                   is current color only.

% Copyright (C) 2001 Arnaud Delorme & Colin Humphries, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Note for programmers - Internal variable structure:
% All in g. except for Eposition and Eg.spacingwhich are inside the boxes
% gcf
%    1 - winlength
%    2 - srate
%    3 - children
% 'backeeg' axis
%    1 - trialtag
%    2 - g.winrej
%    3 - nested call flag
% 'eegaxis'
%    1 - data
%    2 - colorlist
%    3 - submean    % on or off, subtract the mean
%    4 - maxfreq    % empty [] if no gfrequency content
% 'buttons hold other informations' Eposition for instance hold the current postition
%all credit EEGLAB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [outvar1] = eegplot(data, varargin); % p1,p2,p3,p4,p5,p6,p7,p8,p9)

% Defaults (can be re-defined):

DEFAULT_PLOT_COLOR = { [0 0 1], [0.7 0.7 0.7]};         % EEG line color
try, icadefs;
    DEFAULT_FIG_COLOR = BACKCOLOR;
    BUTTON_COLOR = GUIBUTTONCOLOR;
catch
    DEFAULT_FIG_COLOR = [1 1 1];
    BUTTON_COLOR =[0.8 0.8 0.8];
end;
DEFAULT_AXIS_COLOR = 'k';         % X-axis, Y-axis Color, text Color
DEFAULT_GRID_SPACING = 1;         % Grid lines every n seconds
DEFAULT_GRID_STYLE = '-';         % Grid line style
YAXIS_NEG = 'off';                % 'off' = positive up
DEFAULT_NOUI_PLOT_COLOR = 'k';    % EEG line color for noui option
%   0 - 1st color in AxesColorOrder
SPACING_EYE = 'on';               % g.spacingI on/off
SPACING_UNITS_STRING = '';        % '\muV' for microvolt optional units for g.spacingI Ex. uV
%MAXEVENTSTRING = 10;
%DEFAULT_AXES_POSITION = [0.0964286 0.15 0.842 0.75-(MAXEVENTSTRING-5)/100];
% dimensions of main EEG axes
ORIGINAL_POSITION = [50 50 800 500];

if nargin < 1
    help eegplot
    return
end

% %%%%%%%%%%%%%%%%%%%%%%%%
% Setup inputs
% %%%%%%%%%%%%%%%%%%%%%%%%

if ~isstr(data) % If NOT a 'noui' call or a callback from uicontrols
    
    try
        options = varargin;
        if ~isempty( varargin ),
            for i = 1:2:numel(options)
                g.(options{i}) = options{i+1};
            end
        else g= []; end;
    catch
        disp('eegplot() error: calling convention {''key'', value, ... } error'); return;
    end;
    
    % Selection of data range If spectrum plot
    if isfield(g,'freqlimits') || isfield(g,'freqs')
        %        % Check  consistency of freqlimits
        %        % Check  consistency of freqs
        
        % Selecting data and freqs
        [temp, fBeg] = min(abs(g.freqs-g.freqlimits(1)));
        [temp, fEnd] = min(abs(g.freqs-g.freqlimits(2)));
        data = data(:,fBeg:fEnd,:);
        g.freqs     = g.freqs(fBeg:fEnd);
        
        % Updating settings
        if ndims(data) == 2, g.winlength = g.freqs(end) - g.freqs(1); end
        g.srate     = length(g.freqs)/(g.freqs(end)-g.freqs(1));
        g.isfreq    = 1;
    end
    
    % push button: create/remove window
    % ---------------------------------
    defdowncom   = 'eegplot(''defdowncom'',   gcbf);'; % push button: create/remove window
    defmotioncom = 'eegplot(''defmotioncom'', gcbf);'; % motion button: move windows or display current position
    defupcom     = 'eegplot(''defupcom'',     gcbf);';
    defctrldowncom = 'eegplot(''topoplot'',   gcbf);'; % CTRL press and motion -> do nothing by default
    defctrlmotioncom = ''; % CTRL press and motion -> do nothing by default
    defctrlupcom = ''; % CTRL press and up -> do nothing by default
    
    try, g.srate; 		    catch, g.srate		= 256; 	end;
    try, g.spacing; 			catch, g.spacing	= 0; 	end;
    try, g.eloc_file; 		catch, g.eloc_file	= 0; 	end; % 0 mean numbered
    try, g.winlength; 		catch, g.winlength	= 5; 	end; % Number of seconds of EEG displayed
    try, g.position; 	    catch, g.position	= ORIGINAL_POSITION; 	end;
    try, g.title; 		    catch, g.title		= ['Scroll activity -- eegplot()']; 	end;
    try, g.plottitle; 		catch, g.plottitle	= ''; 	end;
    try, g.trialstag; 		catch, g.trialstag	= -1; 	end;
    try, g.winrej; 			catch, g.winrej		= []; 	end;
    try, g.command; 			catch, g.command	= ''; 	end;
    try, g.tag; 				catch, g.tag		= 'EEGPLOT'; end;
    try, g.xgrid;		    catch, g.xgrid		= 'off'; end;
    try, g.ygrid;		    catch, g.ygrid		= 'off'; end;
    try, g.color;		    catch, g.color		= 'off'; end;
    try, g.submean;			catch, g.submean	= 'off'; end;
    try, g.children;			catch, g.children	= 0; end;
    try, g.limits;		    catch, g.limits	    = [0 1000*(size(data,2)-1)/g.srate]; end;
    try, g.freqs;            catch, g.freqs	    = []; end;  % Ramon
    try, g.freqlimits;	    catch, g.freqlimits	= []; end;
    try, g.dispchans; 		catch, g.dispchans  = size(data,1); end;
    try, g.wincolor; 		catch, g.wincolor   = [ 0.7 1 0.9]; end;
    try, g.butlabel; 		catch, g.butlabel   = 'REJECT'; end;
    try, g.colmodif; 		catch, g.colmodif   = { g.wincolor }; end;
    try, g.scale; 		    catch, g.scale      = 'on'; end;
    try, g.events; 		    catch, g.events      = []; end;
    try, g.ploteventdur;     catch, g.ploteventdur = 'off'; end;
    try, g.data2;            catch, g.data2      = []; end;
    try, g.plotdata2;        catch, g.plotdata2 = 'off'; end;
    try, g.mocap;		    catch, g.mocap		= 'off'; end; % nima
    try, g.selectcommand;     catch, g.selectcommand     = { defdowncom defmotioncom defupcom }; end;
    try, g.ctrlselectcommand; catch, g.ctrlselectcommand = { defctrldowncom defctrlmotioncom defctrlupcom }; end;
    try, g.datastd;          catch, g.datastd = []; end; %ozgur
    try, g.normed;            catch, g.normed = 0; end; %ozgur
    try, g.envelope;          catch, g.envelope = 0; end;%ozgur
    try, g.maxeventstring;    catch, g.maxeventstring = 10; end; % JavierLC
    try, g.isfreq;            catch, g.isfreq = 0;    end; % Ramon
    
    if strcmpi(g.ploteventdur, 'on'), g.ploteventdur = 1; else g.ploteventdur = 0; end;
    if ndims(data) > 2
        g.trialstag = size(	data, 2);
    end;
    
    gfields = fieldnames(g);
    for index=1:length(gfields)
        switch gfields{index}
            case {'spacing', 'srate' 'eloc_file' 'winlength' 'position' 'title' 'plottitle' ...
                    'trialstag'  'winrej' 'command' 'tag' 'xgrid' 'ygrid' 'color' 'colmodif'...
                    'freqs' 'freqlimits' 'submean' 'children' 'limits' 'dispchans' 'wincolor' ...
                    'maxeventstring' 'ploteventdur' 'butlabel' 'scale' 'events' 'data2' 'plotdata2' 'mocap' 'selectcommand' 'ctrlselectcommand' 'datastd' 'normed' 'envelope' 'isfreq'},;
            otherwise, error(['eegplot: unrecognized option: ''' gfields{index} '''' ]);
        end;
    end;
    
    % g.data=data; % never used and slows down display dramatically - Ozgur 2010
    
    if length(g.srate) > 1
        disp('Error: srate must be a single number'); return;
    end;
    if length(g.spacing) > 1
        disp('Error: ''spacing'' must be a single number'); return;
    end;
    if length(g.winlength) > 1
        disp('Error: winlength must be a single number'); return;
    end;
    if isstr(g.title) > 1
        disp('Error: title must be is a string'); return;
    end;
    if isstr(g.command) > 1
        disp('Error: command must be is a string'); return;
    end;
    if isstr(g.tag) > 1
        disp('Error: tag must be is a string'); return;
    end;
    if length(g.position) ~= 4
        disp('Error: position must be is a 4 elements array'); return;
    end;
    switch lower(g.xgrid)
        case { 'on', 'off' },;
        otherwise disp('Error: xgrid must be either ''on'' or ''off'''); return;
    end;
    switch lower(g.ygrid)
        case { 'on', 'off' },;
        otherwise disp('Error: ygrid must be either ''on'' or ''off'''); return;
    end;
    switch lower(g.submean)
        case { 'on' 'off' };
        otherwise disp('Error: submean must be either ''on'' or ''off'''); return;
    end;
    switch lower(g.scale)
        case { 'on' 'off' };
        otherwise disp('Error: scale must be either ''on'' or ''off'''); return;
    end;
    
    if ~iscell(g.color)
        switch lower(g.color)
            case 'on', g.color = { 'k', 'm', 'c', 'b', 'g' };
            case 'off', g.color = { [ 0 0 0.4] };
            otherwise
                disp('Error: color must be either ''on'' or ''off'' or a cell array');
                return;
        end;
    end;
    if length(g.dispchans) > size(data,1)
        g.dispchans = size(data,1);
    end;
    if ~iscell(g.colmodif)
        g.colmodif = { g.colmodif };
    end;
    if g.maxeventstring>20 % JavierLC
        disp('Error: maxeventstring must be equal or lesser than 20'); return;
    end;
    
    % max event string;  JavierLC
    % ---------------------------------
    MAXEVENTSTRING = g.maxeventstring;
    DEFAULT_AXES_POSITION = [0.0964286 0.15 0.842 0.75-(MAXEVENTSTRING-5)/100];
    
    % convert color to modify into array of float
    % -------------------------------------------
    for index = 1:length(g.colmodif)
        if iscell(g.colmodif{index})
            tmpcolmodif{index} = g.colmodif{index}{1} ...
                + g.colmodif{index}{2}*10 ...
                + g.colmodif{index}{3}*100;
        else
            tmpcolmodif{index} = g.colmodif{index}(1) ...
                + g.colmodif{index}(2)*10 ...
                + g.colmodif{index}(3)*100;
        end;
    end;
    g.colmodif = tmpcolmodif;
    
    [g.chans,g.frames, tmpnb] = size(data);
    g.frames = g.frames*tmpnb;
    
    if g.spacing == 0
        maxindex = min(1000, g.frames);
        stds = std(data(:,1:maxindex),[],2);
        g.datastd = stds;
        stds = sort(stds);
        if length(stds) > 2
            stds = mean(stds(2:end-1));
        else
            stds = mean(stds);
        end;
        g.spacing = stds*3;
        if g.spacing > 10
            g.spacing = round(g.spacing);
        end
        if g.spacing  == 0 | isnan(g.spacing)
            g.spacing = 1; % default
        end;
    end
    
    % set defaults
    % ------------
    g.incallback = 0;
    g.winstatus = 1;
    g.setelectrode  = 0;
    [g.chans,g.frames,tmpnb] = size(data);
    g.frames = g.frames*tmpnb;
    g.nbdat = 1; % deprecated
    g.time  = 0;
    g.elecoffset = 0;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare figure and axes
    % %%%%%%%%%%%%%%%%%%%%%%%%
    
    figh = figure('UserData', g,... % store the settings here
        'Color',DEFAULT_FIG_COLOR, 'name', g.title,...
        'MenuBar','none','tag', g.tag ,'Position',g.position, ...
        'numbertitle', 'off', 'visible', 'off', 'Units', 'Normalized');
    
    pos = get(figh,'position'); % plot relative to current axes
    q = [pos(1) pos(2) 0 0];
    s = [pos(3) pos(4) pos(3) pos(4)]./100;
    clf;
    
    % Plot title if provided
    if ~isempty(g.plottitle)
        h = findobj('tag', 'eegplottitle');
        if ~isempty(h)
            set(h, 'string',g.plottitle);
        else
            h = textsc(g.plottitle, 'title');
            set(h, 'tag', 'eegplottitle');
        end;
    end
    
    % Background axis
    % ---------------
    ax0 = axes('tag','backeeg','parent',figh,...
        'Position',DEFAULT_AXES_POSITION,...
        'Box','off','xgrid','off', 'xaxislocation', 'top', 'Units', 'Normalized');
    
    % Drawing axis
    % ---------------
    YLabels = num2str((1:g.chans)');  % Use numbers as default
    YLabels = flipud(str2mat(YLabels,' '));
    if isfield(g,'limits')
        gstart = g.limits(1);
        gstop = g.limits(2);
    else
        gstart =0;
        gstop = g.winlength;
    end
    ax1 = axes('Position',DEFAULT_AXES_POSITION,...
        'userdata', data, ...% store the data here
        'tag','eegaxis','parent',figh,...%(when in g, slow down display)
        'Box','on','xgrid', g.xgrid,'ygrid', g.ygrid,...
        'gridlinestyle',DEFAULT_GRID_STYLE,...
        'Xlim',[0 g.winlength*g.srate],...
        'xtick',[0:g.srate*DEFAULT_GRID_SPACING:g.winlength*g.srate],...
        'Ylim',[0 (g.chans+1)*g.spacing],...
        'YTick',[0:g.spacing:g.chans*g.spacing],...
        'YTickLabel', YLabels,...
        'XTickLabel',num2str((gstart:DEFAULT_GRID_SPACING:gstop)'),...
        'TickLength',[.005 .005],...
        'Color','none',...
        'XColor',DEFAULT_AXIS_COLOR,...
        'YColor',DEFAULT_AXIS_COLOR);
    
    if isstr(g.eloc_file) | isstruct(g.eloc_file)  % Read in electrode names
        if isstruct(g.eloc_file) & length(g.eloc_file) > size(data,1)
            g.eloc_file(end) = []; % common reference channel location
        end;
        eegplot('setelect', g.eloc_file, ax1);
    end;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up uicontrols
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % positions of buttons
    posbut(1,:) = [ 0.0464    0.0254    0.0385    0.0339 ]; % <<
    posbut(2,:) = [ 0.0924    0.0254    0.0288    0.0339 ]; % <
    posbut(3,:) = [ 0.1924    0.0254    0.0299    0.0339 ]; % >
    posbut(4,:) = [ 0.2297    0.0254    0.0385    0.0339 ]; % >>
    posbut(5,:) = [ 0.1287    0.0203    0.0561    0.0390 ]; % Eposition
    posbut(6,:) = [ 0.4744    0.0236    0.0582    0.0390 ]; % Espacing
    posbut(7,:) = [ 0.2762    0.01    0.0582    0.0390 ]; % elec
    posbut(8,:) = [ 0.3256    0.01    0.0707    0.0390 ]; % g.time
    posbut(9,:) = [ 0.4006    0.01    0.0582    0.0390 ]; % value
    posbut(14,:) = [ 0.2762    0.05    0.0582    0.0390 ]; % elec tag
    posbut(15,:) = [ 0.3256    0.05    0.0707    0.0390 ]; % g.time tag
    posbut(16,:) = [ 0.4006    0.05    0.0582    0.0390 ]; % value tag
    posbut(10,:) = [ 0.5437    0.0458    0.0275    0.0270 ]; % +
    posbut(11,:) = [ 0.5437    0.0134    0.0275    0.0270 ]; % -
    posbut(12,:) = [ 0.6    0.02    0.14    0.05 ]; % cancel
    posbut(13,:) = [-0.15   0.02    0.07    0.05 ]; % cancel
    posbut(17,:) = [-0.06    0.02    0.09    0.05 ]; % events types
    posbut(20,:) = [-0.17   0.15     0.015    0.8 ]; % slider
    posbut(21,:) = [0.738    0.87    0.06      0.048];%normalize
    posbut(22,:) = [0.738    0.93    0.06      0.048];%stack channels(same offset)
    posbut(:,1) = posbut(:,1)+0.2;
    
    % Five move buttons: << < text > >>
    
    u(1) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position', posbut(1,:), ...
        'Tag','Pushbutton1',...
        'string','<<',...
        'Callback',['global in_callback;', ...
        'if isempty(in_callback);in_callback=1;', ...
        '    try eegplot(''drawp'',1);', ...
        '        clear global in_callback;', ...
        '    catch error_struct;', ...
        '        clear global in_callback;', ...
        '        throw(error_struct);', ...
        '    end;', ...
        'else;return;end;']);%James Desjardins 2013/Jan/22
    u(2) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position', posbut(2,:), ...
        'Tag','Pushbutton2',...
        'string','<',...
        'Callback',['global in_callback;', ...
        'if isempty(in_callback);in_callback=1;', ...
        '    try eegplot(''drawp'',2);', ...
        '        clear global in_callback;', ...
        '    catch error_struct;', ...
        '        clear global in_callback;', ...
        '        throw(error_struct);', ...
        '    end;', ...
        'else;return;end;']);%James Desjardins 2013/Jan/22
    u(5) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',[1 1 1], ...
        'Position', posbut(5,:), ...
        'Style','edit', ...
        'Tag','EPosition',...
        'string', fastif(g.trialstag(1) == -1, '0', '1'),...
        'Callback', 'eegplot(''drawp'',0);' );
    u(3) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(3,:), ...
        'Tag','Pushbutton3',...
        'string','>',...
        'Callback',['global in_callback;', ...
        'if isempty(in_callback);in_callback=1;', ...
        '    try eegplot(''drawp'',3);', ...
        '        clear global in_callback;', ...
        '    catch error_struct;', ...
        '        clear global in_callback;', ...
        '        throw(error_struct);', ...
        '    end;', ...
        'else;return;end;']);%James Desjardins 2013/Jan/22
    u(4) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(4,:), ...
        'Tag','Pushbutton4',...
        'string','>>',...
        'Callback',['global in_callback;', ...
        'if isempty(in_callback);in_callback=1;', ...
        '    try eegplot(''drawp'',4);', ...
        '        clear global in_callback;', ...
        '    catch error_struct;', ...
        '        clear global in_callback;', ...
        '        error(error_struct);', ...
        '    end;', ...
        'else;return;end;']);%James Desjardins 2013/Jan/22
    
    % Text edit fields: ESpacing
    
    u(6) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',[1 1 1], ...
        'Position', posbut(6,:), ...
        'Style','edit', ...
        'Tag','ESpacing',...
        'string',num2str(g.spacing),...
        'Callback', 'eegplot(''draws'',0);' );
    
    % Slider for vertical motion
    u(20) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position', posbut(20,:), ...
        'Style','slider', ...
        'visible', 'off', ...
        'sliderstep', [0.9 1], ...
        'Tag','eegslider', ...
        'callback', [ 'tmpg = get(gcbf, ''userdata'');' ...
        'tmpg.elecoffset = get(gcbo, ''value'')*(tmpg.chans-tmpg.dispchans);' ...
        'set(gcbf, ''userdata'', tmpg);' ...
        'eegplot(''drawp'',0);' ...
        'clear tmpg;' ], ...
        'value', 0);
    
    % Channels, position, value and tag
    
    u(9) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',DEFAULT_FIG_COLOR, ...
        'Position', posbut(7,:), ...
        'Style','text', ...
        'Tag','Eelec',...
        'string',' ');
    u(10) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',DEFAULT_FIG_COLOR, ...
        'Position', posbut(8,:), ...
        'Style','text', ...
        'Tag','Etime',...
        'string','0.00');
    u(11) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',DEFAULT_FIG_COLOR, ...
        'Position',posbut(9,:), ...
        'Style','text', ...
        'Tag','Evalue',...
        'string','0.00');
    
    u(14)= uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',DEFAULT_FIG_COLOR, ...
        'Position', posbut(14,:), ...
        'Style','text', ...
        'Tag','Eelecname',...
        'string','Chan.');
    
    % Values of time/value and freq/power in GUI
    if g.isfreq
        u15_string =  'Freq';
        u16_string  = 'Power';
    else
        u15_string =  'Time';
        u16_string  = 'Value';
    end
    
    u(15) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',DEFAULT_FIG_COLOR, ...
        'Position', posbut(15,:), ...
        'Style','text', ...
        'Tag','Etimename',...
        'string',u15_string);
    
    u(16) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'BackgroundColor',DEFAULT_FIG_COLOR, ...
        'Position',posbut(16,:), ...
        'Style','text', ...
        'Tag','Evaluename',...
        'string',u16_string);
    
    % ESpacing buttons: + -
    u(7) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(10,:), ...
        'Tag','Pushbutton5',...
        'string','+',...
        'FontSize',8,...
        'Callback','eegplot(''draws'',1)');
    u(8) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(11,:), ...
        'Tag','Pushbutton6',...
        'string','-',...
        'FontSize',8,...
        'Callback','eegplot(''draws'',2)');
    
    cb_normalize = ['g = get(gcbf,''userdata'');if g.normed, disp(''Denormalizing...''); else, disp(''Normalizing...''); end;'...
        'hmenu = findobj(gcf, ''Tag'', ''Normalize_menu'');' ...
        'ax1 = findobj(''tag'',''eegaxis'',''parent'',gcbf);' ...
        'data = get(ax1,''UserData'');' ...
        'if isempty(g.datastd), g.datastd = std(data(:,1:min(1000,g.frames),[],2)); end;'...
        'if g.normed, '...
        'for i = 1:size(data,1), '...
        'data(i,:,:) = data(i,:,:)*g.datastd(i);'...
        'if ~isempty(g.data2), g.data2(i,:,:) = g.data2(i,:,:)*g.datastd(i);end;'...
        'end;'...
        'set(gcbo,''string'', ''Norm'');set(findobj(''tag'',''ESpacing'',''parent'',gcbf),''string'',num2str(g.oldspacing));' ...
        'else, for i = 1:size(data,1),'...
        'data(i,:,:) = data(i,:,:)/g.datastd(i);'...
        'if ~isempty(g.data2), g.data2(i,:,:) = g.data2(i,:,:)/g.datastd(i);end;'...
        'end;'...
        'set(gcbo,''string'', ''Denorm'');g.oldspacing = g.spacing;set(findobj(''tag'',''ESpacing'',''parent'',gcbf),''string'',''5'');end;' ...
        'g.normed = 1 - g.normed;' ...
        'eegplot(''draws'',0);'...
        'set(hmenu, ''Label'', fastif(g.normed,''Denormalize channels'',''Normalize channels''));' ...
        'set(gcbf,''userdata'',g);set(ax1,''UserData'',data);clear ax1 g data;' ...
        'eegplot(''drawp'',0);' ...
        'disp(''Done.'')'];
    % Button for Normalizing data
    u(21) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(21,:), ...
        'Tag','Norm',...
        'string','Norm', 'callback', cb_normalize);
    
    cb_envelope = ['g = get(gcbf,''userdata'');'...
        'hmenu = findobj(gcf, ''Tag'', ''Envelope_menu'');' ...
        'g.envelope = ~g.envelope;' ...
        'set(gcbf,''userdata'',g);'...
        'set(gcbo,''string'',fastif(g.envelope,''Spread'',''Stack''));' ...
        'set(hmenu, ''Label'', fastif(g.envelope,''Spread channels'',''Stack channels''));' ...
        'eegplot(''drawp'',0);clear g;'];
    
    % Button to plot envelope of data
    u(22) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(22,:), ...
        'Tag','Envelope',...
        'string','Stack', 'callback', cb_envelope);
    
    
    if isempty(g.command) tmpcom = 'fprintf(''Rejections saved in variable TMPREJ\n'');';
    else tmpcom = g.command;
    end;
    acceptcommand = [ 'g = get(gcbf, ''userdata'');' ...
        'TMPREJ = g.winrej;' ...
        'if g.children, delete(g.children); end;' ...
        'delete(gcbf);' ...
        tmpcom ...
        '; clear g;']; % quitting expression
    if ~isempty(g.command)
        u(12) = uicontrol('Parent',figh, ...
            'Units', 'normalized', ...
            'Position',posbut(12,:), ...
            'Tag','Accept',...
            'string',g.butlabel, 'callback', acceptcommand);
    end;
    u(13) = uicontrol('Parent',figh, ...
        'Units', 'normalized', ...
        'Position',posbut(13,:), ...
        'string',fastif(isempty(g.command),'CLOSE', 'CANCEL'), 'callback', ...
        [	'g = get(gcbf, ''userdata'');' ...
        'if g.children, delete(g.children); end;' ...
        'close(gcbf);'] );
    
    if ~isempty(g.events)
        u(17) = uicontrol('Parent',figh, ...
            'Units', 'normalized', ...
            'Position',posbut(17,:), ...
            'string', 'Event types', 'callback', 'eegplot(''drawlegend'', gcbf)');
    end;
    
    for i = 1: length(u) % Matlab 2014b compatibility
        if isprop(eval(['u(' num2str(i) ')']),'Style')
            set(u(i),'Units','Normalized');
        end
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up uimenus
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Figure Menu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    m(7) = uimenu('Parent',figh,'Label','Figure');
    m(8) = uimenu('Parent',m(7),'Label','Print');
    uimenu('Parent',m(7),'Label','Edit figure', 'Callback', 'eegplot(''noui'');');
    uimenu('Parent',m(7),'Label','Accept and close', 'Callback', acceptcommand );
    uimenu('Parent',m(7),'Label','Cancel and close', 'Callback','delete(gcbf)')
    
    % Portrait %%%%%%%%
    
    timestring = ['[OBJ1,FIG1] = gcbo;',...
        'PANT1 = get(OBJ1,''parent'');',...
        'OBJ2 = findobj(''tag'',''orient'',''parent'',PANT1);',...
        'set(OBJ2,''checked'',''off'');',...
        'set(OBJ1,''checked'',''on'');',...
        'set(FIG1,''PaperOrientation'',''portrait'');',...
        'clear OBJ1 FIG1 OBJ2 PANT1;'];
    
    uimenu('Parent',m(8),'Label','Portrait','checked',...
        'on','tag','orient','callback',timestring)
    
    % Landscape %%%%%%%
    timestring = ['[OBJ1,FIG1] = gcbo;',...
        'PANT1 = get(OBJ1,''parent'');',...
        'OBJ2 = findobj(''tag'',''orient'',''parent'',PANT1);',...
        'set(OBJ2,''checked'',''off'');',...
        'set(OBJ1,''checked'',''on'');',...
        'set(FIG1,''PaperOrientation'',''landscape'');',...
        'clear OBJ1 FIG1 OBJ2 PANT1;'];
    
    uimenu('Parent',m(8),'Label','Landscape','checked',...
        'off','tag','orient','callback',timestring)
    
    % Print command %%%%%%%
    uimenu('Parent',m(8),'Label','Print','tag','printcommand','callback',...
        ['RESULT = inputdlg2( { ''Command:'' }, ''Print'', 1,  { ''print -r72'' });' ...
        'if size( RESULT,1 ) ~= 0' ...
        '  eval ( RESULT{1} );' ...
        'end;' ...
        'clear RESULT;' ]);
    
    % Display Menu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    m(1) = uimenu('Parent',figh,...
        'Label','Display', 'tag', 'displaymenu');
    
    % window grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % userdata = 4 cells : display yes/no, color, electrode yes/no,
    %                      trial boundary adapt yes/no (1/0)
    m(11) = uimenu('Parent',m(1),'Label','Data select/mark', 'tag', 'displaywin', ...
        'userdata', { 1, [0.8 1 0.8], 0, fastif( g.trialstag(1) == -1, 0, 1)});
    
    uimenu('Parent',m(11),'Label','Hide marks','Callback', ...
        ['g = get(gcbf, ''userdata'');' ...
        'if ~g.winstatus' ...
        '  set(gcbo, ''label'', ''Hide marks'');' ...
        'else' ...
        '  set(gcbo, ''label'', ''Show marks'');' ...
        'end;' ...
        'g.winstatus = ~g.winstatus;' ...
        'set(gcbf, ''userdata'', g);' ...
        'eegplot(''drawb''); clear g;'] )
    
    % color %%%%%%%%%%%%%%%%%%%%%%%%%%
    if isunix % for some reasons, does not work under Windows
        uimenu('Parent',m(11),'Label','Choose color', 'Callback', ...
            [ 'g = get(gcbf, ''userdata'');' ...
            'g.wincolor = uisetcolor(g.wincolor);' ...
            'set(gcbf, ''userdata'', g ); ' ...
            'clear g;'] )
    end;
    
    % set channels
    %uimenu('Parent',m(11),'Label','Mark channels', 'enable', 'off', ...
    %'checked', 'off', 'Callback', ...
    %['g = get(gcbf, ''userdata'');' ...
    % 'g.setelectrode = ~g.setelectrode;' ...
    % 'set(gcbf, ''userdata'', g); ' ...
    % 'if ~g.setelectrode setgcbo, ''checked'', ''on''); ...
    % else set(gcbo, ''checked'', ''off''); end;'...
    % ' clear g;'] )
    
    % trials boundaries
    %uimenu('Parent',m(11),'Label','Trial boundaries', 'checked', fastif( g.trialstag(1) == -1, 'off', 'on'), 'Callback', ...
    %['hh = findobj(''tag'',''displaywin'',''parent'', findobj(''tag'',''displaymenu'',''parent'', gcbf ));' ...
    % 'hhdat = get(hh, ''userdata'');' ...
    % 'set(hh, ''userdata'', { hhdat{1},  hhdat{2}, hhdat{3}, ~hhdat{4}} ); ' ...
    %'if ~hhdat{4} set(gcbo, ''checked'', ''on''); else set(gcbo, ''checked'', ''off''); end;' ...
    %' clear hh hhdat;'] )
    
    % plot durations
    % --------------
    if g.ploteventdur & isfield(g.events, 'duration')
        disp(['Use menu "Display > Hide event duration" to hide colored regions ' ...
            'representing event duration']);
    end;
    if isfield(g.events, 'duration')
        uimenu('Parent',m(1),'Label',fastif(g.ploteventdur, 'Hide event duration', 'Plot event duration'),'Callback', ...
            ['g = get(gcbf, ''userdata'');' ...
            'if ~g.ploteventdur' ...
            '  set(gcbo, ''label'', ''Hide event duration'');' ...
            'else' ...
            '  set(gcbo, ''label'', ''Show event duration'');' ...
            'end;' ...
            'g.ploteventdur = ~g.ploteventdur;' ...
            'set(gcbf, ''userdata'', g);' ...
            'eegplot(''drawb''); clear g;'] )
    end;
    
    % X grid %%%%%%%%%%%%
    m(3) = uimenu('Parent',m(1),'Label','Grid');
    timestring = ['FIGH = gcbf;',...
        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
        'if size(get(AXESH,''xgrid''),2) == 2' ... %on
        '  set(AXESH,''xgrid'',''off'');',...
        '  set(gcbo,''label'',''X grid on'');',...
        'else' ...
        '  set(AXESH,''xgrid'',''on'');',...
        '  set(gcbo,''label'',''X grid off'');',...
        'end;' ...
        'clear FIGH AXESH;' ];
    uimenu('Parent',m(3),'Label',fastif(strcmp(g.xgrid, 'off'), ...
        'X grid on','X grid off'), 'Callback',timestring)
    
    % Y grid %%%%%%%%%%%%%
    timestring = ['FIGH = gcbf;',...
        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
        'if size(get(AXESH,''ygrid''),2) == 2' ... %on
        '  set(AXESH,''ygrid'',''off'');',...
        '  set(gcbo,''label'',''Y grid on'');',...
        'else' ...
        '  set(AXESH,''ygrid'',''on'');',...
        '  set(gcbo,''label'',''Y grid off'');',...
        'end;' ...
        'clear FIGH AXESH;' ];
    uimenu('Parent',m(3),'Label',fastif(strcmp(g.ygrid, 'off'), ...
        'Y grid on','Y grid off'), 'Callback',timestring)
    
    % Grid Style %%%%%%%%%
    m(5) = uimenu('Parent',m(3),'Label','Grid Style');
    timestring = ['FIGH = gcbf;',...
        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
        'set(AXESH,''gridlinestyle'',''--'');',...
        'clear FIGH AXESH;'];
    uimenu('Parent',m(5),'Label','- -','Callback',timestring)
    timestring = ['FIGH = gcbf;',...
        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
        'set(AXESH,''gridlinestyle'',''-.'');',...
        'clear FIGH AXESH;'];
    uimenu('Parent',m(5),'Label','_ .','Callback',timestring)
    timestring = ['FIGH = gcbf;',...
        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
        'set(AXESH,''gridlinestyle'','':'');',...
        'clear FIGH AXESH;'];
    uimenu('Parent',m(5),'Label','. .','Callback',timestring)
    timestring = ['FIGH = gcbf;',...
        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
        'set(AXESH,''gridlinestyle'',''-'');',...
        'clear FIGH AXESH;'];
    uimenu('Parent',m(5),'Label','__','Callback',timestring)
    
    % Submean menu %%%%%%%%%%%%%
    cb =       ['g = get(gcbf, ''userdata'');' ...
        'if strcmpi(g.submean, ''on''),' ...
        '  set(gcbo, ''label'', ''Remove DC offset'');' ...
        '  g.submean =''off'';' ...
        'else' ...
        '  set(gcbo, ''label'', ''Do not remove DC offset'');' ...
        '  g.submean =''on'';' ...
        'end;' ...
        'set(gcbf, ''userdata'', g);' ...
        'eegplot(''drawp'', 0); clear g;'];
    uimenu('Parent',m(1),'Label',fastif(strcmp(g.submean, 'on'), ...
        'Do not remove DC offset','Remove DC offset'), 'Callback',cb)
    
    % Scale Eye %%%%%%%%%
    timestring = ['[OBJ1,FIG1] = gcbo;',...
        'eegplot(''scaleeye'',OBJ1,FIG1);',...
        'clear OBJ1 FIG1;'];
    m(7) = uimenu('Parent',m(1),'Label','Show scale','Callback',timestring);
    
    % Title %%%%%%%%%%%%
    uimenu('Parent',m(1),'Label','Title','Callback','eegplot(''title'')')
    
    % Stack/Spread %%%%%%%%%%%%%%%
    cb =       ['g = get(gcbf, ''userdata'');' ...
        'hbutton = findobj(gcf, ''Tag'', ''Envelope'');' ...  % find button
        'if g.envelope == 0,' ...
        '  set(gcbo, ''label'', ''Spread channels'');' ...
        '  g.envelope = 1;' ...
        '  set(hbutton, ''String'', ''Spread'');' ...
        'else' ...
        '  set(gcbo, ''label'', ''Stack channels'');' ...
        '  g.envelope = 0;' ...
        '  set(hbutton, ''String'', ''Stack'');' ...
        'end;' ...
        'set(gcbf, ''userdata'', g);' ...
        'eegplot(''drawp'', 0); clear g;'];
    uimenu('Parent',m(1),'Label',fastif(g.envelope == 0, ...
        'Stack channels','Spread channels'), 'Callback',cb, 'Tag', 'Envelope_menu')
    
    % Normalize/denormalize %%%%%%%%%%%%%%%
    cb_normalize = ['g = get(gcbf,''userdata'');if g.normed, disp(''Denormalizing...''); else, disp(''Normalizing...''); end;'...
        'hbutton = findobj(gcf, ''Tag'', ''Norm'');' ...  % find button
        'ax1 = findobj(''tag'',''eegaxis'',''parent'',gcbf);' ...
        'data = get(ax1,''UserData'');' ...
        'if isempty(g.datastd), g.datastd = std(data(:,1:min(1000,g.frames),[],2)); end;'...
        'if g.normed, '...
        '  for i = 1:size(data,1), '...
        '    data(i,:,:) = data(i,:,:)*g.datastd(i);'...
        '    if ~isempty(g.data2), g.data2(i,:,:) = g.data2(i,:,:)*g.datastd(i);end;'...
        '  end;'...
        '  set(hbutton,''string'', ''Norm'');set(findobj(''tag'',''ESpacing'',''parent'',gcbf),''string'',num2str(g.oldspacing));' ...
        '  set(gcbo, ''label'', ''Normalize channels'');' ...
        'else, for i = 1:size(data,1),'...
        '    data(i,:,:) = data(i,:,:)/g.datastd(i);'...
        '    if ~isempty(g.data2), g.data2(i,:,:) = g.data2(i,:,:)/g.datastd(i);end;'...
        '  end;'...
        '  set(hbutton,''string'', ''Denorm'');'...
        '  set(gcbo, ''label'', ''Denormalize channels'');' ...
        '  g.oldspacing = g.spacing;set(findobj(''tag'',''ESpacing'',''parent'',gcbf),''string'',''5'');end;' ...
        'g.normed = 1 - g.normed;' ...
        'eegplot(''draws'',0);'...
        'set(gcbf,''userdata'',g);set(ax1,''UserData'',data);clear ax1 g data;' ...
        'eegplot(''drawp'',0);' ...
        'disp(''Done.'')'];
    uimenu('Parent',m(1),'Label',fastif(g.envelope == 0, ...
        'Normalize channels','Denormalize channels'), 'Callback',cb_normalize, 'Tag', 'Normalize_menu')
    
    
    % Settings Menu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m(2) = uimenu('Parent',figh,...
        'Label','Settings');
    
    % Window %%%%%%%%%%%%
    uimenu('Parent',m(2),'Label','Time range to display',...
        'Callback','eegplot(''window'')')
    
    % Electrode window %%%%%%%%
    uimenu('Parent',m(2),'Label','Number of channels to display',...
        'Callback','eegplot(''winelec'')')
    
    % Electrodes %%%%%%%%
    m(6) = uimenu('Parent',m(2),'Label','Channel labels');
    
    timestring = ['FIGH = gcbf;',...
        'AXESH = findobj(''tag'',''eegaxis'',''parent'',FIGH);',...
        'YTICK = get(AXESH,''YTick'');',...
        'YTICK = length(YTICK);',...
        'set(AXESH,''YTickLabel'',flipud(str2mat(num2str((1:YTICK-1)''),'' '')));',...
        'clear FIGH AXESH YTICK;'];
    uimenu('Parent',m(6),'Label','Show number','Callback',timestring)
    uimenu('Parent',m(6),'Label','Load .loc(s) file',...
        'Callback','eegplot(''loadelect'');')
    
    % Zooms %%%%%%%%
    zm = uimenu('Parent',m(2),'Label','Zoom off/on');
    if verLessThan('matlab','8.4.0')
        commandzoom = [ 'set(gcbf, ''WindowButtonDownFcn'', [ ''zoom(gcbf,''''down''''); eegplot(''''zoom'''', gcbf, 1);'' ]);' ...
            'tmpg = get(gcbf, ''userdata'');' ...
            'clear tmpg tmpstr;'];
    else
        % Temporary fix to avoid warning when setting a callback and the  mode is active
        % This is failing for us http://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
        commandzoom = [ 'wtemp = warning; warning off;set(gcbf, ''WindowButtonDownFcn'', [ ''zoom(gcbf); eegplot(''''zoom'''', gcbf, 1);'' ]);' ...
            'tmpg = get(gcbf, ''userdata'');' ...
            'warning(wtemp);'...
            'clear wtemp tmpg tmpstr; '];
    end
    
    %uimenu('Parent',zm,'Label','Zoom time', 'callback', ...
    %             [ 'zoom(gcbf, ''xon'');' commandzoom ]);
    %uimenu('Parent',zm,'Label','Zoom channels', 'callback', ...
    %             [ 'zoom(gcbf, ''yon'');' commandzoom ]);
    uimenu('Parent',zm,'Label','Zoom on', 'callback', commandzoom);
    uimenu('Parent',zm,'Label','Zoom off', 'separator', 'on', 'callback', ...
        ['zoom(gcbf, ''off''); tmpg = get(gcbf, ''userdata'');' ...
        'set(gcbf, ''windowbuttondownfcn'', tmpg.commandselect{1});' ...
        'set(gcbf, ''windowbuttonupfcn'', tmpg.commandselect{3});' ...
        'clear tmpg;' ]);
    
    uimenu('Parent',figh,'Label', 'Help', 'callback', 'pophelp(''eegplot'');');
    
    % Events %%%%%%%%
    zm = uimenu('Parent',m(2),'Label','Events');
    complotevent = [ 'tmpg = get(gcbf, ''userdata'');' ...
        'tmpg.plotevent = ''on'';' ...
        'set(gcbf, ''userdata'', tmpg); clear tmpg; eegplot(''drawp'', 0);'];
    comnoevent   = [ 'tmpg = get(gcbf, ''userdata'');' ...
        'tmpg.plotevent = ''off'';' ...
        'set(gcbf, ''userdata'', tmpg); clear tmpg; eegplot(''drawp'', 0);'];
    comeventmaxstring   = [ 'tmpg = get(gcbf, ''userdata'');' ...
        'tmpg.plotevent = ''on'';' ...
        'set(gcbf, ''userdata'', tmpg); clear tmpg; eegplot(''emaxstring'');']; % JavierLC
    comeventleg  = [ 'eegplot(''drawlegend'', gcbf);'];
    
    uimenu('Parent',zm,'Label','Events on'    , 'callback', complotevent, 'enable', fastif(isempty(g.events), 'off', 'on'));
    uimenu('Parent',zm,'Label','Events off'   , 'callback', comnoevent  , 'enable', fastif(isempty(g.events), 'off', 'on'));
    uimenu('Parent',zm,'Label','Events'' string length'   , 'callback', comeventmaxstring, 'enable', fastif(isempty(g.events), 'off', 'on')); % JavierLC
    uimenu('Parent',zm,'Label','Events'' legend', 'callback', comeventleg , 'enable', fastif(isempty(g.events), 'off', 'on'));
    
    
    % %%%%%%%%%%%%%%%%%
    % Set up autoselect
    % NOTE: commandselect{2} option has been moved to a
    %       subfunction to improve speed
    %%%%%%%%%%%%%%%%%%%
    g.commandselect{1} = [ 'if strcmp(get(gcbf, ''SelectionType''),''alt''),' g.ctrlselectcommand{1} ...
        'else '                                            g.selectcommand{1} 'end;' ];
    g.commandselect{3} = [ 'if strcmp(get(gcbf, ''SelectionType''),''alt''),' g.ctrlselectcommand{3} ...
        'else '                                            g.selectcommand{3} 'end;' ];
    
    set(figh, 'windowbuttondownfcn',   g.commandselect{1});
    set(figh, 'windowbuttonmotionfcn', {@defmotion,figh,ax0,ax1,u(10),u(11),u(9)});
    set(figh, 'windowbuttonupfcn',     g.commandselect{3});
    set(figh, 'WindowKeyPressFcn',     @eegplot_readkey);
    set(figh, 'interruptible', 'off');
    set(figh, 'busyaction', 'cancel');
    %  set(figh, 'windowbuttondownfcn', commandpush);
    %  set(figh, 'windowbuttonmotionfcn', commandmove);
    %  set(figh, 'windowbuttonupfcn', commandrelease);
    %  set(figh, 'interruptible', 'off');
    %  set(figh, 'busyaction', 'cancel');
    
    % prepare event array if any
    % --------------------------
    if ~isempty(g.events)
        if ~isfield(g.events, 'type') | ~isfield(g.events, 'latency'), g.events = []; end;
    end;
    
    if ~isempty(g.events)
        if isstr(g.events(1).type)
            [g.eventtypes tmpind indexcolor] = unique_bc({g.events.type}); % indexcolor countinas the event type
        else [g.eventtypes tmpind indexcolor] = unique_bc([ g.events.type ]);
        end;
        g.eventcolors     = { 'r', [0 0.8 0], 'm', 'c', 'k', 'b', [0 0.8 0] };
        g.eventstyle      = { '-' '-' '-'  '-'  '-' '-' '-' '--' '--' '--'  '--' '--' '--' '--'};
        g.eventwidths     = [ 2.5 1 ];
        g.eventtypecolors = g.eventcolors(mod([1:length(g.eventtypes)]-1 ,length(g.eventcolors))+1);
        g.eventcolors     = g.eventcolors(mod(indexcolor-1               ,length(g.eventcolors))+1);
        g.eventtypestyle  = g.eventstyle (mod([1:length(g.eventtypes)]-1 ,length(g.eventstyle))+1);
        g.eventstyle      = g.eventstyle (mod(indexcolor-1               ,length(g.eventstyle))+1);
        
        % for width, only boundary events have width 2 (for the line)
        % -----------------------------------------------------------
        indexwidth = ones(1,length(g.eventtypes))*2;
        if iscell(g.eventtypes)
            for index = 1:length(g.eventtypes)
                if strcmpi(g.eventtypes{index}, 'boundary'), indexwidth(index) = 1; end;
            end;
        end;
        g.eventtypewidths = g.eventwidths (mod(indexwidth([1:length(g.eventtypes)])-1 ,length(g.eventwidths))+1);
        g.eventwidths     = g.eventwidths (mod(indexwidth(indexcolor)-1               ,length(g.eventwidths))+1);
        
        % latency and duration of events
        % ------------------------------
        g.eventlatencies  = [ g.events.latency ]+1;
        if isfield(g.events, 'duration')
            durations = { g.events.duration };
            durations(cellfun(@isempty, durations)) = { NaN };
            g.eventlatencyend   = g.eventlatencies + [durations{:}]+1;
        else g.eventlatencyend   = [];
        end;
        g.plotevent       = 'on';
    end;
    if isempty(g.events)
        g.plotevent      = 'off';
    end;
    
    set(figh, 'userdata', g);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot EEG Data
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(ax1)
    hold on
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot Spacing I
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    YLim = get(ax1,'Ylim');
    A = DEFAULT_AXES_POSITION;
    axes('Position',[A(1)+A(3) A(2) 1-A(1)-A(3) A(4)],'Visible','off','Ylim',YLim,'tag','eyeaxes')
    axis manual
    if strcmp(SPACING_EYE,'on'),  set(m(7),'checked','on')
    else set(m(7),'checked','off');
    end
    eegplot('scaleeye', [], gcf);
    if strcmp(lower(g.scale), 'off')
        eegplot('scaleeye', 'off', gcf);
    end;
    
    eegplot('drawp', 0);
    %   eegplot('drawp', 0);
    if g.dispchans ~= g.chans
        eegplot('zoom', gcf);
    end;
    eegplot('scaleeye', [], gcf);
    
    h = findobj(gcf, 'style', 'pushbutton');
    set(h, 'backgroundcolor', BUTTON_COLOR);
    h = findobj(gcf, 'tag', 'eegslider');
    set(h, 'backgroundcolor', BUTTON_COLOR);
    set(figh, 'visible', 'on');
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End Main Function
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else
    try, p1 = varargin{1}; p2 = varargin{2}; p3 = varargin{3}; catch, end;
    switch data
        case 'drawp' % Redraw EEG and change position
            
            % this test help to couple eegplot windows
            if exist('p3', 'var')
                figh = p3;
                figure(p3);
            else
                figh = gcf;                          % figure handle
            end;
            
            if strcmp(get(figh,'tag'),'dialog')
                figh = get(figh,'UserData');
            end
            ax0 = findobj('tag','backeeg','parent',figh); % axes handle
            ax1 = findobj('tag','eegaxis','parent',figh); % axes handle
            g = get(figh,'UserData');
            data = get(ax1,'UserData');
            ESpacing = findobj('tag','ESpacing','parent',figh);   % ui handle
            EPosition = findobj('tag','EPosition','parent',figh); % ui handle
            if ~isempty(EPosition) && ~isempty(ESpacing)
                if g.trialstag(1) == -1
                    g.time    = str2num(get(EPosition,'string'));
                else
                    g.time    = str2num(get(EPosition,'string'));
                    g.time    = g.time - 1;
                end;
                g.spacing = str2num(get(ESpacing,'string'));
            end;
            
            if p1 == 1
                g.time = g.time-g.winlength;     % << subtract one window length
            elseif p1 == 2
                g.time = g.time-fastif(g.winlength>=1, 1, g.winlength/5);             % < subtract one second
            elseif p1 == 3
                g.time = g.time+fastif(g.winlength>=1, 1, g.winlength/5);             % > add one second
            elseif p1 == 4
                g.time = g.time+g.winlength;     % >> add one window length
            end
            
            if g.trialstag ~= -1 % time in second or in trials
                multiplier = g.trialstag;
            else
                multiplier = g.srate;
            end;
            
            % Update edit box
            % ---------------
            g.time = max(0,min(g.time,ceil((g.frames-1)/multiplier)-g.winlength));
            if g.trialstag(1) == -1
                set(EPosition,'string',num2str(g.time));
            else
                set(EPosition,'string',num2str(g.time+1));
            end;
            set(figh, 'userdata', g);
            
            lowlim = round(g.time*multiplier+1);
            highlim = round(min((g.time+g.winlength)*multiplier+2,g.frames));
            
            % Plot data and update axes
            % -------------------------
            if ~isempty(g.data2)
                switch lower(g.submean) % subtract the mean ?
                    case 'on',
                        meandata = mean(g.data2(:,lowlim:highlim)');
                        if any(isnan(meandata))
                            meandata = nan_mean(g.data2(:,lowlim:highlim)');
                        end;
                    otherwise, meandata = zeros(1,g.chans);
                end;
            else
                switch lower(g.submean) % subtract the mean ?
                    case 'on',
                        meandata = mean(data(:,lowlim:highlim)');
                        if any(isnan(meandata))
                            meandata = nan_mean(data(:,lowlim:highlim)');
                        end;
                    otherwise, meandata = zeros(1,g.chans);
                end;
            end;
            if strcmpi(g.plotdata2, 'off')
                axes(ax1)
                cla
            end;
            
            oldspacing = g.spacing;
            if g.envelope
                g.spacing = 0;
            end
            % plot data
            % ---------
            axes(ax1)
            hold on
            
            % plot channels whose "badchan" field is set to 1.
            % Bad channels are plotted first so that they appear behind the good
            % channels in the eegplot figure window.
            for i = 1:g.chans
                if strcmpi(g.plotdata2, 'on')
                    tmpcolor = [ 1 0 0 ];
                else tmpcolor = g.color{mod(i-1,length(g.color))+1};
                end;
                
                if isfield(g, 'eloc_file') & ...
                        isfield(g.eloc_file, 'badchan') & ...
                        g.eloc_file(g.chans-i+1).badchan;
                    tmpcolor = [ .85 .85 .85 ];
                    plot(data(g.chans-i+1,lowlim:highlim) -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), ...
                        'color', tmpcolor, 'clipping','on')
                    plot(1,mean(data(g.chans-i+1,lowlim:highlim) -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing),2),'<r','MarkerFaceColor','r','MarkerSize',6);
                end
                
            end
            
            % plot good channels on top of bad channels (if g.eloc_file(i).badchan = 0... or there is no bad channel information)
            for i = 1:g.chans
                if strcmpi(g.plotdata2, 'on')
                    tmpcolor = [ 1 0 0 ];
                else tmpcolor = g.color{mod(g.chans-i,length(g.color))+1};
                end;
                
                %        keyboard;
                if (isfield(g, 'eloc_file') & ...
                        isfield(g.eloc_file, 'badchan') & ...
                        ~g.eloc_file(g.chans-i+1).badchan) | ...
                        (~isfield(g, 'eloc_file')) | ...
                        (~isfield(g.eloc_file, 'badchan'));
                    plot(data(g.chans-i+1,lowlim:highlim) -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), ...
                        'color', tmpcolor, 'clipping','on')
                end
                
            end
            
            % draw selected channels
            % ------------------------
            if ~isempty(g.winrej) & size(g.winrej,2) > 2
                for tpmi = 1:size(g.winrej,1) % scan rows
                    if (g.winrej(tpmi,1) >= lowlim & g.winrej(tpmi,1) <= highlim) | ...
                            (g.winrej(tpmi,2) >= lowlim & g.winrej(tpmi,2) <= highlim)
                        abscmin = max(1,round(g.winrej(tpmi,1)-lowlim));
                        abscmax = round(g.winrej(tpmi,2)-lowlim);
                        maxXlim = get(gca, 'xlim');
                        abscmax = min(abscmax, round(maxXlim(2)-1));
                        for i = 1:g.chans
                            if g.winrej(tpmi,g.chans-i+1+5)
                                plot(abscmin+1:abscmax+1,data(g.chans-i+1,abscmin+lowlim:abscmax+lowlim) ...
                                    -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), 'color','r','clipping','on')
                            end;
                        end
                    end;
                end;
            end;
            g.spacing = oldspacing;
            set(ax1, 'Xlim',[1 g.winlength*multiplier+1],...
                'XTick',[1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1]);
            %          if g.isfreq % Ramon
            %              set(ax1, 'XTickLabel', num2str((g.freqs(1):DEFAULT_GRID_SPACING:g.freqs(end))'));
            %          else
            
            if isfield(g,'limits')
                g.time = g.limits(1);
                set(gcf,'UserData',g);
            else
                g.time = 0;
            end
            set(ax1, 'XTickLabel', num2str((g.time:DEFAULT_GRID_SPACING:g.time+g.winlength)'));
            %          end
            
            % ordinates: even if all elec are plotted, some may be hidden
            set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );
            
            if g.children ~= 0
                if ~exist('p2', 'var')
                    p2 =[];
                end;
                eegplot( 'drawp', p1, p2, g.children);
                figure(figh);
            end;
            
            % draw second data if necessary
            if ~isempty(g.data2)
                tmpdata = data;
                set(ax1, 'userdata', g.data2);
                g.data2 = [];
                g.plotdata2 = 'on';
                set(figh, 'userdata', g);
                eegplot('drawp', 0);
                g.plotdata2 = 'off';
                g.data2 = get(ax1, 'userdata');
                set(ax1, 'userdata', tmpdata);
                set(figh, 'userdata', g);
            else
                eegplot('drawb');
            end;
            
        case 'drawb' % Draw background ******************************************************
            % Redraw EEG and change position
            
            ax0 = findobj('tag','backeeg','parent',gcf); % axes handle
            ax1 = findobj('tag','eegaxis','parent',gcf); % axes handle
            
            g = get(gcf,'UserData');  % Data (Note: this could also be global)
            
            % Plot data and update axes
            axes(ax0);
            cla;
            hold on;
            % plot rejected windows
            if g.trialstag ~= -1
                multiplier = g.trialstag;
            else
                multiplier = g.srate;
            end;
            
            % draw rejection windows
            % ----------------------
            lowlim = round(g.time*multiplier+1);
            highlim = round(min((g.time+g.winlength)*multiplier+1));
            displaymenu = findobj('tag','displaymenu','parent',gcf);
            if ~isempty(g.winrej) & g.winstatus
                if g.trialstag ~= -1 % epoched data
                    indices = find((g.winrej(:,1)' >= lowlim & g.winrej(:,1)' <= highlim) | ...
                        (g.winrej(:,2)' >= lowlim & g.winrej(:,2)' <= highlim));
                    if ~isempty(indices)
                        tmpwins1 = g.winrej(indices,1)';
                        tmpwins2 = g.winrej(indices,2)';
                        if size(g.winrej,2) > 2
                            tmpcols  = g.winrej(indices,3:5);
                        else tmpcols  = g.wincolor;
                        end;
                        try, eval('[cumul indicescount] = histc(tmpwins1, (min(tmpwins1)-1):g.trialstag:max(tmpwins2));');
                        catch, [cumul indicescount] = myhistc(tmpwins1, (min(tmpwins1)-1):g.trialstag:max(tmpwins2));
                        end;
                        count = zeros(size(cumul));
                        %if ~isempty(find(cumul > 1)), find(cumul > 1), end;
                        for tmpi = 1:length(tmpwins1)
                            poscumul = indicescount(tmpi);
                            heightbeg = count(poscumul)/cumul(poscumul);
                            heightend = heightbeg + 1/cumul(poscumul);
                            count(poscumul) = count(poscumul)+1;
                            h = patch([tmpwins1(tmpi)-lowlim tmpwins2(tmpi)-lowlim ...
                                tmpwins2(tmpi)-lowlim tmpwins1(tmpi)-lowlim], ...
                                [heightbeg heightbeg heightend heightend], ...
                                tmpcols(tmpi,:));  % this argument is color
                            set(h, 'EdgeColor', get(h, 'facecolor'))
                        end;
                    end;
                else
                    event2plot1 = find ( g.winrej(:,1) >= lowlim & g.winrej(:,1) <= highlim );
                    event2plot2 = find ( g.winrej(:,2) >= lowlim & g.winrej(:,2) <= highlim );
                    event2plot3 = find ( g.winrej(:,1) <  lowlim & g.winrej(:,2) >  highlim );
                    event2plot  = union_bc(union(event2plot1, event2plot2), event2plot3);
                    
                    for tpmi = event2plot(:)'
                        if size(g.winrej,2) > 2
                            tmpcols  = g.winrej(tpmi,3:5);
                        else tmpcols  = g.wincolor;
                        end;
                        h = patch([g.winrej(tpmi,1)-lowlim g.winrej(tpmi,2)-lowlim ...
                            g.winrej(tpmi,2)-lowlim g.winrej(tpmi,1)-lowlim], ...
                            [0 0 1 1], tmpcols);
                        set(h, 'EdgeColor', get(h, 'facecolor'))
                    end;
                end;
            end;
            
            % plot tags
            % ---------
            %if trialtag(1) ~= -1 & displaystatus % put tags at arbitrary places
            % 	for tmptag = trialtag
            %		if tmptag >= lowlim & tmptag <= highlim
            %			plot([tmptag-lowlim tmptag-lowlim], [0 1], 'b--');
            %		end;
            %	end;
            %end;
            
            % draw events if any
            % ------------------
            if strcmpi(g.plotevent, 'on')
                
                % JavierLC ###############################
                MAXEVENTSTRING = g.maxeventstring;
                if MAXEVENTSTRING<0
                    MAXEVENTSTRING = 0;
                elseif MAXEVENTSTRING>75
                    MAXEVENTSTRING=75;
                end
                AXES_POSITION = [0.0964286 0.15 0.842 0.75-(MAXEVENTSTRING-5)/100];
                % JavierLC ###############################
                
                % find event to plot
                % ------------------
                event2plot    = find ( g.eventlatencies >=lowlim & g.eventlatencies <= highlim );
                if ~isempty(g.eventlatencyend)
                    event2plot2 = find ( g.eventlatencyend >= lowlim & g.eventlatencyend <= highlim );
                    event2plot3 = find ( g.eventlatencies  <  lowlim & g.eventlatencyend >  highlim );
                    event2plot  = union_bc(union(event2plot, event2plot2), event2plot3);
                end;
                for index = 1:length(event2plot)
                    %Just repeat for the first one
                    if index == 1
                        EVENTFONT = ' \fontsize{10} ';
                        ylims=ylim;
                    end
                    
                    % draw latency line
                    % -----------------
                    tmplat = g.eventlatencies(event2plot(index))-lowlim-1;
                    tmph   = plot([ tmplat tmplat ], ylims, 'color', g.eventcolors{ event2plot(index) }, ...
                        'linestyle', g.eventstyle { event2plot(index) }, ...
                        'linewidth', g.eventwidths( event2plot(index) ) );
                    
                    % schtefan: add Event types text above event latency line
                    % -------------------------------------------------------
                    %             EVENTFONT = ' \fontsize{10} ';
                    %             ylims=ylim;
                    evntxt = strrep(num2str(g.events(event2plot(index)).type),'_','-');
                    if length(evntxt)>MAXEVENTSTRING, evntxt = [ evntxt(1:MAXEVENTSTRING-1) '...' ]; end; % truncate
                    try,
                        tmph2 = text([tmplat], ylims(2)-0.005, [EVENTFONT evntxt], ...
                            'color', g.eventcolors{ event2plot(index) }, ...
                            'horizontalalignment', 'left',...
                            'rotation',90);
                    catch, end;
                    
                    % draw duration is not 0
                    % ----------------------
                    if g.ploteventdur & ~isempty(g.eventlatencyend) ...
                            & g.eventwidths( event2plot(index) ) ~= 2.5 % do not plot length of boundary events
                        tmplatend = g.eventlatencyend(event2plot(index))-lowlim-1;
                        if tmplatend ~= 0,
                            tmplim = ylims;
                            tmpcol = g.eventcolors{ event2plot(index) };
                            h = patch([ tmplat tmplatend tmplatend tmplat ], ...
                                [ tmplim(1) tmplim(1) tmplim(2) tmplim(2) ], ...
                                tmpcol );  % this argument is color
                            set(h, 'EdgeColor', 'none')
                        end;
                    end;
                end;
            else % JavierLC
                MAXEVENTSTRING = 10; % default
                AXES_POSITION = [0.0964286 0.15 0.842 0.75-(MAXEVENTSTRING-5)/100];
            end;
            
            if g.trialstag(1) ~= -1
                
                % plot trial limits
                % -----------------
                tmptag = [lowlim:highlim];
                tmpind = find(mod(tmptag-1, g.trialstag) == 0);
                for index = tmpind
                    plot([tmptag(index)-lowlim-1 tmptag(index)-lowlim-1], [0 1], 'b--');
                end;
                alltag = tmptag(tmpind);
                
                % compute Xticks
                % --------------
                tagnum = (alltag-1)/g.trialstag+1;
                set(ax0,'XTickLabel', tagnum,'YTickLabel', [],...
                    'Xlim',[0 g.winlength*multiplier],...
                    'XTick',alltag-lowlim+g.trialstag/2, 'YTick',[], 'tag','backeeg');
                
                axes(ax1);
                tagpos  = [];
                tagtext = [];
                if ~isempty(alltag)
                    alltag = [alltag(1)-g.trialstag alltag alltag(end)+g.trialstag]; % add border trial limits
                else
                    alltag = [ floor(lowlim/g.trialstag)*g.trialstag ceil(highlim/g.trialstag)*g.trialstag ]+1;
                end;
                
                nbdiv = 20/g.winlength; % approximative number of divisions
                divpossible = [ 100000./[1 2 4 5] 10000./[1 2 4 5] 1000./[1 2 4 5] 100./[1 2 4 5 10 20]]; % possible increments
                [tmp indexdiv] = min(abs(nbdiv*divpossible-(g.limits(2)-g.limits(1)))); % closest possible increment
                incrementpoint = divpossible(indexdiv)/1000*g.srate;
                
                % tag zero below is an offset used to be sure that 0 is included
                % in the absicia of the data epochs
                if g.limits(2) < 0, tagzerooffset  = (g.limits(2)-g.limits(1))/1000*g.srate+1;
                else                tagzerooffset  = -g.limits(1)/1000*g.srate;
                end;
                if tagzerooffset < 0, tagzerooffset = 0; end;
                
                for i=1:length(alltag)-1
                    if ~isempty(tagpos) & tagpos(end)-alltag(i)<2*incrementpoint/3
                        tagpos  = tagpos(1:end-1);
                    end;
                    if ~isempty(g.freqlimits)
                        tagpos  = [ tagpos linspace(alltag(i),alltag(i+1)-1, nbdiv) ];
                    else
                        if tagzerooffset ~= 0
                            tmptagpos = [alltag(i)+tagzerooffset:-incrementpoint:alltag(i)];
                        else
                            tmptagpos = [];
                        end;
                        tagpos  = [ tagpos [tmptagpos(end:-1:2) alltag(i)+tagzerooffset:incrementpoint:(alltag(i+1)-1)]];
                    end;
                end;
                
                % find corresponding epochs
                % -------------------------
                if ~g.isfreq
                    tmplimit = g.limits;
                    tpmorder = 1E-3;
                else
                    tmplimit = g.freqlimits;
                    tpmorder = 1;
                end
                tagtext = eeg_point2lat(tagpos, floor((tagpos)/g.trialstag)+1, g.srate, tmplimit,tpmorder);
                set(ax1,'XTickLabel', tagtext,'XTick', tagpos-lowlim);
            else
                set(ax0,'XTickLabel', [],'YTickLabel', [],...
                    'Xlim',[0 g.winlength*multiplier],...
                    'XTick',[], 'YTick',[], 'tag','backeeg');
                
                axes(ax1);
                if g.isfreq
                    set(ax1, 'XTickLabel', num2str((g.freqs(1):DEFAULT_GRID_SPACING:g.freqs(end))'),...
                        'XTick',[1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1]);
                else
                    set(ax1,'XTickLabel', num2str((g.time:DEFAULT_GRID_SPACING:g.time+g.winlength)'),...
                        'XTick',[1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1]);
                end
                
                set(ax1, 'Position', AXES_POSITION) % JavierLC
                set(ax0, 'Position', AXES_POSITION) % JavierLC
            end;
            
            % ordinates: even if all elec are plotted, some may be hidden
            set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );
            
            axes(ax1)
            
        case 'draws'
            % Redraw EEG and change scale
            
            ax1 = findobj('tag','eegaxis','parent',gcf);         % axes handle
            g = get(gcf,'UserData');
            data = get(ax1, 'userdata');
            ESpacing = findobj('tag','ESpacing','parent',gcf);   % ui handle
            EPosition = findobj('tag','EPosition','parent',gcf); % ui handle
            if g.trialstag(1) == -1
                g.time    = str2num(get(EPosition,'string'));
            else
                g.time    = str2num(get(EPosition,'string'))-1;
            end;
            g.spacing = str2num(get(ESpacing,'string'));
            
            orgspacing= g.spacing;
            if p1 == 1
                g.spacing= g.spacing+ 0.1*orgspacing; % increase g.spacing(5%)
            elseif p1 == 2
                g.spacing= max(0,g.spacing-0.1*orgspacing); % decrease g.spacing(5%)
            end
            if round(g.spacing*100) == 0
                maxindex = min(10000, g.frames);
                g.spacing = 0.01*max(max(data(:,1:maxindex),[],2),[],1)-min(min(data(:,1:maxindex),[],2),[],1);  % Set g.spacingto max/min data
            end;
            
            % update edit box
            % ---------------
            set(ESpacing,'string',num2str(g.spacing,4))
            set(gcf, 'userdata', g);
            eegplot('drawp', 0);
            set(ax1,'YLim',[0 (g.chans+1)*g.spacing],'YTick',[0:g.spacing:g.chans*g.spacing])
            set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );
            
            % update scaling eye (I) if it exists
            % -----------------------------------
            eyeaxes = findobj('tag','eyeaxes','parent',gcf);
            if ~isempty(eyeaxes)
                eyetext = findobj('type','text','parent',eyeaxes,'tag','thescalenum');
                set(eyetext,'string',num2str(g.spacing,4))
            end
            
            return;
            
        case 'window'  % change window size
            % get new window length with dialog box
            % -------------------------------------
            g = get(gcf,'UserData');
            result       = inputdlg2( { fastif(g.trialstag==-1,'New window length (s):', 'Number of epoch(s):') }, 'Change window length', 1,  { num2str(g.winlength) });
            if size(result,1) == 0 return; end;
            
            g.winlength = eval(result{1});
            set(gcf, 'UserData', g);
            eegplot('drawp',0);
            return;
            
        case 'winelec'  % change channel window size
            % get new window length with dialog box
            % -------------------------------------
            fig = gcf;
            g = get(gcf,'UserData');
            result = inputdlg2( ...
                { 'Number of channels to display:' } , 'Change number of channels to display', 1,  { num2str(g.dispchans) });
            if size(result,1) == 0 return; end;
            
            g.dispchans = eval(result{1});
            if g.dispchans<0 | g.dispchans>g.chans
                g.dispchans =g.chans;
            end;
            set(gcf, 'UserData', g);
            eegplot('updateslider', fig);
            eegplot('drawp',0);
            eegplot('scaleeye', [], fig);
            return;
            
        case 'emaxstring'  % change events' string length  ;  JavierLC
            % get dialog box
            % -------------------------------------
            g = get(gcf,'UserData');
            result = inputdlg2({ 'Max events'' string length:' } , 'Change events'' string length to display', 1,  { num2str(g.maxeventstring) });
            if size(result,1) == 0 return; end;
            g.maxeventstring = eval(result{1});
            set(gcf, 'UserData', g);
            eegplot('drawb');
            return;
            
        case 'loadelect' % load channels
            [inputname,inputpath] = uigetfile('*','Channel locations file');
            if inputname == 0 return; end;
            if ~exist([ inputpath inputname ])
                error('no such file');
            end;
            
            AXH0 = findobj('tag','eegaxis','parent',gcf);
            eegplot('setelect',[ inputpath inputname ],AXH0);
            return;
            
        case 'setelect'
            % Set channels
            eloc_file = p1;
            axeshand = p2;
            outvar1 = 1;
            if isempty(eloc_file)
                outvar1 = 0;
                return
            end
            
            tmplocs = readlocs(eloc_file);
            YLabels = { tmplocs.labels };
            YLabels = strvcat(YLabels);
            
            YLabels = flipud(str2mat(YLabels,' '));
            set(axeshand,'YTickLabel',YLabels)
            
        case 'title'
            % Get new title
            h = findobj('tag', 'eegplottitle');
            
            if ~isempty(h)
                result       = inputdlg2( { 'New title:' }, 'Change title', 1,  { get(h(1), 'string') });
                if ~isempty(result), set(h, 'string', result{1}); end;
            else
                result       = inputdlg2( { 'New title:' }, 'Change title', 1,  { '' });
                if ~isempty(result), h = textsc(result{1}, 'title'); set(h, 'tag', 'eegplottitle');end;
            end;
            
            return;
            
        case 'scaleeye'
            % Turn scale I on/off
            obj = p1;
            figh = p2;
            g = get(figh,'UserData');
            % figh = get(obj,'Parent');
            
            if ~isempty(obj)
                eyeaxes = findobj('tag','eyeaxes','parent',figh);
                children = get(eyeaxes,'children');
                if isstr(obj)
                    if strcmp(obj, 'off')
                        set(children, 'visible', 'off');
                        set(eyeaxes, 'visible', 'off');
                        return;
                    else
                        set(children, 'visible', 'on');
                        set(eyeaxes, 'visible', 'on');
                    end;
                else
                    toggle = get(obj,'checked');
                    if strcmp(toggle,'on')
                        set(children, 'visible', 'off');
                        set(eyeaxes, 'visible', 'off');
                        set(obj,'checked','off');
                        return;
                    else
                        set(children, 'visible', 'on');
                        set(eyeaxes, 'visible', 'on');
                        set(obj,'checked','on');
                    end;
                end;
            end;
            
            eyeaxes = findobj('tag','eyeaxes','parent',figh);
            ax1 = findobj('tag','eegaxis','parent',gcf); % axes handle
            YLim = double(get(ax1, 'ylim'));
            
            ESpacing = findobj('tag','ESpacing','parent',figh);
            g.spacing= str2num(get(ESpacing,'string'));
            
            axes(eyeaxes); cla; axis off;
            set(eyeaxes, 'ylim', YLim);
            
            Xl = double([.35 .65; .5 .5; .35 .65]);
            Yl = double([ g.spacing g.spacing; g.spacing 0; 0 0] + YLim(1));
            plot(Xl(1,:),Yl(1,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline'); hold on;
            plot(Xl(2,:),Yl(2,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline');
            plot(Xl(3,:),Yl(3,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline');
            text(.5,(YLim(2)-YLim(1))/23+Yl(1),num2str(g.spacing,4),...
                'HorizontalAlignment','center','FontSize',10,...
                'tag','thescalenum')
            text(Xl(2)+.1,Yl(1),'+','HorizontalAlignment','left',...
                'verticalalignment','middle', 'tag', 'thescale')
            text(Xl(2)+.1,Yl(4),'-','HorizontalAlignment','left',...
                'verticalalignment','middle', 'tag', 'thescale')
            if ~isempty(SPACING_UNITS_STRING)
                text(.5,-YLim(2)/23+Yl(4),SPACING_UNITS_STRING,...
                    'HorizontalAlignment','center','FontSize',10, 'tag', 'thescale')
            end
            text(.5,(YLim(2)-YLim(1))/10+Yl(1),'Scale',...
                'HorizontalAlignment','center','FontSize',10, 'tag', 'thescale')
            set(eyeaxes, 'tag', 'eyeaxes');
            
        case 'noui'
            if ~isempty(varargin)
                eegplot( varargin{:} ); fig = gcf;
            else
                fig = findobj('tag', 'EEGPLOT');
            end;
            set(fig, 'menubar', 'figure');
            
            % find button and text
            obj = findobj(fig, 'style', 'pushbutton'); delete(obj);
            obj = findobj(fig, 'style', 'edit'); delete(obj);
            obj = findobj(fig, 'style', 'text');
            %objscale = findobj(obj, 'tag', 'thescale');
            %delete(setdiff(obj, objscale));
            obj = findobj(fig, 'tag', 'Eelec');delete(obj);
            obj = findobj(fig, 'tag', 'Etime');delete(obj);
            obj = findobj(fig, 'tag', 'Evalue');delete(obj);
            obj = findobj(fig, 'tag', 'Eelecname');delete(obj);
            obj = findobj(fig, 'tag', 'Etimename');delete(obj);
            obj = findobj(fig, 'tag', 'Evaluename');delete(obj);
            obj = findobj(fig, 'type', 'uimenu');delete(obj);
            
        case 'zoom' % if zoom
            fig = varargin{1};
            ax1 = findobj('tag','eegaxis','parent',fig);
            ax2 = findobj('tag','backeeg','parent',fig);
            tmpxlim  = get(ax1, 'xlim');
            tmpylim  = get(ax1, 'ylim');
            tmpxlim2 = get(ax2, 'xlim');
            set(ax2, 'xlim', get(ax1, 'xlim'));
            g = get(fig,'UserData');
            
            % deal with abscissa
            % ------------------
            if g.trialstag ~= -1
                Eposition = str2num(get(findobj('tag','EPosition','parent',fig), 'string'));
                g.winlength = (tmpxlim(2) - tmpxlim(1))/g.trialstag;
                Eposition = Eposition + (tmpxlim(1) - tmpxlim2(1)-1)/g.trialstag;
                Eposition = round(Eposition*1000)/1000;
                set(findobj('tag','EPosition','parent',fig), 'string', num2str(Eposition));
            else
                Eposition = str2num(get(findobj('tag','EPosition','parent',fig), 'string'))-1;
                g.winlength = (tmpxlim(2) - tmpxlim(1))/g.srate;
                Eposition = Eposition + (tmpxlim(1) - tmpxlim2(1)-1)/g.srate;
                Eposition = round(Eposition*1000)/1000;
                set(findobj('tag','EPosition','parent',fig), 'string', num2str(Eposition+1));
            end;
            
            % deal with ordinate
            % ------------------
            g.elecoffset = tmpylim(1)/g.spacing;
            g.dispchans  = round(1000*(tmpylim(2)-tmpylim(1))/g.spacing)/1000;
            
            set(fig,'UserData', g);
            eegplot('updateslider', fig);
            eegplot('drawp', 0);
            eegplot('scaleeye', [], fig);
            
            % reactivate zoom if 3 arguments
            % ------------------------------
            if exist('p2', 'var') == 1
                if verLessThan('matlab','8.4.0')
                    set(gcbf, 'windowbuttondownfcn', [ 'zoom(gcbf,''down''); eegplot(''zoom'', gcbf, 1);' ]);
                else
                    % This is failing for us: http://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
                    %               hManager = uigetmodemanager(gcbf);
                    %               [hManager.WindowListenerHandles.Enabled] = deal(false);
                    
                    % Temporary fix
                    wtemp = warning; warning off;
                    set(gcbf, 'WindowButtonDownFcn', [ 'zoom(gcbf); eegplot(''zoom'', gcbf, 1);' ]);
                    warning(wtemp);
                end
            end;
            
        case 'updateslider' % if zoom
            fig = varargin{1};
            g = get(fig,'UserData');
            sliider = findobj('tag','eegslider','parent',fig);
            if g.elecoffset < 0
                g.elecoffset = 0;
            end;
            if g.dispchans >= g.chans
                g.dispchans = g.chans;
                g.elecoffset = 0;
                set(sliider, 'visible', 'off');
            else
                set(sliider, 'visible', 'on');
                set(sliider, 'value', g.elecoffset/g.chans, ...
                    'sliderstep', [1/(g.chans-g.dispchans) g.dispchans/(g.chans-g.dispchans)]);
                %'sliderstep', [1/(g.chans-1) g.dispchans/(g.chans-1)]);
            end;
            if g.elecoffset < 0
                g.elecoffset = 0;
            end;
            if g.elecoffset > g.chans-g.dispchans
                g.elecoffset = g.chans-g.dispchans;
            end;
            set(fig,'UserData', g);
            eegplot('scaleeye', [], fig);
            
        case 'drawlegend'
            fig = varargin{1};
            g = get(fig,'UserData');
            
            if ~isempty(g.events) % draw vertical colored lines for events, add event name text above
                nleg = length(g.eventtypes);
                fig2 = figure('numbertitle', 'off', 'name', '', 'visible', 'off', 'menubar', 'none', 'color', DEFAULT_FIG_COLOR);
                pos = get(fig2, 'position');
                set(fig2, 'position', [ pos(1) pos(2) 130 14*nleg+20]);
                
                for index = 1:nleg
                    plot([10 30], [(index-0.5) * 10 (index-0.5) * 10], 'color', g.eventtypecolors{index}, 'linestyle', ...
                        g.eventtypestyle{ index }, 'linewidth', g.eventtypewidths( index )); hold on;
                    if iscell(g.eventtypes)
                        th=text(35, (index-0.5)*10, g.eventtypes{index}, ...
                            'color', g.eventtypecolors{index});
                    else
                        th=text(35, (index-0.5)*10, num2str(g.eventtypes(index)), ...
                            'color', g.eventtypecolors{index});
                    end;
                end;
                xlim([0 130]);
                ylim([0 nleg*10]);
                axis off;
                set(fig2, 'visible', 'on');
            end;
            
            
            % motion button: move windows or display current position (channel, g.time and activation)
            % ----------------------------------------------------------------------------------------
            % case moved as subfunction
            % add topoplot
            % ------------
        case 'topoplot'
            fig = varargin{1};
            g = get(fig,'UserData');
            if ~isstruct(g.eloc_file) || ~isfield(g.eloc_file, 'theta') || isempty( [ g.eloc_file.theta ])
                return;
            end;
            ax1 = findobj('tag','backeeg','parent',fig);
            tmppos = get(ax1, 'currentpoint');
            ax1 = findobj('tag','eegaxis','parent',fig); % axes handle
            % plot vertical line
            yl = ylim;
            plot([ tmppos tmppos ], yl, 'color', [0.8 0.8 0.8]);
            
            if g.trialstag ~= -1,
                lowlim = round(g.time*g.trialstag+1);
            else, lowlim = round(g.time*g.srate+1);
            end;
            data = get(ax1,'UserData');
            datapos = max(1, round(tmppos(1)+lowlim));
            datapos = min(datapos, g.frames);
            
            figure; topoplot(data(:,datapos), g.eloc_file);
            if g.trialstag == -1,
                latsec = (datapos-1)/g.srate;
                title(sprintf('Latency of %d seconds and %d milliseconds', floor(latsec), round(1000*(latsec-floor(latsec)))));
            else
                trial = ceil((datapos-1)/g.trialstag);
                
                latintrial = eeg_point2lat(datapos, trial, g.srate, g.limits, 0.001);
                title(sprintf('Latency of %d ms in trial %d', round(latintrial), trial));
            end;
            return;
            
            % release button: check window consistency, add to trial boundaries
            % -------------------------------------------------------------------
        case 'defupcom'
            fig = varargin{1};
            g = get(fig,'UserData');
            ax1 = findobj('tag','backeeg','parent',fig);
            g.incallback = 0;
            set(fig,'UserData', g);  % early save in case of bug in the following
            if strcmp(g.mocap,'on'), g.winrej = g.winrej(end,:);end; % nima
            if ~isempty(g.winrej)', ...
                    if g.winrej(end,1) == g.winrej(end,2) % remove unitary windows
                    g.winrej = g.winrej(1:end-1,:);
                    else
                        if g.winrej(end,1) > g.winrej(end,2) % reverse values if necessary
                            g.winrej(end, 1:2) = [g.winrej(end,2) g.winrej(end,1)];
                        end;
                        g.winrej(end,1) = max(1, g.winrej(end,1));
                        g.winrej(end,2) = min(g.frames, g.winrej(end,2));
                        if g.trialstag == -1 % find nearest trials boundaries if necessary
                            I1 = find((g.winrej(end,1) >= g.winrej(1:end-1,1)) & (g.winrej(end,1) <= g.winrej(1:end-1,2)) );
                            if ~isempty(I1)
                                g.winrej(I1,2) = max(g.winrej(I1,2), g.winrej(end,2)); % extend epoch
                                g.winrej = g.winrej(1:end-1,:); % remove if empty match
                            else,
                                I2 = find((g.winrej(end,2) >= g.winrej(1:end-1,1)) & (g.winrej(end,2) <= g.winrej(1:end-1,2)) );
                                if ~isempty(I2)
                                    g.winrej(I2,1) = min(g.winrej(I2,1), g.winrej(end,1)); % extend epoch
                                    g.winrej = g.winrej(1:end-1,:); % remove if empty match
                                else,
                                    I2 = find((g.winrej(end,1) <= g.winrej(1:end-1,1)) & (g.winrej(end,2) >= g.winrej(1:end-1,1)) );
                                    if ~isempty(I2)
                                        g.winrej(I2,:) = []; % remove if empty match
                                    end;
                                end;
                            end;
                        end;
                    end;
            end;
            set(fig,'UserData', g);
            eegplot('drawp', 0);
            if strcmp(g.mocap,'on'), show_mocap_for_eegplot(g.winrej); g.winrej = g.winrej(end,:); end; % nima
            
            % push button: create/remove window
            % ---------------------------------
        case 'defdowncom'
            show_mocap_timer = timerfind('tag','mocapDisplayTimer'); if ~isempty(show_mocap_timer),  end; % nima
            fig = varargin{1};
            g = get(fig,'UserData');
            
            ax1 = findobj('tag','backeeg','parent',fig);
            tmppos = get(ax1, 'currentpoint');
            if strcmp(get(fig, 'SelectionType'),'normal');
                
                fig = varargin{1};
                g = get(fig,'UserData');
                ax1 = findobj('tag','backeeg','parent',fig);
                tmppos = get(ax1, 'currentpoint');
                g = get(fig,'UserData'); % get data of backgroung image {g.trialstag g.winrej incallback}
                if g.incallback ~= 1 % interception of nestest calls
                    if g.trialstag ~= -1,
                        lowlim = round(g.time*g.trialstag+1);
                        highlim = round(g.winlength*g.trialstag);
                    else,
                        lowlim  = round(g.time*g.srate+1);
                        highlim = round(g.winlength*g.srate);
                    end;
                    if (tmppos(1) >= 0) & (tmppos(1) <= highlim),
                        if isempty(g.winrej) Allwin=0;
                        else Allwin = (g.winrej(:,1) < lowlim+tmppos(1)) & (g.winrej(:,2) > lowlim+tmppos(1));
                        end;
                        if any(Allwin) % remove the mark or select electrode if necessary
                            lowlim = find(Allwin==1);
                            if g.setelectrode  % select electrode
                                ax2 = findobj('tag','eegaxis','parent',fig);
                                tmppos = get(ax2, 'currentpoint');
                                tmpelec = g.chans + 1 - round(tmppos(1,2) / g.spacing);
                                tmpelec = min(max(tmpelec, 1), g.chans);
                                g.winrej(lowlim,tmpelec+5) = ~g.winrej(lowlim,tmpelec+5); % set the electrode
                            else  % remove mark
                                g.winrej(lowlim,:) = [];
                            end;
                        else
                            if g.trialstag ~= -1 % find nearest trials boundaries if epoched data
                                alltrialtag = [0:g.trialstag:g.frames];
                                I1 = find(alltrialtag < (tmppos(1)+lowlim) );
                                if ~isempty(I1) & I1(end) ~= length(alltrialtag),
                                    g.winrej = [g.winrej' [alltrialtag(I1(end)) alltrialtag(I1(end)+1) g.wincolor zeros(1,g.chans)]']';
                                end;
                            else,
                                g.incallback = 1;  % set this variable for callback for continuous data
                                if size(g.winrej,2) < 5
                                    g.winrej(:,3:5) = repmat(g.wincolor, [size(g.winrej,1) 1]);
                                end;
                                if size(g.winrej,2) < 5+g.chans
                                    g.winrej(:,6:(5+g.chans)) = zeros(size(g.winrej,1),g.chans);
                                end;
                                g.winrej = [g.winrej' [tmppos(1)+lowlim tmppos(1)+lowlim g.wincolor zeros(1,g.chans)]']';
                            end;
                        end;
                        set(fig,'UserData', g);
                        eegplot('drawp', 0);  % redraw background
                    end;
                end;
            elseif strcmp(get(fig, 'SelectionType'),'normal');
                
                
            end;
        otherwise
            error(['Error - invalid eegplot() parameter: ',data])
    end
end
% Function to show the value and electrode at mouse position
function defmotion(varargin)
fig = varargin{3};
ax1 = varargin{4};
tmppos = get(ax1, 'currentpoint');

if  all([tmppos(1,1) >= 0,tmppos(1,2)>= 0])
    g = get(fig,'UserData');
    if g.trialstag ~= -1,
        lowlim = round(g.time*g.trialstag+1);
    else, lowlim = round(g.time*g.srate+1);
    end;
    if g.incallback
        g.winrej = [g.winrej(1:end-1,:)' [g.winrej(end,1) tmppos(1)+lowlim g.winrej(end,3:end)]']';
        set(fig,'UserData', g);
        eegplot('drawb');
    else
        hh = varargin{6}; % h = findobj('tag','Etime','parent',fig);
        if g.trialstag ~= -1,
            tmpval = mod(tmppos(1)+lowlim-1,g.trialstag)/g.trialstag*(g.limits(2)-g.limits(1)) + g.limits(1);
            if g.isfreq, tmpval = tmpval/1000 + g.freqs(1); end
            set(hh, 'string', num2str(tmpval));
        else
            tmpval = (tmppos(1)+lowlim-1)/g.srate;
            if g.isfreq, tmpval = tmpval+g.freqs(1); end
            set(hh, 'string', num2str(tmpval)); % put g.time in the box
        end;
        ax1 = varargin{5};% ax1 = findobj('tag','eegaxis','parent',fig);
        tmppos = get(ax1, 'currentpoint');
        tmpelec = round(tmppos(1,2) / g.spacing);
        tmpelec = min(max(double(tmpelec), 1),g.chans);
        labls = get(ax1, 'YtickLabel');
        hh = varargin{8}; % hh = findobj('tag','Eelec','parent',fig);  % put electrode in the box
        if ~g.envelope
            set(hh, 'string', labls(tmpelec+1,:));
        else
            set(hh, 'string', ' ');
        end
        hh = varargin{7}; % hh = findobj('tag','Evalue','parent',fig);
        if ~g.envelope
            eegplotdata = get(ax1, 'userdata');
            set(hh, 'string', num2str(eegplotdata(g.chans+1-tmpelec, min(g.frames,max(1,double(round(tmppos(1)+lowlim)))))));  % put value in the box
        else
            set(hh,'string',' ');
        end
    end;
end

% function not supported under Mac
% --------------------------------
function [reshist, allbin] = myhistc(vals, intervals);

reshist = zeros(1, length(intervals));
allbin = zeros(1, length(vals));

for index=1:length(vals)
    minvals = vals(index)-intervals;
    bintmp  = find(minvals >= 0);
    [mintmp indextmp] = min(minvals(bintmp));
    bintmp = bintmp(indextmp);
    
    allbin(index) = bintmp;
    reshist(bintmp) = reshist(bintmp)+1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function infoBV = read_vhdr_BV(file)
fid = fopen(file);
ilabel = 1;
while ~feof(fid)
    fline = fgetl(fid);
    if ~isempty(fline)
        [rem,tok]=strtok(fline,'=');
        if numel(rem)<2
        elseif  strcmp(rem,'NumberOfChannels')
            infoBV.NumberOfChannels = str2num(tok(2:end));
        elseif strcmp(rem,'DataType')
            infoBV.DataType = tok(2:end);
        elseif strcmp(rem,'DataFormat')
            infoBV.dataformat = (tok(2:end));
        elseif strcmp(rem,'DataOrientation')
            infoBV.DataOrientation = (tok(2:end));
        elseif strcmp(rem,'BinaryFormat')
            infoBV.BinaryFormat = tok(2:end);
        elseif strcmp(rem,'DataPoints')
            infoBV.DataPoints = str2num(tok(2:end));
        elseif strcmp(rem,'SamplingInterval')
            infoBV.SamplingInterval=str2num(tok(2:end));
        elseif strcmp(rem,'Layers')
            infoBV.Layers=str2num(tok(2:end));
        elseif strcmp(rem,'SegmentDataPoints')
            infoBV.SegmentDataPoints = str2num(tok(2:end));
        elseif strcmp(rem,'SegmentDataPointsPre')
            infoBV.SegmentDataPointsPre = str2num(tok(2:end));
        elseif strcmp(rem,'SegmentDataPointsPost')
            infoBV.SegmentDataPointsPost = str2num(tok(2:end));
        elseif strcmp(rem(1:2),'Ch') %& strcmp(tok(1:5),'=SEEG'
            infoBV.label{ilabel} = tok; %(7:12);
            ilabel = ilabel+1;
        end
    else
    end
end
fclose(fid);



% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = 1;

nch = size(PMI{currentsub}.data(cf).HRF.AvgC,2);
MeasListAct=PMI{currentsub}.data(cf).MeasListAct;
MeasListAct(:)=1;
plotLst= find( PMI{currentsub}.data(cf).MeasListAct);
if ~isempty(plotLst)
    PMI{currentsub}.plotLst = plotLst;
    set(handles.edit_plotLst,'string',mat2str(plotLst))
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
    
end
handles.newlist=1;
set(guiHOMER,'UserData',PMI)
guidata(hObject, handles);
updatedisplay(handles)


% --- Executes on button press in radio_LPF.
function radio_LPF_Callback(hObject, eventdata, handles)
% hObject    handle to radio_LPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_LPF


% --- Executes on button press in btn_substractICA.
function btn_substractICA_Callback(hObject, eventdata, handles)
% hObject    handle to btn_substractICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_addICA.
function btn_addICA_Callback(hObject, eventdata, handles)
% hObject    handle to btn_addICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
% currentsub=1;
% PMI = get(guiHOMER,'UserData');
% cf = PMI{currentsub}.currentFile;
% d = PMI{currentsub}.data(cf).HRF.AvgC;
% indt = [PMI{currentsub}.tmpICA.indt(1):PMI{currentsub}.tmpICA.indt(2)];%Time indice
% intensnorm = d(indt,:);
% %Detrent DATA segment for centrering
% X = 1:1:size(intensnorm,1);
% Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
% Mb2 =  intensnorm(1,:)'; %offset
% A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
% spar = intensnorm - A;
%
%     A = PMI{currentsub}.tmpICA.Factors{1};
%     B = PMI{currentsub}.tmpICA.Factors{2};
%     C = PMI{currentsub}.tmpICA.Factors{3};
%     listgood = PMI{currentsub}.tmpICA.listgood;
%     ComponentToKeep = PMI{1}.tmpICA.selected;
%     Ac = A(ComponentToKeep,:); Bc = B(ComponentToKeep,:); Cc = C(:,ComponentToKeep);
%    % figure;plot(reshape(Xm,[numel(indt),size(d,2)]))
%     Xm = (Ac')*Cc';
%
%     %data = d(indt,:)%,d(indt,end/2+1:end));
%     d(indt,listgood) = d(indt,listgood)+Xm;
%     %d(indt,listgood) = d(indt,listgood)-reshape(Xm,[numel(indt),size(d,2)]);
%
%
%     PMI{currentsub}.data(cf).HRF.AvgC= d;
%     set(guiHOMER,'UserData',PMI);
%     updatedisplay(handles);



% --- Executes on button press in btn_wavelet.
function btn_wavelet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_wavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.radio_WAVELET,'value',1)

time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
if isempty(time_start)|isempty(time_stop)
    set(handles.text_Advice,'string','Please enter time start and time stop to define time window')
    return
end
module = get(handles.popupmenu_file,'value')
filelist = handles.NIRS.Dt.EEG.pp.p{module}
timelist = 0
EEG = read_eeg_multiplex_besa(filelist,timelist,time_start, time_stop )


gui_MorletTG(time_start,time_stop,EEG);


% --- Executes on button press in btn_Video.
function btn_Video_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
cf = 1
PMI = get(guiHOMER,'UserData');
   
topbar = max(max(PMI{1}.data.HRF.AvgC(:,PMI{1}.plotLst)));
val = get(handles.btn_Video,'string')
if str2num(val(end)) %,'on') %si déjà en court le fermer ne pas redémaré attendre
    set(handles.btn_Video,'string','Video 0')
        return 
%     close(h)
%     guidata(hObject,handles);
%     
else
    time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
    time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
    id = get(handles.popupmenu_file,'value');
    filenamevideo = handles.NIRS.Dt.Video.pp(end).p{id};
   
    if isfield(handles.NIRS.Dt.Video.pp(end),'sync_timesec')
         offset =  handles.NIRS.Dt.Video.pp(end).sync_timesec{id};
         if isnan(offset)
             msgbox('No synchronised by trig')
             return
         end
    end
    
     try
        [heure,m,s] = hms(seconds(offset+time_start));
         set(handles.text_Advice,'string',['Time ',sprintf('%02.0f',heure),':',sprintf('%02.0f',m),':',sprintf('%02.0f',s),' ', filenamevideo]);
    catch
         set(handles.text_Advice,'string',['Time ', num2str(round(offset+time_start)),'s']);
     end
     if ~isfield(PMI{1},'videoview') 
         PMI{1}.videoview = 1
     end
    if PMI{1}.videoview
        h= figure
        text(0.2, 0.5, 'Please wait video is loading','fontsize',14)
    else
        return
    end
 
    
     if isempty(offset)
        msgbox('Segment first to synchronised video with NIRS')
        return
     end
     if isempty(time_start)|isempty(time_stop)
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        msgbox({'\fontsize{14} Enter time start and time stop'}, '','help', CreateStruct)
        return
    end
    %set(handles.uipanel9,'visible','off');
    %set(handles.uipanel_manuelNorm,'visible','off');
    if (  time_stop- time_start  )>30
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        msgbox({'\fontsize{14} Please select an interval time' ; 'start and stop shorter than 30 seconds'}, '','help', CreateStruct)
        return
    end
    
    if isfield(PMI{1},'videooption')
        if ~isfield( PMI{1}.videooption,'brighness')
            PMI{1}.videooption.brighness = 0
        end
        if ~isfield( PMI{1}.videooption,'codec')
             PMI{1}.videooption.codec = 'VideoReader'; 
        end
    else
        PMI{1}.videooption.brighness = 0
        PMI{1}.videooption.codec = 'VideoReader'; 
    end

%     %Test to use mmread codec or defaut computer install codec 
%     try
% %         [video_example, audio_example] = mmread(filenamevideo,[],[1,1.04]);
% %         if ~isempty(video_example)
% %             readoption=1 %use mmread codec
% %         else
%             [video_example,audio_example]=ReadMyAudioVideo(filenamevideo,[1,1.04]);
%             readoption=2
%         end  
%     catch
%          readoption=1
%     end
              

    %load the frame
    if strcmp( PMI{1}.videooption.codec,'mmread')
        try
        [video_example,audio_example]=mmread(filenamevideo,[],[offset+time_start,offset+time_stop]);
        catch
            msgbox(['Unable to read video ', filenamevideo])
            return
        end
    elseif strcmp( PMI{1}.videooption.codec,'VideoReader') 
        [video_example,audio_example]=ReadMyAudioVideo(filenamevideo,[offset+time_start,offset+time_stop]);
    end
    
    set(handles.btn_Video,'string','Video 1')
     nbframe = numel(video_example.frames)
    
    if ~isempty(audio_example)%& iframe==1
         sound(audio_example.data(:,1),audio_example.rate);
    end    
     set(gca,'xticklabel','')
%     set(gca,'yticklabel','')
    for iframe=1:nbframe        
     tic
     ini = video_example.frames(iframe).cdata+PMI{1}.videooption.brighness;
     image(ini);%PMI{1}.videooption.brighness
     drawnow      
     a=toc;
      pause(1/video_example.rate-a)
    end
    close(h)
    set(handles.btn_Video,'string','Video 0')    
end    
    
    
    
    %play the audio video stream
%       set(handles.axes_video,'visible','on');
% 	  haxesvideo = handles.axes_video;
%       set(gca,'xticklabel','')
%       set(gca,'yticklabel','')
%       haxesbar = handles.display_axes;
%       axes(haxesbar);  
%       hold on      
%       hline = plot([time_start time_start],[0,topbar],'k')
%       pause(0.001)
%     nbframe = numel(video_example.frames)
%     if ~isempty(audio_example)
%     nbbyframe = round(size(audio_example.data,1)/nbframe)
%     end
   % sound(audio_example.data,audio_example.rate)

   %  p = audioplayer(audio_example.data,
   %  audio_example.rate);
    %       play(p, [1 (get(p, 'SampleRate') * 3)]);
%   nbframe = nbframe - 1  
%    
%         if ~isempty(audio_example)%& iframe==1
%        % sound(audio_example.data((1+nbbyframe*(iframe-1):nbbyframe*iframe),1),audio_example.rate);
%          sound(audio_example.data(:,1),audio_example.rate);
%         end
%     h= figure
%     for iframe=1:nbframe        
%      %axes(haxesvideo);
%      tic
%      image(video_example.frames(iframe).cdata);
%      drawnow
%      a=toc;
%       pause(1/video_example.rate-a)
%     
% %     
% %          if ~mod(iframe,20)   
% %              axes(haxesbar);
% %              t = (iframe*1/video_example.rate)+time_start;
% %              set( hline, 'xdata',[t t])
% %              drawnow
% %          end
% %     %    pause(1/video_example.rate)
% %        %   h = plot([(iframe*1/video_example.rate)+time_start (iframe*1/video_example.rate)+time_start],[0,1],'k') 
%     end
%     close(h)
 %   delete(hline);
    
%      [video_example,audio_example]=ReadMyAudioVideo(filenamevideo,[offset+time_start,offset+time_stop]);
%     
%     
%     
%     v=VideoReader(filenamevideo,'CurrentTime',offset+time_start);
%     nbframe = round((time_stop - time_start)*v.FrameRate);
%     set(handles.axes_video,'visible','on');
%     axes(handles.axes_video);
%     filenameaudio = [filenamevideo(1:end-4),'.wav'];
%     try 
%         info = audioinfo( filenameaudio);
%         samples = [round((offset+time_start)*info.SampleRate):round((offset+time_stop)*info.SampleRate)];
%         [y,Fs] = audioread(filenameaudio,[samples(1),samples(end)]); 
%         player = audioplayer(y,Fs)
%         flagaudio =1;
%     catch 
%         flagaudio =0;
%     end
%     %currAxes = axes;
%     id = 1;
% 
%     
%     set(handles.axes_video,'visible','off')
%     set(handles.axes_video,'xticklabel','')
%     set(handles.axes_video,'yticklabel','')
%      if flagaudio
%             play(player)
%      end
%     for i = 1:nbframe
%        
%         val = get(handles.btn_Video,'string');
%         if ~str2num(val(end))
%             if flagaudio
%             stop(player)
%             end
%             break           
%         end
%         if i==1
%             h = plot( [time_start, time_start],[0,1],'k');
%         end
%        
%         vidFrame = readFrame(v); 
%         axes( handles.axes_video)
%         image(vidFrame, 'Parent', handles.axes_video);   
%           set(gca,'xticklabel','')
%           set(gca,'yticklabel','')
%         idelta = i * 1/v.FrameRate;
%        % if id ==15
%             axes(handles.display_axes);            
%             hold on;
%             delete(h);
%             h = plot( [time_start+idelta, time_start+idelta],[0,1],'k');
%             id = 1;
%       %  end
%       %   pause(1/v.FrameRate);
%         
%         id = id+1;
%     end
   % axes(handles.axes_video)
   % cla
   % set(handles.axes_video,'visible','off')
   % set(handles.uipanel9,'visible','on')
   % set(handles.uipanel_manuelNorm,'visible','on')

% --- Executes on button press in btn_SNR.
function btn_SNR_Callback(hObject, eventdata, handles)
% hObject    handle to btn_SNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
cf = 1

PMI = get(guiHOMER,'UserData');
time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
tstart=find(PMI{currentsub}.data(cf).HRF.tHRF<time_start);
tstop=find(PMI{currentsub}.data(cf).HRF.tHRF<time_stop);
indt = tstart(end):tstop(end);
plotLst =  PMI{currentsub}.plotLst;
ML = PMI{currentsub}.data(cf).MeasList;
for ich=1:numel(plotLst)
    idnoise = find(PMI{currentsub}.data.HRF.noise(indt,plotLst(ich))==1);
    idsignal = find(PMI{currentsub}.data.HRF.noise(indt,plotLst(ich))==0);
    
    SNR(ich) = log10(nanvar(PMI{currentsub}.data.HRF.AvgC(indt(idsignal),plotLst(ich)))./nanvar(PMI{currentsub}.data.HRF.AvgC(indt(idnoise),plotLst(ich))));
    strDet = SDDet2strboxy_ISS(ML(plotLst(ich),2));
    strSrs = SDPairs2strboxy_ISS(ML(plotLst(ich),1));
    label{ich} = [strDet, ' ',strSrs ];
end
NIRSMAT  = get(handles.edit_nirsmat, 'String');
val = get(handles.popupmenu_module, 'String');
id = get(handles.popupmenu_module, 'value');
module = val{id};
id = get(handles.popupmenu_file,'value')
val = get(handles.popupmenu_file, 'String');
idfile = val{id};
id = get(handles.popupmenu_file,'value')
val = get(handles.popupmenu_file, 'String');
idfile = val{id};
[file,path]=uiputfile(['SNR',sprintf('%03.1f',time_start),sprintf('_%03.1f',time_stop),module,idfile,'.xls']);
matxls = [label;num2cell(SNR)]
xlswrite([path,file],matxls);
SNR = []
for ich=size(PMI{currentsub}.data.HRF.AvgC,2)
    idnoise = find(PMI{currentsub}.data.HRF.noise(indt,ich)==1);
    idsignal = find(PMI{currentsub}.data.HRF.noise(indt,ich)==0);
    SNR(ich) = log10(nanvar(PMI{currentsub}.data.HRF.AvgC(indt(idsignal),ich))./nanvar(PMI{currentsub}.data.HRF.AvgC(indt(idnoise),ich)));
    strDet = SDDet2strboxy_ISS(ML(ich,2));
    strSrs = SDPairs2strboxy_ISS(ML(ich,1));
    label{ich} = [strDet, ' ',strSrs ];
end

save([path,file(1:end-3),'mat'],'SNR','label','NIRSMAT','module','idfile','time_start','time_stop','-mat' );




% --- Executes on selection change in listbox_CorrectionDecomposition.
function listbox_CorrectionDecomposition_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_CorrectionDecomposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_CorrectionDecomposition contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_CorrectionDecomposition

try
    [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
    load(fullfile(pathstr,'CorrectionApply.mat'))
catch
    set(handles.listbox_CorrectionDecomposition,'string','')
    return
end

global currentsub
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf  = 1;

%set(handles.listbox_CorrectionApply,'string','new')
id = get(handles.listbox_CorrectionDecomposition,'value');

% PARCORR(1,id).file;
% PARCORR(1,id).module;
actualmodule = get(handles.popupmenu_module,'value');
actualfile = get(handles.popupmenu_file,'value');
set(handles.popupmenu_module, 'value',PARCORR(1,id).module);
set(handles.popupmenu_file, 'value',PARCORR(1,id).file);
PMI{currentsub}.plotLst = PARCORR(1,id).listgood;
plotLst = PMI{currentsub}.plotLst;
PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];

set(handles.edit_time_start,'string', num2str(PMI{currentsub}.data(cf).HRF.tHRF(PARCORR(1,id).indt(1))));
set(handles.edit_time_stop,'string', num2str(PMI{currentsub}.data(cf).HRF.tHRF(PARCORR(1,id).indt(end))));
if isfield(PARCORR(1,id),'topo')
    PMI{currentsub}.tmptopo =  PARCORR(1,id).topo;
    PMI{currentsub}.tmplistgood =  PARCORR(1,id).listgood;
else
    PMI{currentsub}.tmptopo = [0 ,0];
    PMI{currentsub}.tmplistgood = [1,2];
end
set(guiHOMER,'UserData',PMI);
guidata(hObject, handles);

if PARCORR(1,id).module == actualmodule & actualfile ==PARCORR(1,id).file
else
        updatedata(handles,1,1);
end
step = str2num(get(handles.edit_stepvalue,'string'));

if strcmp(get(handles.context_enable_autozoom,'Checked'),'on')	
tstart = str2num(get(handles.edit_time_start,'string'));
tstop = str2num(get(handles.edit_time_stop,'string'));
set(handles.edit_start,'string', num2str(tstart-3));
set(handles.edit_duration,'string', num2str(tstop-tstart+6));
end
updatedisplay(handles);
guidata(handles.figure1, handles);
if strcmp(get(handles.context_listbox_Correction_newfigure,'checked'),'on')
    figure;hold on
    set(handles.context_listbox_Correction_newfigure,'checked','off')
    guidata(hObject, handles);
else
    axes(handles.display_axes);
end

for i=1:numel(PARCORR(1,id).listgood)
 idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, PARCORR(1,id).listgood(i),...
            numel(PMI{currentsub}.color)/3);
h=plot(PMI{currentsub}.data(cf).HRF.tHRF(PARCORR(1,id).indt(1):PARCORR(1,id).indt(end)),PARCORR(1,id).Xm(:,i,1)+step,'color', PMI{currentsub}.color(idxc,:));
       srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(PARCORR(1,id).listgood(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(PARCORR(1,id).listgood(i),2));
set(h,'displayname',['ch',num2str(PARCORR(1,id).listgood(i)),'_',srs, '_', det]);
set(h,'linewidth',2);
%set(gca,'fontsize',12)
xlabel('Time (s)');


end


%set(handles.

% --- Executes during object creation, after setting all properties.
function listbox_CorrectionDecomposition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_CorrectionDecomposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_restoreInterval.
function btn_restoreInterval_Callback(hObject, eventdata, handles)
% hObject    handle to btn_restoreInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button = questdlg(['Are you sure you want to restore the interval all correction will be deleted'],'Restore','Yes','No','Cancel','No');
switch button
    case 'No'
        return
    case 'Cancel'
        return
end
global currentsub
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf  = 1;

time_start = str2num(get(handles.edit_time_start,'string'));
time_stop = str2num(get(handles.edit_time_stop,'string'));
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
end
if isempty(time_stop); indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);else;
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
end
if isempty(indstart);
    indstart=1;
else
    indstart = indstart(end);
end

indstop = indstop(end);

idmodule = get(handles.popupmenu_module, 'value');
idfile = get(handles.popupmenu_file,'value');
NC = handles.NIRS.Cf.H.C.N;
xname =handles.NIRS.Dt.fir.pp(end-1).p{idfile};
%[pathstr, name, ext] = fileparts(handles.NIRS.Dt.fir.pp(idmodule).p{idfile})
dprevious= fopen_NIR(xname,NC)';
PMI{currentsub}.data(cf).HRF.AvgC(indstart:indstop,:) = dprevious(indstart:indstop,:);
outfile = handles.NIRS.Dt.fir.pp(end).p{idfile};
fwrite_NIR(outfile,PMI{currentsub}.data(cf).HRF.AvgC');
set(guiHOMER,'UserData',PMI)

%Look at the CorrectionPARAFAC to remove from the list
%Correction item will be remove if they overlap in time with the restore process.
[pathstr, name, ext] = fileparts(handles.NIRSpath{1});
try
    load(fullfile(pathstr,'CorrectionApply.mat'));
catch
    return
end
%find the time overlap
idremove = [];

binidremove= zeros(numel(PMI{currentsub}.data(cf).HRF.tHRF),1);
binidremove(indstart:indstop) = 1;
for icorr=1:numel(PARCORR)
    if PARCORR(icorr).file == idfile
        binindt = zeros(numel(PMI{currentsub}.data(cf).HRF.tHRF),1);
        binindt(PARCORR(icorr).indt)=1;
        if sum(binindt.*binidremove)
            if strcmp('Offset Ajustement',PARCORR(icorr).type)
                idkeep = icorr;
            else
                idremove = [idremove,icorr];
            end
        end
    end
end
PARCORR.type
if ~isempty(idremove)
    PARCORR(idremove) =[];
    if isempty(PARCORR)
        delete(fullfile(pathstr,'CorrectionApply.mat'));
        set(handles.listbox_CorrectionDecomposition,'string','');
        set(handles.listbox_CorrectionDecomposition,'value',1);
    else
        save(fullfile(pathstr,'CorrectionApply.mat'),'PARCORR');
        for i=1:numel(PARCORR)
            CORRlist{i} = PARCORR(i).label;
        end
        set(handles.listbox_CorrectionDecomposition,'string',CORRlist);
        set(handles.listbox_CorrectionDecomposition,'value',1);
    end
end
updatedisplay(handles);


% --- Executes on selection change in popupmenu_ROIEEG.
function popupmenu_ROIEEG_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ROIEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_ROIEEG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ROIEEG


% --- Executes during object creation, after setting all properties.
function popupmenu_ROIEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ROIEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_Component.
function listbox_Component_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_Component contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_Component
[pathstr, name, ext] = fileparts(handles.NIRSpath{1});
try
    load(fullfile(pathstr,'SelectedFactors.mat'))
catch
    return
end
global currentsub
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf  = 1;

%set(handles.listbox_CorrectionApply,'string','new')
id = get(handles.listbox_Component,'value');
if numel(PARCOMP)<id
    id = numel(PARCOMP);
    set(handles.listbox_Component,'value',id)
end
% PARCOMP(1,id).file
%PARCOMP(1,id).module
actualmodule = get(handles.popupmenu_module,'value');
actualfile = get(handles.popupmenu_file,'value');
set(handles.popupmenu_module, 'value',PARCOMP(1,id).module);
set(handles.popupmenu_file, 'value',PARCOMP(1,id).file);
PMI{currentsub}.plotLst = PARCOMP(1,id).listgood;
plotLst = PMI{currentsub}.plotLst;
PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];

set(handles.edit_time_start,'string', num2str(PMI{currentsub}.data(cf).HRF.tHRF(PARCOMP(1,id).indt(1))));
set(handles.edit_time_stop,'string', num2str(PMI{currentsub}.data(cf).HRF.tHRF(PARCOMP(1,id).indt(end))));
if isfield(PARCOMP(1,id),'topo')
    PMI{currentsub}.tmptopo =  PARCOMP(1,id).topo;
    PMI{currentsub}.tmplistgood =  PARCOMP(1,id).listgood;
else
    PMI{currentsub}.tmptopo =  zeros(numel(PARCOMP(1,id).listgood),1);
    PMI{currentsub}.tmplistgood =  PARCOMP(1,id).listgood;
end
set(guiHOMER,'UserData',PMI)
guidata(hObject, handles);
if PARCOMP(1,id).module == actualmodule & actualfile ==PARCOMP(1,id).file
else
        updatedata(handles,1,1);
        disp('update')
end
step = str2num(get(handles.edit_stepvalue,'string'));
if strcmp(get(handles.context_enable_autozoom,'Checked'),'on')	
tstart = str2num(get(handles.edit_time_start,'string'));
tstop = str2num(get(handles.edit_time_stop,'string'));
set(handles.edit_start,'string', num2str(tstart-3));
set(handles.edit_duration,'string', num2str(tstop-tstart+6));
end
    
updatedisplay(handles)


        
if strcmp(get(handles.context_listbox_Component_newfigure,'checked'),'on')
    figure;hold on
    set(handles.context_listbox_Component_newfigure,'checked','off')
    guidata(hObject, handles);
else
    axes(handles.display_axes)
end
try
for i=1:numel(PARCOMP(1,id).listgood)
 idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, PARCOMP(1,id).listgood(i),...
            numel(PMI{currentsub}.color)/3);
         h=plot(PMI{currentsub}.data(cf).HRF.tHRF(PARCOMP(1,id).indt(1):PARCOMP(1,id).indt(end)),PARCOMP(1,id).Xm(:,i,1)+step,'color',PMI{currentsub}.color(idxc,:));
        srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(PARCOMP(1,id).listgood(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(PARCOMP(1,id).listgood(i),2));
        set(h,'displayname',['ch',num2str(PARCOMP(1,id).listgood(i)),'_',srs, '_', det]);
end
catch
end
guidata(handles.figure1, handles);
% % %PLOTALEJANDRA
% figure
% lstSV = 1
%  npos = numel(PARCOMP(1,id).v(:,lstSV))/2;
% plot(mean(PARCOMP(1,id).Xm(:,:,1),1))
%  plot(PMI{currentsub}.data(cf).HRF.tHRF(PARCOMP(1,id).indt),PARCOMP(1,id).Xm(:,40,1))
% %
%         tempPCA=reshape(PARCOMP(1,id).v(:,1),npos,2);
%         figure
% plot(tempPCA')

% --- Executes during object creation, after setting all properties.
function listbox_Component_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_Component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_GetDecomposition.
function btn_GetDecomposition_Callback(hObject, eventdata, handles)
% hObject    handle to btn_GetDecomposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Left arrow the component will be store in the list of component,
%PARCORR.Xm
%PARCORR.indt same for all structure
%to all structure

%same as + but just to keep in the memory for futher analysis
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
d = PMI{currentsub}.data(cf).HRF.AvgC;
[pathstr, name, ext] = fileparts(handles.NIRSpath{1});
idcurrent = get(handles.listbox_Component,'value');
idval = get(handles.popupmethodselected,'value');
listmethod = get(handles.popupmethodselected,'string');
if strcmp(listmethod{idval},'Parafac') %get(handles.popupmethodselected,'value')==1  %extract PARAFAC
    if isfield (PMI{currentsub}, 'tmpPARAFAC')
        indt = [PMI{currentsub}.tmpPARAFAC.indt(1):PMI{currentsub}.tmpPARAFAC.indt(2)];%Time indice
        intensnorm = d(indt,:);
        %Detrent DATA segment for centrering
        PMI{currentsub}.tmpPARAFAC
        X = 1:1:size(intensnorm,1);
        Mb1 = ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
        Mb2 = intensnorm(1,:)'; %offset
        A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
        spar = intensnorm - A;
        A = PMI{currentsub}.tmpPARAFAC.Factors{1};
        B = PMI{currentsub}.tmpPARAFAC.Factors{2};
        C = PMI{currentsub}.tmpPARAFAC.Factors{3};
        listgood = PMI{currentsub}.tmpPARAFAC.listgood;
        ComponentToKeep = PMI{1}.tmpPARAFAC.selected;
        Ac = A(:,ComponentToKeep); Bc = B(:,ComponentToKeep); Cc = C(:,ComponentToKeep);
%         figure;plot(reshape(Xm,[numel(indt),size(d,2)]))
        [Xm]=nmodel(({Ac,Bc,Cc}));
        data = cat(3,d(indt,1:end/2),d(indt,end/2+1:end));
        data(:,listgood,:) = data(:,listgood,:)-Xm;
%         if sum(data(:))==PMI{currentsub}.tmpPARAFAC.checksumd    %check sum protection
%             data(:,listgood,:) = data(:,listgood,:)-Xm;
%         else
%             msgbox('Please update the decomposition, this is not the one for the current data')
%             return
%         end
    %
    %save in a structure all apply correction
        try
            load(fullfile(pathstr,'SelectedFactors.mat'))
            newfile = 0;
        catch
            %donot exist create the stucture
            PARCOMP.file= get(handles.popupmenu_file,'value');
            fileall = get(handles.popupmenu_file,'string');
            PARCOMP.filestr =  fileall{get(handles.popupmenu_file,'value')};
            PARCOMP.module  =get(handles.popupmenu_module, 'value');
            moduleall = get(handles.popupmenu_module,'string');
            PARCOMP.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
            PARCOMP.listgood =  listgood
            PARCOMP.indt = indt %indice de temps.
            PARCOMP.data = data(:,listgood,:);
            PARCOMP.Xm = Xm;
            PARCOMP.FacA = PMI{currentsub}.tmpPARAFAC.Factors{1};
            PARCOMP.FacB = PMI{currentsub}.tmpPARAFAC.Factors{2};
            PARCOMP.FacC = PMI{currentsub}.tmpPARAFAC.Factors{3};
            PARCOMP.ComponentToKeep = ComponentToKeep;
            labelid  = get(handles.edit_selectedlabel,'string');
            PARCOMP.label= [labelid,'PARAFAC' sprintf('%03.0f',size(PARCOMP,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
            PARCOMP.type = 'PARAFAC';
            FacSpatial = PMI{currentsub}.tmpPARAFAC.Factors{2}
            selected = PMI{currentsub}.tmpPARAFAC.selected;
            PARCOMP.topo = FacSpatial(:,selected);
                
            newfile = 1;
        end
        if newfile == 0
            id = numel(PARCOMP);
            PARCOMP(id+1).file= get(handles.popupmenu_file,'value');
            fileall = get(handles.popupmenu_file,'string');
            PARCOMP(id+1).filestr =  fileall{get(handles.popupmenu_file,'value')};
            PARCOMP(id+1).module  =get(handles.popupmenu_module, 'value');
            moduleall = get(handles.popupmenu_module,'string');
            PARCOMP(id+1).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
            PARCOMP(id+1).listgood =  listgood;
            PARCOMP(id+1).indt = indt; %indice de temps.
            PARCOMP(id+1).data = data(:,listgood,:);
            PARCOMP(id+1).Xm = Xm;
            PARCOMP(id+1).FacA = PMI{currentsub}.tmpPARAFAC.Factors{1};
            PARCOMP(id+1).FacB = PMI{currentsub}.tmpPARAFAC.Factors{2};
            PARCOMP(id+1).FacC = PMI{currentsub}.tmpPARAFAC.Factors{3};
            PARCOMP(id+1).ComponentToKeep = ComponentToKeep;
            labelid  = get(handles.edit_selectedlabel,'string');
            PARCOMP(id+1).label= [labelid,'PARAFAC' sprintf('%03.0f',size(PARCOMP,2)),' ',fileall{get(handles.popupmenu_file,'value')}]
            PARCOMP(id+1).type = 'PARAFAC';
            FacSpatial = PMI{currentsub}.tmpPARAFAC.Factors{2}
            selected = PMI{currentsub}.tmpPARAFAC.selected;
            PARCOMP(id+1).topo = FacSpatial(:,selected);;

        
            %ajoute a la fin et replace au centre
            PARCOMP = [PARCOMP(1,1:idcurrent-1),PARCOMP(1,id+1), PARCOMP(1,idcurrent:id)]
        end
        
        for i=1:(numel(PARCOMP))
            CORRlist{i} = PARCOMP(i).label
        end
        set(handles.listbox_Component,'string',CORRlist)
%         PMI{currentsub}.data(cf).HRF.AvgC(indt,:) = reshape(data,[numel(indt),size(d,2)]);
    else
        msgbox('Press the run button before the component could be add to the list');
        return
    end
elseif strcmp(listmethod{idval},'PCA')  %extract PCA get(handles.radio_PCA,'value') %substract
    if isfield (PMI{currentsub}, 'tmpPCA')
        %Detrent DATA segment for centrering
        indt = [PMI{currentsub}.tmpPCA.indt(1):PMI{currentsub}.tmpPCA.indt(2)];%Time indice
        intensnorm = d(indt,:);
        X = 1:1:size(intensnorm,1);
        Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
        Mb2 =  intensnorm(1,:)'; %offset
        A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
        spar = intensnorm - A;
        listgood = PMI{currentsub}.tmpPCA.listgood;
        lstSV = PMI{1}.tmpPCA.selected;
        u =PMI{currentsub}.tmpPCA.u ;
        s =PMI{currentsub}.tmpPCA.s ;
        v = PMI{currentsub}.tmpPCA.v ;
        temp = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
        data = d(indt,:);
        data(:,listgood) = data(:,listgood)- temp;
%         data = d(indt,:);
%         if sum( data(:) )==PMI{currentsub}.tmpPCA.checksumd    %check sum protection
%             data(:,listgood) = data(:,listgood)- temp;
%         else
%             msgbox('Please update the decomposition, this is not the one for the current data')
%             return
%         end
%         
        %save in a structure all apply correction
        [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
        try
            load(fullfile(pathstr,'SelectedFactors.mat'))
            newfile = 0
        catch
        %donot exist create the stucture
            PARCOMP.file= get(handles.popupmenu_file,'value');
            fileall = get(handles.popupmenu_file,'string');
            PARCOMP.filestr =  fileall{get(handles.popupmenu_file,'value')};
            PARCOMP.module  =get(handles.popupmenu_module, 'value');
            moduleall = get(handles.popupmenu_module,'string');
            PARCOMP.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
            PARCOMP.listgood =  listgood
            PARCOMP.indt = indt %indice de temps.
            PARCOMP.data = data(:,listgood,:);
            %PARCOMP.Xm = temp;
            PARCOMP.u =  PMI{currentsub}.tmpPCA.u;
            PARCOMP.s = PMI{currentsub}.tmpPCA.s;
            PARCOMP.v = PMI{currentsub}.tmpPCA.v;
            lstSV = PMI{1}.tmpPCA.selected;
            PARCOMP.Xm = PARCOMP.u(:,lstSV)*PARCOMP.s(lstSV,lstSV)*PARCOMP.v(:,lstSV)';
            PARCOMP.ComponentToKeep = PMI{1}.tmpPCA.selected;
            labelid  = get(handles.edit_selectedlabel,'string');
            PARCOMP.label= [labelid,'PCA' sprintf('%03.0f',size(PARCOMP,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
            PARCOMP.type = 'PCA';
            PARCOMP.topo = s(lstSV,lstSV)*v(:,lstSV)';
            newfile = 1;
        end
        if newfile == 0
            id = numel(PARCOMP);
            PARCOMP(id+1).file= get(handles.popupmenu_file,'value');
            fileall = get(handles.popupmenu_file,'string');
            PARCOMP(id+1).filestr =  fileall{get(handles.popupmenu_file,'value')};
            PARCOMP(id+1).module  =get(handles.popupmenu_module, 'value');
            moduleall = get(handles.popupmenu_module,'string');
            PARCOMP(id+1).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
            PARCOMP(id+1).listgood =  listgood;
            PARCOMP(id+1).indt = indt; %indice de temps.
            PARCOMP(id+1).data = data(:,listgood,:);
            PARCOMP(id+1).u =  PMI{currentsub}.tmpPCA.u;
            PARCOMP(id+1).s = PMI{currentsub}.tmpPCA.s;
            PARCOMP(id+1).v = PMI{currentsub}.tmpPCA.v;
            PARCOMP(id+1).ComponentToKeep = PMI{1}.tmpPCA.selected;
            lstSV = PMI{1}.tmpPCA.selected;
            PARCOMP(id+1).Xm = PARCOMP(id+1).u(:,lstSV)*PARCOMP(id+1).s(lstSV,lstSV)*PARCOMP(id+1).v(:,lstSV)';
            labelid  = get(handles.edit_selectedlabel,'string');
            PARCOMP(id+1).label= [labelid,'PCA' sprintf('%03.0f',size(PARCOMP,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
            PARCOMP(id+1).type = 'PCA';
            PARCOMP(id+1).topo = s(lstSV,lstSV)*v(:,lstSV)';
        end
    
        for i=1:numel(PARCOMP)
            CORRlist{i} = PARCOMP(i).label;
        end
        set(handles.listbox_Component,'string',CORRlist);
        %PMI{currentsub}.data(cf).HRF.AvgC(indt,:) = data;
    else
        msgbox('Press the run button before the component could be add to the list');
        return
    end
    
elseif strcmp(listmethod{idval},'ICA')   %extract ICA Component
    if isfield (PMI{currentsub}, 'tmpICA')
        indt = [PMI{currentsub}.tmpICA.indt(1):PMI{currentsub}.tmpICA.indt(2)];%Time indice
        listgood = PMI{currentsub}.tmpICA.listgood;
        selected = PMI{currentsub}.tmpICA.selected;
        A =  PMI{currentsub}.tmpICA.Factors{1};
        C = PMI{currentsub}.tmpICA.Factors{3};
        Xm = (C(:,PMI{currentsub}.tmpICA.selected)*A(PMI{currentsub}.tmpICA.selected,:))';
        try
            load(fullfile(pathstr,'SelectedFactors.mat'));
            newfile = 0;
        catch
            %donot exist create the stucture
            PARCOMP.file= get(handles.popupmenu_file,'value');
            fileall = get(handles.popupmenu_file,'string');
            PARCOMP.filestr =  fileall{get(handles.popupmenu_file,'value')};
            PARCOMP.module  =get(handles.popupmenu_module, 'value');
            moduleall = get(handles.popupmenu_module,'string');
            PARCOMP.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
            PARCOMP.listgood =  listgood;
            PARCOMP.Factors =  PMI{currentsub}.tmpICA.Factors;
            PARCOMP.indt = indt; %indice de temps.
            PARCOMP.data = PMI{currentsub}.tmpICA.spar ;    
            PARCOMP.ComponentToKeep = PMI{currentsub}.tmpICA.selected;     
            labelid  = get(handles.edit_selectedlabel,'string');
            PARCOMP.Xm = Xm;
            PARCOMP.label= [labelid,'ICA', sprintf('%03.0f',size(PARCOMP,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
            PARCOMP.type = 'ICA';
            PARCOMP.topo =  C(:,PMI{currentsub}.tmpICA.selected);
          
            newfile = 1;
        end
        if newfile == 0
            id = numel(PARCOMP);
            %donot exist create the stucture
            PARCOMP(id+1).file= get(handles.popupmenu_file,'value');
            fileall = get(handles.popupmenu_file,'string');
            PARCOMP(id+1).filestr =  fileall{get(handles.popupmenu_file,'value')};
            PARCOMP(id+1).module  =get(handles.popupmenu_module, 'value');
            moduleall = get(handles.popupmenu_module,'string');
            PARCOMP(id+1).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
            PARCOMP(id+1).listgood =  listgood;
            PARCOMP(id+1).Factors =  PMI{currentsub}.tmpICA.Factors;
            PARCOMP(id+1).indt = indt; %indice de temps.
            PARCOMP(id+1).data = PMI{currentsub}.tmpICA.spar ;    
            PARCOMP(id+1).ComponentToKeep = PMI{currentsub}.tmpICA.selected;
            A =  PMI{currentsub}.tmpICA.Factors{1};
            C = PMI{currentsub}.tmpICA.Factors{3};
            Xm = C(:,PMI{currentsub}.tmpICA.selected)*A(PMI{currentsub}.tmpICA.selected,:);
            labelid  = get(handles.edit_selectedlabel,'string');
            PARCOMP(id+1).Xm = Xm;
            PARCOMP(id+1).label= [labelid,'ICA', sprintf('%03.0f',size(PARCOMP,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
            PARCOMP(id+1).type = 'ICA';
            PARCOMP(id+1).topo =  C(:,PMI{currentsub}.tmpICA.selected);
            newfile = 1;
        end
        for i=1:numel(PARCOMP)
            if iscell( PARCOMP(i).label)
                CORRlist{i} = cat(2,PARCOMP(i).label{:});
            else
                CORRlist{i} = PARCOMP(i).label;
            end
        end
        set(handles.listbox_Component,'string',CORRlist);
    else
        msgbox('Press the run button before the component could be add to the list');
        return
    end
    set(handles.listbox_Component,'string',CORRlist);
  
elseif strcmp(listmethod{idval},'GLM')   %extract GLM Component
    if isfield (PMI{currentsub}, 'tmpGLM')      
        indt = [PMI{currentsub}.tmpGLM.indt(1):PMI{currentsub}.tmpGLM.indt(2)];%Time indice
        listgood = PMI{currentsub}.tmpGLM.listgood;
        lstSV = PMI{1}.tmpGLM.selected;
        beta = PMI{currentsub}.tmpGLM.beta ;
        selected = PMI{currentsub}.tmpGLM.selected;
        idreg = PMI{currentsub}.tmpGLM.idreg;
        Xm = PMI{1}.tmpGLM.AUX.data{idreg(selected)}*beta(selected,:);
        label =  PMI{1}.tmpGLM.AUX.label{idreg(selected)};
    
        %save in a structure all apply correction
        [pathstr, name, ext] = fileparts(handles.NIRSpath{1});
        try
            load(fullfile(pathstr,'SelectedFactors.mat'));
            newfile = 0;
        catch
            %donot exist create the stucture
            PARCOMP.file= get(handles.popupmenu_file,'value');
            fileall = get(handles.popupmenu_file,'string');
            PARCOMP.filestr =  fileall{get(handles.popupmenu_file,'value')};
            PARCOMP.module  =get(handles.popupmenu_module, 'value');
            moduleall = get(handles.popupmenu_module,'string');
            PARCOMP.modulestr = moduleall{get(handles.popupmenu_module, 'value')};
            PARCOMP.listgood =  listgood;
            PARCOMP.beta =  PMI{1}.tmpGLM.beta;
            PARCOMP.std =  PMI{1}.tmpGLM.std;
            PARCOMP.AUX = PMI{1}.tmpGLM.AUX;
            PARCOMP.indt = indt; %indice de temps.
            PARCOMP.data = d(indt,listgood) ;%data(:,listgood,:);
            PARCOMP.Xm = Xm;
            PARCOMP.ComponentToKeep = PMI{1}.tmpGLM.selected;
            PARCOMP.idreg = PMI{1}.tmpGLM.idreg;
            labelid  = get(handles.edit_selectedlabel,'string');
            PARCOMP.label= [labelid,'GLM',label , sprintf('%03.0f',size(PARCOMP,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
            PARCOMP.type = 'GLM';
            PARCOMP.topo =  beta(selected,:);
            
            newfile = 1;
        end
        if newfile == 0
            id = numel(PARCOMP);
            %donot exist create the stucture
            PARCOMP(id+1).file= get(handles.popupmenu_file,'value');
            fileall = get(handles.popupmenu_file,'string');
            PARCOMP(id+1).filestr =  fileall{get(handles.popupmenu_file,'value')};
            PARCOMP(id+1).module  =get(handles.popupmenu_module, 'value');
            moduleall = get(handles.popupmenu_module,'string');
            PARCOMP(id+1).modulestr = moduleall{get(handles.popupmenu_module, 'value')};
            PARCOMP(id+1).listgood =  listgood;
            PARCOMP(id+1).beta =  PMI{currentsub}.tmpGLM.beta;
            PARCOMP(id+1).std =  PMI{currentsub}.tmpGLM.std;
            PARCOMP(id+1).AUX = PMI{currentsub}.tmpGLM.AUX;
            PARCOMP(id+1).indt = indt; %indice de temps.
            PARCOMP(id+1).data = d(indt,listgood) ;%data(:,listgood,:);
            PARCOMP(id+1).Xm = Xm;
            PARCOMP(id+1).ComponentToKeep = PMI{1}.tmpGLM.selected;
            PARCOMP(id+1).idreg = PMI{1}.tmpGLM.idreg;
            labelid  = get(handles.edit_selectedlabel,'string');
            PARCOMP(id+1).label= [labelid,'GLM',label , sprintf('%03.0f',size(PARCOMP,2)),' ',fileall{get(handles.popupmenu_file,'value')}];
            PARCOMP(id+1).type = 'GLM';
            PARCOMP(id+1).topo =  beta(selected,:);
            newfile = 1;
        end
        for i=1:numel(PARCOMP)
            if iscell( PARCOMP(i).label)
                CORRlist{i} = cat(2,PARCOMP(i).label{:});
            else
                CORRlist{i} = PARCOMP(i).label;
            end
        end
        set(handles.listbox_Component,'string',CORRlist);
    else
        msgbox('Press the run button before the component could be add to the list')
        return
    end
    
elseif strcmp(listmethod{idval},'Offset Adjustment')
    msgbox('No offset ajustement could be add in the component list');
    return
end


% --- Executes on button press in radiobutton23.
function radiobutton23_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton23



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_selectedlabel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selectedlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selectedlabel as text
%        str2double(get(hObject,'String')) returns contents of edit_selectedlabel as a double


% --- Executes during object creation, after setting all properties.
function edit_selectedlabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selectedlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_Detrend.
function btn_Detrend_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Detrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.radio_Detrend,'value',1)


guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart = 1
    end
end

if isempty(time_stop)
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
    
end
indstart = indstart(end);
indstop = indstop(end);
% if isempty(time_start)|isempty(time_stop)
%     set(handles.text_Advice,'string','Please enter time start and time stop to define time window')
%     return
% end

gui_DetrendIO(indstart,indstop);


% --- Executes on selection change in popupmethodselected.
function popupmethodselected_Callback(hObject, eventdata, handles)
% hObject    handle to popupmethodselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmethodselected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmethodselected


% --- Executes during object creation, after setting all properties.
function popupmethodselected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmethodselected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_run.
function btn_run_Callback(hObject, eventdata, handles)
% hObject    handle to btn_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Change individual bouton in popup menu
idval= get(handles.popupmethodselected,'val');
listmethod = get(handles.popupmethodselected,'string');
if strcmp(listmethod{idval},'Parafac') %get(handles.popupmethodselected,'value')==1 %PARAFAC
    btn_PARAFAC_Callback(hObject, eventdata, handles)
elseif strcmp(listmethod{idval},'PCA') %get(handles.popupmethodselected,'value')==2 %PCA
    btn_PCA_Callback(hObject, eventdata, handles)
elseif strcmp(listmethod{idval},'GLM')
    btn_GLM_Callback(hObject, eventdata, handles)
elseif strcmp(listmethod{idval},'ICA') %ICA
    btn_ICA_Callback(hObject, eventdata, handles)
    % elseif get(handles.popupmethodselected,'value')==4 %Wavelet
    %
    % elseif get(handles.popupmethodselected,'value')==5 %Linear detrend
    %     btn_Detrend_Callback(hObject, eventdata, handles)
    % elseif get(handles.popupmethodselected,'value')==6 %spline
    %
    % elseif get(handles.popupmethodselected,'value')==9 %PARAFAC WAVELET
    %     btn_PARAFACWAVELET_Callback(hObject, eventdata, handles)
    
end
function  btn_GLM_Callback(hObject, eventdata, handles)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart = 1;
    end
end

if isempty(time_stop);
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
    
end
indstart = indstart(end);
indstop = indstop(end);
% if isempty(time_start)|isempty(time_stop)
%     set(handles.text_Advice,'string','Please enter time start and time stop to define time window')
%     return
% end
GUI_GLMidentification(indstart,indstop);

function btn_PARAFACWAVELET_Callback(hObject, eventdata, handles)
%
%set(handles.radio_Parafac,'value',1)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
if isempty(time_start); indstart = 1;else;
    indstart= find(time_start>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstart)
        indstart = 1;
    end
end

if isempty(time_stop);
    indstop =numel(PMI{currentsub}.data(cf).HRF.tHRF);
else
    indstop= find(time_stop>PMI{currentsub}.data(cf).HRF.tHRF);
    if isempty(indstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    end
    
end
indstart = indstart(end);
indstop = indstop(end);
% if isempty(time_start)|isempty(time_stop)
%     set(handles.text_Advice,'string','Please enter time start and time stop to define time window')
%     return
% end

gui_ParafacWaveletIO(indstart,indstop);


%


% --- Executes on button press in btn_deleteComponent.
function btn_deleteComponent_Callback(hObject, eventdata, handles)
% hObject    handle to btn_deleteComponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[pathstr, name, ext] = fileparts(handles.NIRSpath{1});
try
    load(fullfile(pathstr,'SelectedFactors.mat'))
catch
    msgbox('Empty component, could not be deleted')
    return
end
id = get(handles.listbox_Component,'value');

%update list CorrectionApply
if numel(PARCOMP)>1
    PARCOMP(id) = [];
    save(fullfile(pathstr,'SelectedFactors.mat'),'PARCOMP');
    for i=1:numel(PARCOMP)
        CORRlist{i} = PARCOMP(i).label;
    end
    set(handles.listbox_Component,'string',CORRlist);
else
    delete(fullfile(pathstr,'SelectedFactors.mat'));
    set(handles.listbox_Component,'string',' ');
end
listbox_Component_Callback(hObject, eventdata, handles)


% --- Executes on button press in bnt_editzonecolor.
function bnt_editzonecolor_Callback(hObject, eventdata, handles)
% hObject    handle to bnt_editzonecolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
num = get(handles.popupmenu_zone,'value');
newcolor = uisetcolor([0,0,0]);
DOT{currentsub}.zone.color(num,:) =  newcolor;

set(handles.radio_colorzone,'value',1);
set(guiHOMER,'UserData',DOT);
radio_colorzone_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function btn_Restore_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btn_Restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in listbox_AUXCH.
function listbox_AUXCH_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_AUXCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_AUXCH contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_AUXCH


% --- Executes during object creation, after setting all properties.
function listbox_AUXCH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_AUXCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function menu_AUXlistview_Callback(hObject, eventdata, handles)
% hObject    handle to menu_AUXlistview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Auxlistview_select_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Auxlistview_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
AUXCH = get(handles.listbox_AUXCH,'string');
iAUXCH = get(handles.listbox_AUXCH,'value');
a = AUXCH{iAUXCH};
AUXCH{iAUXCH} = [a(1:end-4),'  on'];
PMI{currentsub}.data(cf).AUX.view{iAUXCH} = 1;
set(guiHOMER,'UserData',PMI);
set(handles.listbox_AUXCH,'string', AUXCH)
updatedisplay(handles);

% --------------------------------------------------------------------
function menu_Auxlistview_unselect_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Auxlistview_unselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

AUXCH = get(handles.listbox_AUXCH,'string');
iAUXCH = get(handles.listbox_AUXCH,'value');
a = AUXCH{iAUXCH};
AUXCH{iAUXCH} = [a(1:end-4),' off'];
global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
PMI{currentsub}.data(cf).AUX.view{iAUXCH} = 0;
set(guiHOMER,'UserData',PMI);
set(handles.listbox_AUXCH,'string', AUXCH)
updatedisplay(handles);


% --------------------------------------------------------------------
function menu_axesPhysiologie_newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to menu_axesPhysiologie_newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotAxes_physiologie(handles,1)
% --------------------------------------------------------------------
function menuAxes_Physiologie_Callback(hObject, eventdata, handles)
% hObject    handle to menuAxes_Physiologie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_SortComponent.
function btn_SortComponent_Callback(hObject, eventdata, handles)
% hObject    handle to btn_SortComponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[pathstr, name, ext] = fileparts(handles.NIRSpath{1});
load(fullfile(pathstr,'SelectedFactors.mat'));
for i=1:numel(PARCOMP)
    filelist(i) = PARCOMP(i).file;
    timelist(i) = PARCOMP(i).indt(1);
end
[val,id]=sort(filelist);
diff = [1, val(2:end) - val(1:end-1) ]; %find each similar file
idnewfile=find(diff);
idfinalorder = [];
for ival = 1:numel(idnewfile) %sort event subfile
    idfile = find(val(idnewfile(ival))==filelist);
    [valt,idsort] = sort(timelist(idfile));
    idfinalorder = [idfinalorder,idfile(idsort)];
end


%change order here
PARCOMP  = PARCOMP( idfinalorder);

for i=1:numel(PARCOMP)
    CORRlist{i} = PARCOMP(i).label;
end
set(handles.listbox_Component,'string',CORRlist);
save(fullfile(pathstr,'SelectedFactors.mat'),'PARCOMP');


% --- Executes on selection change in popupauxnb.
function popupauxnb_Callback(hObject, eventdata, handles)
% hObject    handle to popupauxnb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupauxnb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupauxnb
set(handles.listbox_AUXCH,'value',1)
updatedata(handles,1,1)
updatedisplay(handles)
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function popupauxnb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupauxnb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function context_listbox_Component_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_listbox_Component_newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Component_newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.context_listbox_Component_newfigure,'checked'),'on')
    set(handles.context_listbox_Component_newfigure,'checked','off')
else
    set(handles.context_listbox_Component_newfigure,'checked','on')
end
listbox_Component_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function context_triggercopy_Callback(hObject, eventdata, handles)
% hObject    handle to context_triggercopy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
idfile = get(handles.popupmenu_file,'value');
idmodulehold = get(handles.popupmenu_module_hold, 'value');
% if strcmp(handles.NIRS.Dt.fir.pp(idmodulehold).pre,'READ_RAW_NIRS')&get(handles.radiobutton_hold,'value')
%     selected_trig=handles.NIRS.Dt.fir.aux5ini{idfile};
% else
%     selected_trig = handles.selected_trig;
% end
dt = PMI{currentsub}.data(cf).HRF.tHRF(2)-PMI{currentsub}.data(cf).HRF.tHRF(1);
val = get(handles.popupmenu8,'value');
trigall = get(handles.popupmenu8,'string');
strtrig = trigall{val};
numtrig = str2num(strtrig(1:end-4));
selected_trig = handles.triggers(ismember(handles.triggers(:,1), numtrig),:);

data = selected_trig(:,2)*dt;
label = [];
for i = 1:numel(data)
    label = [label, sprintf('%f\n',data(i))];
end
clipboard('copy', label);
%
%     if get(handles.radio_triggers,'value')
%
%         for i = 1:size(selected_trig,1)
%             try
%             trig_ind = selected_trig(i,2);
%             if trig_ind <=numel(PMI{currentsub}.data(cf).HRF.tHRF)
%                 trig_ind = PMI{currentsub}.data(cf).HRF.tHRF(trig_ind);
%                 h = plot([trig_ind, trig_ind],[yli(1), yli(2)], 'g');
%             end
%             catch
%             end
%         end
%     end

% --------------------------------------------------------------------
function context_trigger_Callback(hObject, eventdata, handles)
% hObject    handle to context_trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_HbO.
function btn_HbO_Callback(hObject, eventdata, handles)
% hObject    handle to btn_HbO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = 1;

nch = size(PMI{currentsub}.data(cf).HRF.AvgC,2);
PMI{currentsub}.data(cf).MeasList(:,4)==1; % TODO: Verify if this is the intended operation as it looks like an assignment.
MeasListAct=PMI{currentsub}.data(cf).MeasListAct;
MeasListAct(:)=1;
plotLst= find( PMI{currentsub}.data(cf).MeasListAct & PMI{currentsub}.data(cf).MeasList(:,4)==1);
if ~isempty(plotLst)
    PMI{currentsub}.plotLst = plotLst;
    set(handles.edit_plotLst,'string',mat2str(plotLst))
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
    
end
handles.newlist=1;
set(guiHOMER,'UserData',PMI)
guidata(hObject, handles);
updatedisplay(handles)


% --- Executes on button press in btn_HBR.
function btn_HBR_Callback(hObject, eventdata, handles)
% hObject    handle to btn_HBR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = 1;

nch = size(PMI{currentsub}.data(cf).HRF.AvgC,2);
PMI{currentsub}.data(cf).MeasList(:,4)==1; % TODO: See if this is intented. Looks like bebug code.
MeasListAct=PMI{currentsub}.data(cf).MeasListAct;
MeasListAct(:)=1;
plotLst= find( PMI{currentsub}.data(cf).MeasListAct & PMI{currentsub}.data(cf).MeasList(:,4)==2);
if ~isempty(plotLst)
    PMI{currentsub}.plotLst = plotLst;
    set(handles.edit_plotLst,'string',mat2str(plotLst))
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
    
end
handles.newlist=1;
set(guiHOMER,'UserData',PMI)
guidata(hObject, handles);
updatedisplay(handles)


% --------------------------------------------------------------------
function context_module_deletedata_Callback(hObject, eventdata, handles)
% hObject    handle to context_module_deletedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idmodule = get(handles.popupmenu_module,'value');
rDtp = handles.NIRS.Dt.fir.pp(idmodule).p;
modulename = handles.NIRS.Dt.fir.pp(idmodule).pre;
button = questdlg(['Are you sure you want to delete file from', modulename],'Delete','yes','no','no');

switch button
    case 'yes'
        for f = 1:numel(rDtp)
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            deletenir  = fullfile(dir1,[fil1 '.nir']);
            deletevmrk = fullfile(dir1,[fil1 '.vmrk']);
            deletevhdr = fullfile(dir1,[fil1 '.vhdr']);
            delete(deletenir );
            delete(deletevmrk);
            delete(deletevhdr);
        end
    case 'no'
end


% --- Executes during object creation, after setting all properties.
function btn_substractPARAFAC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btn_substractPARAFAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in btn_setaux.
function btn_setaux_Callback(hObject, eventdata, handles)
% hObject    handle to btn_setaux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
currentsub=1;
cf = PMI{currentsub}.currentFile;
if isfield(PMI{currentsub}.data(cf), 'AUX') % Verify if the aux channel is available.
    PMI{currentsub}.data(cf).AUX;
    idfile=PMI{1}.idfile;
    time_start = str2num(get(handles.edit_time_start,'string')); %Debut segment
    time_stop = str2num(get(handles.edit_time_stop,'string'));  %Fin segment
    iAUX = get(handles.popupauxnb,'value');
    nameAUX = handles.NIRS.Dt.AUX(iAUX).pp(end).p{idfile};
    [pathtmp,filetmp,exttmp]=fileparts(nameAUX);
    [AUX.data,AUX.infoBV,AUX.marker,AUX.ind_dur_ch] = fopen_EEG(nameAUX);

    if 1==get(handles.btn_setaux,'value')
        pstart = round((PMI{currentsub}.data(cf).AUX.fs{1})*time_start);
        pstop = round((PMI{currentsub}.data(cf).AUX.fs{1})*time_stop);
        AUX.data(pstart:pstop,:)=0;
    end
    nb=size(AUX.data,1);
    fwrite_EEG(nameAUX,AUX(iAUX),1,nb);
else
    errordlg('The auxiliary channel could not be set: Aux not found.', 'Configuration error') % Inform user.
end

% --- Executes on selection change in popup_setaux.
function popup_setaux_Callback(hObject, eventdata, handles)
% hObject    handle to popup_setaux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_setaux contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_setaux


% --- Executes during object creation, after setting all properties.
function popup_setaux_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_setaux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function context_listbox_Component_Copy_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Component_Copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idComp = get(handles.listbox_Component,'value');
strComp = get(handles.listbox_Component,'string');
label = strComp{idComp};
clipboard('copy', label);


% --------------------------------------------------------------------
function menu_context_trigger_select_Callback(hObject, eventdata, handles)
% hObject    handle to menu_context_trigger_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

triggerinfo = get(handles.popupmenu8,'string');
itriggerinfo = get(handles.popupmenu8,'value');
a = triggerinfo{itriggerinfo};
triggerinfo{itriggerinfo} = [a(1:end-4),'  on'];
set(handles.popupmenu8,'string', triggerinfo);
updatedisplay(handles);


% --------------------------------------------------------------------
function menu_context_trigger_unselect_Callback(hObject, eventdata, handles)
% hObject    handle to menu_context_trigger_unselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
triggerinfo = get(handles.popupmenu8,'string');
itriggerinfo = get(handles.popupmenu8,'value');
a = triggerinfo{itriggerinfo};
triggerinfo{itriggerinfo} = [a(1:end-4),' off'];
global currentsub;
set(handles.popupmenu8,'string', triggerinfo)
updatedisplay(handles);


% --------------------------------------------------------------------
function remove_bad_ch_Unselect_Callback(hObject, eventdata, handles)
% hObject    handle to remove_bad_ch_Unselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
plotLstrm = str2double(get(gco, 'Tag'));
id = find(PMI{currentsub}.plotLst==plotLstrm);

PMI{currentsub}.plotLst(id)=[];
PMI{currentsub}.plot(id,:)=[];
% = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
set(handles.edit_plotLst,'string',num2str(PMI{currentsub}.plotLst));
set(guiHOMER,'UserData',PMI);
handles.newlist=1;
guidata(hObject, handles);
updatedisplay(handles)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_exportmeanzonecurve.
function btn_exportmeanzonecurve_Callback(hObject, eventdata, handles)
% hObject    handle to btn_exportmeanzonecurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
 num = get(handles.popupmenu_zone,'value');
  plotLst = PMI{currentsub}.zone.plotLst{num};
   label = PMI{currentsub}.zone.label{num};
   try
     colorid = PMI{currentsub}.zone.color(num,:);
   catch
       colorid = [0 0 0]
   end
     listCH = zeros(size( PMI{currentsub}.data.HRF.AvgC,2),1);
   listCH(plotLst) = 1;
  both =  reshape(listCH,numel(listCH)/2,2);
  chHbO = find(sum(both,2));
  chHbR =   chHbO + numel(listCH)/2;
  
 if strcmp(get(handles.Context_MenuExportZone_All,'checked'),'off') %only one curve
      figure; hold on
      plotLst = PMI{currentsub}.zone.plotLst{num};
      label = PMI{currentsub}.zone.label{num};
      colorid = PMI{currentsub}.zone.color(num,:);
      listCH = zeros(size( PMI{currentsub}.data.HRF.AvgC,2),1);
      listCH(plotLst) = 1;
      both =  reshape(listCH,numel(listCH)/2,2);
      chHbO = find(sum(both,2));
      chHbR =   chHbO + numel(listCH)/2;
      % plot(PMI{currentsub}.data.HRF.AvgC(:,chHbO))
      plot(PMI{currentsub}.data.HRF.tHRF,nanmean(PMI{currentsub}.data.HRF.AvgC(:,chHbO),2),'displayname',['HbO',label],'color',colorid,'linewidth',4);
      plot(PMI{currentsub}.data.HRF.tHRF,nanmean(PMI{currentsub}.data.HRF.AvgC(:,chHbR),2),'displayname',['HbR',label],'color',colorid,'linestyle','--','linewidth',4);
      set(gca,'fontsize',20)
      xlabel('Time (s)')
 else
     figure; hold on
     for num = 1:numel(PMI{currentsub}.zone.label)           
          plotLst = PMI{currentsub}.zone.plotLst{num};
          label = PMI{currentsub}.zone.label{num};
          colorid = PMI{currentsub}.zone.color(num,:);
          listCH = zeros(size( PMI{currentsub}.data.HRF.AvgC,2),1);
          listCH(plotLst) = 1;
          both =  reshape(listCH,numel(listCH)/2,2);
          chHbO = find(sum(both,2));
          chHbR =   chHbO + numel(listCH)/2;
          % plot(PMI{currentsub}.data.HRF.AvgC(:,chHbO))
          plot(PMI{currentsub}.data.HRF.tHRF,nanmean(PMI{currentsub}.data.HRF.AvgC(:,chHbO),2),'displayname',['HbO',label],'color',colorid,'linewidth',4);
          plot(PMI{currentsub}.data.HRF.tHRF,nanmean(PMI{currentsub}.data.HRF.AvgC(:,chHbR),2),'displayname',['HbR',label],'color',colorid,'linestyle','--','linewidth',4);
          set(gca,'fontsize',20)
          xlabel('Time (s)')
     end
 end
tstart = str2num(get(handles.edit_start,'string'));
duration =  str2num(get(handles.edit_duration,'string'));
if tstart==0
else
    xlim([tstart,tstart+duration]);
end
    
function edit_videoname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_videoname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_videoname as text
%        str2double(get(hObject,'String')) returns contents of edit_videoname as a double


% --- Executes during object creation, after setting all properties.
function edit_videoname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_videoname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_AUDIO.
function btn_AUDIO_Callback(hObject, eventdata, handles)
% hObject    handle to btn_AUDIO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
time_start = str2num(get(handles.edit_time_start,'string'));
time_stop = str2num(get(handles.edit_time_stop,'string')); 
id = get(handles.popupmenu_file,'value');
filenameaudio = handles.NIRS.Dt.Audio.pp(end).p{id};
offset =  handles.NIRS.Dt.Audio.pp(end).sync_timesec{id};
try
  [h,m,s] = hms(seconds(offset+time_start))
   set(handles.text_Advice,'string',['Time ',sprintf('%02.0f',h),':',sprintf('%02.0f',m),':',sprintf('%02.0f',s),' ', filenameaudio])
catch
   set(handles.text_Advice,'string',['Time ', num2str(round(offset+time_start)),'s ', filenameaudio])
end
 if isempty(offset)
        msgbox('Segment first to synchronised video with NIRS')
        return
 end   
  if isempty(time_start)|isempty(time_stop)
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        msgbox({'\fontsize{14} Enter time start and time stop'}, '','help', CreateStruct)
        return
  end
   if (  time_stop- time_start  )>60
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        msgbox({'\fontsize{14} Please select an interval time' ; 'start and stop shorter than 30 seconds'}, '','help', CreateStruct)
        return
   end
    
   info = audioinfo(filenameaudio);
   [y,Fs] = audioread(filenameaudio,round(([time_start time_stop] +offset)*info.SampleRate));
   sound(y,Fs)
   axes(handles.display_axes)
   for i=time_start:0.25:time_stop
         t = (i);
         if i==time_start
             hline = plot([time_start time_start],[0,1],'k')
         else
            set(hline, 'xdata',[t t])
         end
         pause(0.25)
   end
  delete(hline);
   


% --- Executes on button press in btn_zoom_timestart_timestop.
function btn_zoom_timestart_timestop_Callback(hObject, eventdata, handles)
% hObject    handle to btn_zoom_timestart_timestop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tstart = str2num(get(handles.edit_time_start,'string'));
tstop = str2num(get(handles.edit_time_stop,'string'));
set(handles.edit_start,'string', num2str(tstart-3));
set(handles.edit_duration,'string', num2str(tstop-tstart+6));
updatedisplay(handles)


% --------------------------------------------------------------------
function menu_settings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_settings_video_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

optionVideo


% --- Executes on selection change in listbox_normlist.
function listbox_artefact_identify_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_normlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_normlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_normlist


% --- Executes during object creation, after setting all properties.
function listbox_artefact_identify_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_normlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_autozoom.
function radiobutton_autozoom_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_autozoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_autozoom


% --------------------------------------------------------------------
function context_autozoom_Callback(hObject, eventdata, handles)
% hObject    handle to context_autozoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_enable_autozoom_Callback(hObject, eventdata, handles)
% hObject    handle to context_enable_autozoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if  strcmp(get(handles.context_enable_autozoom,'Checked'),'off')	
    set(handles.context_enable_autozoom,'Checked','on')
elseif strcmp(get(handles.context_enable_autozoom,'Checked'),'on')
     set(handles.context_enable_autozoom,'Checked','off')
end


% --- Executes on button press in btn_sort_correction.
function btn_sort_correction_Callback(hObject, eventdata, handles)
% hObject    handle to btn_sort_correction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[pathstr, name, ext] = fileparts(handles.NIRSpath{1});
load(fullfile(pathstr,'CorrectionApply.mat'));
for i=1:numel(PARCORR)
    filelist(i) = PARCORR(i).file;
    timelist(i) = PARCORR(i).indt(1);
end
[val,id]=sort(filelist);
diff = [1, val(2:end) - val(1:end-1) ]; %find each similar file
idnewfile=find(diff);
idfinalorder = [];
for ival = 1:numel(idnewfile) %sort event subfile
    idfile = find(val(idnewfile(ival))==filelist);
    [valt,idsort] = sort(timelist(idfile));
    idfinalorder = [idfinalorder,idfile(idsort)];
end


%change order here
PARCORR  = PARCORR( idfinalorder);

for i=1:numel(PARCORR)
    CORRlist{i} = PARCORR(i).label;
end
set(handles.listbox_CorrectionDecomposition,'string',CORRlist);
save(fullfile(pathstr,'CorrectionApply.mat'),'PARCORR');


% --- Executes on button press in btn_zoomout.
function btn_zoomout_Callback(hObject, eventdata, handles)
% hObject    handle to btn_zoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tstart = str2num(get(handles.edit_start,'string'));
duration = str2num(get(handles.edit_duration,'string'));
set(handles.edit_start,'string', num2str(tstart-3));
set(handles.edit_duration,'string', num2str(duration+6));
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
tstop = tstart + duration
PMI = get(guiHOMER,'UserData');
if tstart<0 | tstop > PMI{1}.data.HRF.tHRF(end);
    
set(handles.edit_start,'string', '0');
set(handles.edit_duration,'string', '1');
end
updatedisplay(handles)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function context_listbox_Component_rename_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Component_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
label = get(handles.listbox_Component,'string');
id = get(handles.listbox_Component,'value');
msgbox('FINISH CODE');


% --------------------------------------------------------------------
function Context_MenuExportZone_Callback(hObject, eventdata, handles)
% hObject    handle to Context_MenuExportZone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Context_MenuExportZone_All_Callback(hObject, eventdata, handles)
% hObject    handle to Context_MenuExportZone_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.Context_MenuExportZone_All,'checked'),'on')
    set(handles.Context_MenuExportZone_All,'checked','off')
else
  set(handles.Context_MenuExportZone_All,'checked','on')
end


% --- Executes on button press in radio_multiplezone.
function radio_multiplezone_Callback(hObject, eventdata, handles)
% hObject    handle to radio_multiplezone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_multiplezone


% --------------------------------------------------------------------
function context_listbox_Correction_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Correction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_listbox_Correction_newfigure_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Correction_newfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.context_listbox_Correction_newfigure,'checked'),'on')
    set(handles.context_listbox_Correction_newfigure,'checked','off')
else
    set(handles.context_listbox_Correction_newfigure,'checked','on')
end
listbox_CorrectionDecomposition_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function context_listbox_Correction_Copy_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Correction_Copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idComp = get(handles.listbox_CorrectionDecomposition,'value');
strComp = get(handles.listbox_CorrectionDecomposition,'string');
label = strComp{idComp};
clipboard('copy', label);

% --------------------------------------------------------------------
function context_listbox_Correction_rename_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Correction_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
label = get(handles.listbox_CorrectionDecomposition,'string');
id = get(handles.listbox_CorrectionDecomposition,'value');
msgbox('FINISH CODE');


% --------------------------------------------------------------------
function context_addnewtrigger_Callback(hObject, eventdata, handles)
% hObject    handle to context_addnewtrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
position = get(gca,  'CurrentPoint')
idfile = get(handles.popupmenu_file,'value');
aux5 = handles.NIRS.Dt.fir.aux5{idfile};
newsample = round(handles.NIRS.Cf.dev.fs*position(1))
handles.NIRS.Dt.fir.aux5{idfile} = [aux5;[1, newsample]]
defaulttrig = '1';
defaulttrig = EnterANewTrigger(defaulttrig);
trigval = str2num(defaulttrig);

newsample = round(handles.NIRS.Cf.dev.fs*position(1))
handles.NIRS.Dt.fir.aux5{idfile} = [aux5;[trigval, newsample]]

%add the trig in the list
if isfield(handles.NIRS.Dt.fir,'aux5')
    handles.triggers = handles.NIRS.Dt.fir.aux5{idfile};
    
    unique_trig = unique(handles.triggers(:,1));
    %for i=1:numel(unique_trig)
    trig_list = [];
    for i=1:numel(unique_trig)
        trig_list = [trig_list;  {[num2str(unique_trig(i)),'  on']}];
    end
    set(handles.popupmenu8, 'String', trig_list);
    set(handles.popupmenu8, 'value', numel(trig_list));
    
end


NIRS = handles.NIRS;
NIRSpath = get(handles.edit_nirsmat,'string');
subjectnb = get(handles.edit_nirsmat,'value');
save(NIRSpath{subjectnb,1},'NIRS','-mat');
guidata(hObject, handles);
updatedisplay(handles)


% --------------------------------------------------------------------
function Context_removetrigger_Callback(hObject, eventdata, handles)
% hObject    handle to Context_removetrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

position = get(gco,  'xdata');
idfile = get(handles.popupmenu_file,'value');
aux5 = handles.NIRS.Dt.fir.aux5{idfile};
idselect = find(aux5(:,2)==round(handles.NIRS.Cf.dev.fs*position(1)));
trigval = aux5(idselect ,1);
onlytrigval = find(aux5(:,1)==trigval);
nbtrig = find(onlytrigval == idselect);
set(handles.Context_infotrig,'text',['Trig ',num2str(trigval),' Position=',num2str(nbtrig),'st']);

% --------------------------------------------------------------------
function Context_removetrig_Callback(hObject, eventdata, handles)
% hObject    handle to Context_removetrig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1

position = get(gco,  'xdata')
idfile = get(handles.popupmenu_file,'value');
aux5 = handles.NIRS.Dt.fir.aux5{idfile};
idremove = find(aux5(:,2)==round(handles.NIRS.Cf.dev.fs*position(1)))
aux5(idremove,:)=[];
handles.NIRS.Dt.fir.aux5{idfile} = aux5

%add the trig in the list
if isfield(handles.NIRS.Dt.fir,'aux5')
    handles.triggers = handles.NIRS.Dt.fir.aux5{idfile};
    
    unique_trig = unique(handles.triggers(:,1));
    %for i=1:numel(unique_trig)
    trig_list = [];
    for i=1:numel(unique_trig)
        trig_list = [trig_list;  {[num2str(unique_trig(i)),'  on']}];
    end
    set(handles.popupmenu8, 'String', trig_list);
    set(handles.popupmenu8, 'value', numel(trig_list));
    
end


NIRS = handles.NIRS;
NIRSpath = get(handles.edit_nirsmat,'string');
subjectnb = get(handles.edit_nirsmat,'value');
save(NIRSpath{subjectnb,1},'NIRS','-mat');
guidata(hObject, handles);
updatedisplay(handles)


% --------------------------------------------------------------------
function menu_settings_trigger_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings_trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_trigger_color_Callback(hObject, eventdata, handles)
% hObject    handle to context_trigger_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%c = uisetcolor
idfile = get(handles.popupmenu_file,'value');
triggerinfo = get(handles.popupmenu8,'string');
itriggerinfo = get(handles.popupmenu8,'value');
a = triggerinfo{itriggerinfo};
colorRGV = uisetcolor();

eval(['handles.NIRS.Dt.fir.auxtrig', num2str(a(1:end-4)),'.color = [',num2str(colorRGV(1)),',',num2str(colorRGV(2)),',',num2str(colorRGV(3)),'];']);


% a(1:end-4)
% aux5 = handles.NIRS.Dt.fir.aux5{idfile};
% unique_trig = unique(aux5(:,1)) 
% for i=1:numel(unique_trig)
%     eval(['handles.NIRS.Dt.fir.auxtrig', num2str(a(1:end-4)),'.label = ''',num2str(a(1:end-4)),''';']);
%     colorRGV = [0, 1 , 0]
%     eval(['handles.NIRS.Dt.fir.auxtrig', num2str(a(1:end-4)),'.color = [',num2str(colorRGV(1)),',',num2str(colorRGV(2)),',',num2str(colorRGV(3)),'];']);
% end
NIRS = handles.NIRS;
NIRSpath = get(handles.edit_nirsmat,'string');
subjectnb = get(handles.edit_nirsmat,'value');

save(NIRSpath{subjectnb,1},'NIRS','-mat');
guidata(hObject, handles);
updatedisplay(handles)
% --------------------------------------------------------------------
function context_trigger_label_Callback(hObject, eventdata, handles)
% hObject    handle to context_trigger_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idfile = get(handles.popupmenu_file,'value');

triggerinfo = get(handles.popupmenu8,'string');
itriggerinfo = get(handles.popupmenu8,'value');
a = triggerinfo{itriggerinfo};
newname = Zonename();
 eval(['handles.NIRS.Dt.fir.auxtrig', num2str(a(1:end-4)),'.label = ''',newname,''';']);
NIRS = handles.NIRS;
NIRSpath = get(handles.edit_nirsmat,'string');
subjectnb = get(handles.edit_nirsmat,'value');

save(NIRSpath{subjectnb,1},'NIRS','-mat');
guidata(hObject, handles);
updatedisplay(handles)
% --------------------------------------------------------------------
function menu_multimodalfiles_Callback(hObject, eventdata, handles)
% hObject    handle to menu_multimodalfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NIRSpath = get(handles.edit_nirsmat,'string');
h = ModifyMultimodalName(NIRSpath);
waitfor(h);

 load(NIRSpath{1}) ;
 handles.NIRS = NIRS;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Context_Zone_Rename_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Zone_Rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
num = get(handles.popupmenu_zone,'value');

newname = Zonename(DOT{currentsub}.zone.label(num));
DOT{currentsub}.zone.label(num) = newname;
set(handles.popupmenu_zone,'string',DOT{currentsub}.zone.label);
set(guiHOMER,'UserData',DOT);

% --------------------------------------------------------------------
function Context_Zone_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Context_Zone_ViewChOrder_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Zone_ViewChOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
num = get(handles.popupmenu_zone,'value');
zonelabel.plotlst = num2str(DOT{currentsub}.zone.plotLst{num});
newname = Zoneplotlst(zonelabel);
DOT{currentsub}.zone.plotLst{num} =  str2num(newname);
set(handles.popupmenu_zone,'string',DOT{currentsub}.zone.label);
set(guiHOMER,'UserData',DOT);

% --------------------------------------------------------------------
function Context_Zone_EditColor_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Zone_EditColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentsub;
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
DOT = get(guiHOMER,'UserData');
cf = DOT{currentsub}.currentFile;
num = get(handles.popupmenu_zone,'value');
newcolor = uisetcolor([0,0,0]);
DOT{currentsub}.zone.color(num,:) = newcolor;

set(handles.radio_colorzone,'value',1);
set(guiHOMER,'UserData',DOT);
radio_colorzone_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function menu_Auxlistview_copy_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Auxlistview_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


AUXCH = get(handles.listbox_AUXCH,'string');
iAUXCH = get(handles.listbox_AUXCH,'value');
a = AUXCH{iAUXCH};
clipboard('copy', a);


% --------------------------------------------------------------------
function Context_infotrig_Callback(hObject, eventdata, handles)
% hObject    handle to Context_infotrig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_clearartifact_Callback(hObject, eventdata, handles)
% hObject    handle to context_clearartifact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    currentsub=1;
    PMI = get(guiHOMER,'UserData');
    cf = PMI{currentsub}.currentFile;
    NChalf = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
    plotLst = PMI{currentsub}.plotLst;
    plotLstAll = [];
    for i = 1:numel(plotLst)
        chlist = plotLst(i);
        if isempty(find(chlist - NChalf < 1 ))
            plotLstAll = [plotLstAll,chlist, chlist-NChalf]; %add 830 nm channels
        else
            plotLstAll = [plotLstAll,chlist,chlist+NChalf]; %add 690 nm channels
        end
    end
    plotLst = plotLstAll;
    posxstart = str2num(get(handles.edit_time_start,'string'));
    posxstop = str2num(get(handles.edit_time_stop,'string'));
    if isempty(posxstart)
        indstart = 1;
    else
        indstart = find(posxstart>PMI{currentsub}.data(cf).HRF.tHRF);
    end
    if isempty(posxstop)
        indstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    else
        indstop= find(posxstop>PMI{currentsub}.data(cf).HRF.tHRF);
    end
    PMI{currentsub}.data(cf).HRF.noise(:,plotLst)= 0;
    set(guiHOMER,'UserData',PMI);
    updatedisplay(handles);


% --- Executes on button press in radio_3dmontageupdate.
function radio_3dmontageupdate_Callback(hObject, eventdata, handles)
% hObject    handle to radio_3dmontageupdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_3dmontageupdate


% --------------------------------------------------------------------
function context_listbox_Correction_CopyData_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Correction_CopyData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentsub;

idComp = get(handles.listbox_CorrectionDecomposition,'value');
strComp = get(handles.listbox_CorrectionDecomposition,'string');
label = strComp{idComp};
[pathstr, ~, ~] = fileparts(handles.NIRSpath{1});
strTopo = sprintf('%s\t %s\t %s\n', 'Detector', 'Source', label);

try
    load(fullfile(pathstr,'CorrectionApply.mat'))
    load(fullfile(pathstr, 'NIRS.mat'));
catch
    return
end
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
ML = PMI{currentsub}.data.MeasList;
topo = PARCORR(idComp).topo;


% if strcmpi((PARCOMP(idComp).type),'GLM')
    ntopo = numel(topo);
% elseif strcmpi((PARCOMP(idComp).type),'Parafac')
%     ntopo = ;
% elseif strcmpi((PARCOMP(idComp).type),'PCA')
%     ntopo = ;
% end
if ntopo == 0
    msgbox('There is no data associated to this label')
    return
end
for j = 1:ntopo
    ichannel = PARCORR(idComp).listgood(j);
    switch  NIRS.Cf.dev.n
        case 'ISS Imagent'
            strDet = SDDet2strboxy_ISS(ML(ichannel,2));
            strSrs = SDPairs2strboxy_ISS(ML(ichannel,1));
                     
        case 'NIRx'                          
            strDet = SDDet2strboxy(ML(ichannel,2));
            strSrs = SDPairs2strboxy(ML(ichannel,1));
                                                                 
        otherwise                            
            strDet = SDDet2strboxy_ISS(ML(ichannel,2));
            strSrs = SDPairs2strboxy_ISS(ML(ichannel,1));
                       
    end
    strTopo = [strTopo,sprintf('%s\t %s\t %0.5g\n', strDet, strSrs, topo(j))]; 
end

clipboard('copy', strTopo);

msgbox('The data of the selected component is copied in the clipboard');




% --------------------------------------------------------------------
function context_listbox_Component_CopyData_Callback(hObject, eventdata, handles)
% hObject    handle to context_listbox_Component_CopyData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentsub;

idComp = get(handles.listbox_Component,'value');
strComp = get(handles.listbox_Component,'string');
label = strComp{idComp};
strTopo = sprintf('%s\t %s\t %s\n', 'Detector', 'Source', label);
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
PMI = get(guiHOMER,'UserData');
[pathstr, ~ , ~ ] = fileparts(handles.NIRSpath{1});
try
    load(fullfile(pathstr,'SelectedFactors.mat'));
    load(fullfile(pathstr, 'NIRS.mat'));
catch
%     prelabel = get(handles.edit_selectedlabel, 'string');
%     if strncmpi((label),[prelabel,'PARAFAC'],numel(prelabel)+numel('PARAFAC'))
%         FacSpatial = PMI{currentsub}.tmpPARAFAC.Factors{2};
%         selected = PMI{currentsub}.tmpPARAFAC.selected;
%         PARCOMP.topo = FacSpatial(:,selected);
%         listgood = PMI{currentsub}.tmpPARAFAC.listgood;
%         PARCOMP.listgood = listgood;
%         
%     elseif strncmpi((label),[prelabel,'PCA'],numel(prelabel)+numel('PCA'))
%         listgood = PMI{currentsub}.tmpPCA.listgood;
%         u =PMI{currentsub}.tmpPCA.u ;
%         s =PMI{currentsub}.tmpPCA.s ;
%         v = PMI{currentsub}.tmpPCA.v ;
%         lstSV = PMI{1}.tmpPCA.selected;
%         PARCOMP.topo = s(lstSV,lstSV)*v(:,lstSV)';
%         PARCOMP.listgood = listgood;
%         
%     elseif strncmpi((label),[prelabel,'GLM'],numel(prelabel)+numel('GLM'))
%         beta = PMI{currentsub}.tmpGLM.beta ;
%         selected = PMI{currentsub}.tmpGLM.selected;
%         listgood = PMI{currentsub}.tmpGLM.listgood;
%         PARCOMP.topo =  beta(selected,:);
%         PARCOMP.listgood = listgood;
%        
%     end
%     load(fullfile(pathstr, 'NIRS.mat'));
    msgbox('The data couldn''t be copied.');
    return
end

ML = PMI{currentsub}.data.MeasList;
topo = PARCOMP(idComp).topo;


% if strcmpi((PARCOMP(idComp).type),'GLM')
    ntopo = numel(topo);
% elseif strcmpi((PARCOMP(idComp).type),'Parafac')
%     ntopo = ;
% elseif strcmpi((PARCOMP(idComp).type),'PCA')
%     ntopo = ;
% end

if ntopo == 0
    msgbox('There is no data associated to this label')
    return
end

for j = 1:ntopo
    ichannel = PARCOMP(idComp).listgood(j);
    switch  NIRS.Cf.dev.n
        case 'ISS Imagent'
            strDet = SDDet2strboxy_ISS(ML(ichannel,2));
            strSrs = SDPairs2strboxy_ISS(ML(ichannel,1));
                     
        case 'NIRx'                          
            strDet = SDDet2strboxy(ML(ichannel,2));
            strSrs = SDPairs2strboxy(ML(ichannel,1));
                                                                 
        otherwise                            
            strDet = SDDet2strboxy_ISS(ML(ichannel,2));
            strSrs = SDPairs2strboxy_ISS(ML(ichannel,1));
                       
    end
    strTopo = [strTopo,sprintf('%s\t %s\t %0.5g\n', strDet, strSrs, topo(j))]; 
end

clipboard('copy', strTopo);

msgbox('The data of the selected component is copied in the clipboard');


% --------------------------------------------------------------------
function context_reportnoise_Callback(hObject, eventdata, handles)
% hObject    handle to context_reportnoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[pathstr, ~ , ~ ] = fileparts(handles.NIRSpath{1});
try
    load(fullfile(pathstr, 'NIRS.mat'));
catch
    return
end

guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
NC = NIRS.Cf.H.C.N;  

for f=1:size(rDtp,1)
    d = fopen_NIR(rDtp{f,1},NC);           
    samp_length = size(d,2);
    time = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:samp_length*1/NIRS.Cf.dev.fs;
    dim_time = size(time);
    dim_noise = size(PMI{currentsub}.data(1).HRF.noise);
    perc_artifactedChannel = sum((PMI{currentsub}.data(1).HRF.noise')/size(PMI{currentsub}.data(1).HRF.noise',1)*100);
    perc_artifactPerChannel = sum((PMI{currentsub}.data(1).HRF.noise)/size(PMI{currentsub}.data(1).HRF.noise,1)*100);
    
    if sum((PMI{currentsub}.data(1).HRF.noise)) == 0
        msgbox('The data is not artifacted. The noise report may be insignificant', 'Noise Report');
    end
    
    figure;
    subplot(5,5,[1:3,6:8,11:13]);
    imagesc(time,1:size(d,1),PMI{currentsub}.data(1).HRF.noise');
    ylabel('Channel id number','fontsize', 8);
    title('Channels in function of time','fontsize', 10);
    set(gca,'fontsize',8);
    xlim([0,time(dim_time(2))]);

    subplot(5,5,[4:5,9:10,14:15]);
    channel = linspace(dim_noise(2),1,dim_noise(2));
    plot(perc_artifactPerChannel,channel);
    set(gca,'Ydir','reverse');
    set(gca,'fontsize',8);
    xlabel('Percentage of artifacted time (%)','fontsize',8);    
    title('Channels in function of the percentage of artifacted time','fontsize', 10);
    ylim([1,dim_noise(2)]);
    xlim([0,100]);

    subplot(5,5,[16:18,21:23]);
    plot(time,perc_artifactedChannel);
    ylabel('Percentage of artifacted channels (%)','fontsize', 8);
    xlabel('Time (s)','fontsize', 8)
    title({'Percentage of artifacted channels in function of time'},'fontsize', 10);
    xlim([1,time(dim_time(2))]);
    ylim([0,100]);
    set(gca,'fontsize',8)
end
