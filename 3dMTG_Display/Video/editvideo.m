function varargout = editvideo(varargin)
% EDITVIDEO M-file for editvideo.fig
%      EDITVIDEO, by itself, creates a new EDITVIDEO or raises the existing
%      singleton*.
%
%      H = EDITVIDEO returns the handle to a new EDITVIDEO or the handle to
%      the existing singleton*.
%
%      EDITVIDEO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDITVIDEO.M with the given input arguments.
%
%      EDITVIDEO('Property','Value',...) creates a new EDITVIDEO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before editvideo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to editvideo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help editvideo

% Last Modified by GUIDE v2.5 02-Apr-2013 13:54:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editvideo_OpeningFcn, ...
                   'gui_OutputFcn',  @editvideo_OutputFcn, ...
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


% --- Executes just before editvideo is made visible.
function editvideo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to editvideo (see VARARGIN)

% Choose default command line output for editvideo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes editvideo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = editvideo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_view.
function btn_view_Callback(hObject, eventdata, handles)
% hObject    handle to btn_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mov = aviread(get(handles.edit_video,'string'),1);
figure;
imagesc(mov.cdata);
titre = get(handles.edit_video,'string')
title(titre);

function edit_video_Callback(hObject, eventdata, handles)
% hObject    handle to edit_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_video as text
%        str2double(get(hObject,'String')) returns contents of edit_video as a double


% --- Executes during object creation, after setting all properties.
function edit_video_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_browse.
function btn_browse_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path]= uigetfile('.avi')
set(handles.edit_video,'string',[path,file])



% --- Executes on button press in btn_crop.
function btn_crop_Callback(hObject, eventdata, handles)
% hObject    handle to btn_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dy = str2num(get(handles.edit_cropx,'string'));
dx = str2num(get(handles.edit_cropy,'string'));
fileinfo = aviinfo(get(handles.edit_video,'string'));
dframe =  fileinfo.NumFrames;
try 
[pathstr, name, ext] = fileparts(get(handles.edit_video,'string'));
filein = [pathstr,'\',name,ext];
fileout = [pathstr,'\','c',name,ext];
aviobj = avifile(fileout,'fps',1);
if fileinfo.Width > numel(dx)
    disp(['witdh : ',num2str(numel(dx))])
end
if fileinfo.Height > numel(dy)
    disp(['Height : ',num2str(numel(dy))])
end
for i=1:dframe
    mov = aviread(filein,i); %initial à couper  
    mov2.cdata = mov.cdata(dx,dy,:); 
    mov2.colormap = [];      
    aviobj = addframe(aviobj,mov2);  %ajouter la frame au film
end
aviobj = close(aviobj);
disp(['Done Crop file : ', fileout, ' created'] )
catch
    aviobj = close(aviobj);
end


function edit_cropx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cropx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cropx as text
%        str2double(get(hObject,'String')) returns contents of edit_cropx as a double


% --- Executes during object creation, after setting all properties.
function edit_cropx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cropx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cropy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cropy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cropy as text
%        str2double(get(hObject,'String')) returns contents of edit_cropy as a double


% --- Executes during object creation, after setting all properties.
function edit_cropy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cropy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_browse2.
function btn_browse2_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browse2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file2,path2] = uigetfile('.avi');
set(handles.edit_file2,'string',[path2,file2]);


function edit_file2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_file2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_file2 as text
%        str2double(get(hObject,'String')) returns contents of edit_file2 as a double


% --- Executes during object creation, after setting all properties.
function edit_file2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_file2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_file1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_file1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_file1 as text
%        str2double(get(hObject,'String')) returns contents of edit_file1 as a double


% --- Executes during object creation, after setting all properties.
function edit_file1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_file1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_browse1.
function btn_browse1_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browse1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file1,path1]= uigetfile('.avi')
set(handles.edit_file1,'string',[path1,file1])

% --- Executes on button press in btn_merge_Horizontal.
function btn_merge_Horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to btn_merge_Horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file1 = get(handles.edit_file1,'string')
[pathstr, name1, ext] = fileparts(file1)
file2 = get(handles.edit_file2,'string')
[pathstr2, name2, ext] = fileparts(file2)
fileout = [pathstr,'\',name1,'and',name2,ext]
fileinfo1 = aviinfo(file1)
fileinfo2 = aviinfo(file2)
badinfo = 0
if fileinfo1.NumFrames ~= fileinfo2.NumFrames
    disp(['File 1 Frames = ',num2str(fileinfo1.NumFrames), '           File 2 Frames',num2str(fileinfo2.NumFrames)])
    badinfo = 1;
end
if fileinfo1.Width ~=fileinfo2.Width
    disp(['File 1 width = ',num2str(fileinfo1.width), '              File 2 width',num2str(fileinfo2.width)])
    badinfo = 1
end
if badinfo
    return
end
try
aviobj = avifile(fileout,'fps',1);
for i=1:dframe
    mov = aviread(file1,i); %initial à couper  
    mov2 = aviread(file2,i); %
    movall.cdata =  cat(2,mov.cdata , mov2.cdata);
    movall.colormap = [];
    aviobj = addframe(aviobj,movall);  %ajouter la frame au film
end

catch
    aviobj = close(aviobj);
end

% --- Executes on button press in btn_MergeVertical.
function btn_MergeVertical_Callback(hObject, eventdata, handles)
% hObject    handle to btn_MergeVertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file1 = get(handles.edit_file1,'string')
[pathstr, name1, ext] = fileparts(file1)
file2 = get(handles.edit_file2,'string')
[pathstr2, name2, ext] = fileparts(file2)
fileout = [pathstr,'\',name1,'and',name2,ext]
fileinfo1 = aviinfo(file1)
fileinfo2 = aviinfo(file2)
badinfo = 0
if fileinfo1.NumFrames ~= fileinfo2.NumFrames
    disp(['File 1 Frames = ',num2str(fileinfo1.NumFrames), '           File 2 Frames',num2str(fileinfo2.NumFrames)])
    badinfo = 1;
end
if fileinfo1.Width ~=fileinfo2.Width
    disp(['File 1 width = ',num2str(fileinfo1.width), '              File 2 width',num2str(fileinfo2.width)])
    badinfo = 1
end
if badinfo
    return
end
try
aviobj = avifile(fileout,'fps',1);
for i=1:dframe
    mov = aviread(file1,i); %initial à couper  
    mov2 = aviread(file2,i); %
    movall.cdata =  cat(1,mov.cdata , mov2.cdata);
    movall.colormap = [];
    aviobj = addframe(aviobj,movall);  %ajouter la frame au film
end

catch
    aviobj = close(aviobj);
end


function edit_Time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Time as text
%        str2double(get(hObject,'String')) returns contents of edit_Time as a double


% --- Executes during object creation, after setting all properties.
function edit_Time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_labelidentification_Callback(hObject, eventdata, handles)
% hObject    handle to edit_labelidentification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_labelidentification as text
%        str2double(get(hObject,'String')) returns contents of edit_labelidentification as a double


% --- Executes during object creation, after setting all properties.
function edit_labelidentification_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_labelidentification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_title.
function btn_title_Callback(hObject, eventdata, handles)
% hObject    handle to btn_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileinfo = aviinfo(get(handles.edit_video,'string'));
[pathstr, name, ext] = fileparts(get(handles.edit_video,'string'));
filein = [pathstr,'\',name,ext];
fileout = [pathstr,'\','title',name,ext];
dframe =  fileinfo.NumFrames;
titlewidth = str2num(get(handles.edit_width_title,'string'));
xpos_time = str2num(get(handles.edit_xpos_time,'string'));

xpos_title = str2num(get(handles.edit_xpos_label,'string'));
ypos_time = str2num(get(handles.edit_yalign,'string'));
ypos_title = str2num(get(handles.edit_yalign,'string'));
aviobj = avifile(fileout,'fps',1) %fichier à creer
time = str2num(get(handles.edit_Time,'string'));
texttoadd = get(handles.edit_labelidentification,'string');
try
for i=1:dframe
    mov = aviread(filein,i); %initial à ajouter une bande blanche en haut avec un titre
    [dx,dy,dz]=size(mov.cdata);
    addupperband = uint8(ones(titlewidth,dy,dz)*255);
    mov2.cdata = cat(1,addupperband,mov.cdata);
    mov2.colormap = [];    
    %ADD time en haut
    h = figure    
    if numel(time)<i
        labelt=sprintf('%01.2f',time(end));
    else
        labelt=sprintf('%01.2f',time(i));
    end
    set(gca,'visible','off')
    text(0,0,['Time : ', labelt,'s'],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
    saveas(h,'temp','tif');
    close(h) 
    imgall = imread('temp','tif');     
    indy = find(sum(sum(imgall,3)~=765,1));
    indx = find(sum(sum(imgall,3)~=765,2));
    indx = indx(1):indx(end);
    indy = indy(1):indy(end);
    label_time = imgall(indx,indy,:); 
    %FIN DE ADD TIME (label_time)
    [dy,dx,c] = size(label_time);
    mov2.cdata(ypos_time:ypos_time+dy-1, xpos_time:xpos_time+dx-1,:) = label_time;   
    %%%
    % Add title
    %%%    
    h = figure    
    set(gca,'visible','off')
    text(0,0,[texttoadd],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
    saveas(h,'temp','tif');
    close(h)    
    imgall = imread('temp','tif');     
    indy = find(sum(sum(imgall,3)~=765,1));
    indx = find(sum(sum(imgall,3)~=765,2));
    indx = indx(1):indx(end);
    indy = indy(1):indy(end);
    titre = imgall(indx,indy,:); 
    [dy,dx,c] = size(titre);    
    mov2.cdata(ypos_title:ypos_title+dy-1, xpos_title:xpos_title+dx-1,:) =titre;
    if i==1
        figure
        imagesc(label_time)
        title(['Height = ',num2str(numel(indx)), ' Width = ', num2str(numel(indy))])
        figure
        imagesc(titre)
        title(['Height = ',num2str(numel(indx)), ' Width = ', num2str(numel(indy))])
        figure
        imagesc(mov2.cdata)    
    end
    
    aviobj = addframe(aviobj,mov2);  %ajouter la frame au film
end
aviobj = close(aviobj);
catch
    'error Rescale'
    aviobj = close(aviobj);
end



function edit_width_title_Callback(hObject, eventdata, handles)
% hObject    handle to edit_width_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_width_title as text
%        str2double(get(hObject,'String')) returns contents of edit_width_title as a double


% --- Executes during object creation, after setting all properties.
function edit_width_title_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_width_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xpos_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xpos_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xpos_time as text
%        str2double(get(hObject,'String')) returns contents of edit_xpos_time as a double


% --- Executes during object creation, after setting all properties.
function edit_xpos_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xpos_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xpos_label_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xpos_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xpos_label as text
%        str2double(get(hObject,'String')) returns contents of edit_xpos_label as a double


% --- Executes during object creation, after setting all properties.
function edit_xpos_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xpos_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_inserttitlename.
function btn_inserttitlename_Callback(hObject, eventdata, handles)
% hObject    handle to btn_inserttitlename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileinfo = aviinfo(get(handles.edit_video,'string'));
[pathstr, name, ext] = fileparts(get(handles.edit_video,'string'));
filein = [pathstr,'\',name,ext];
dframe =  1;
texttoadd=get(handles.edit_labelidentification,'string')
posx = 100;
posy = 20;
h = figure    
set(gca,'visible','off')
text(0,0,[texttoadd],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
titre = imgall(indx,indy,:); 
figure
imagesc(titre)
title(['Height = ',num2str(numel(indx)), ' Width = ', num2str(numel(indy))])
[dy,dx,c] = size(titre);
for i =dframe
    mov = aviread(filein,i); %initial 
    mov.colormap = [];        
    mov.cdata(posy:posy+dy-1, posx:posx+dx-1,:) =titre;
    figure
    imagesc(mov.cdata)
end




function edit_font_Callback(hObject, eventdata, handles)
% hObject    handle to edit_font (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_font as text
%        str2double(get(hObject,'String')) returns contents of edit_font as a double


% --- Executes during object creation, after setting all properties.
function edit_font_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_font (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yalign_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yalign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yalign as text
%        str2double(get(hObject,'String')) returns contents of edit_yalign as a double


% --- Executes during object creation, after setting all properties.
function edit_yalign_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yalign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 [file,path]= uigetfile('.tif')
set(handles.edit_legendfile,'string',[path,file])


function edit_legendfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_legendfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_legendfile as text
%        str2double(get(hObject,'String')) returns contents of edit_legendfile as a double


% --- Executes during object creation, after setting all properties.
function edit_legendfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_legendfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_addlegend.
function btn_addlegend_Callback(hObject, eventdata, handles)
% hObject    handle to btn_addlegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileinfo = aviinfo(get(handles.edit_video,'string'));
[pathstr, name, ext] = fileparts(get(handles.edit_video,'string'));
filein = [pathstr,'\',name,ext];
fileout = [pathstr,'\','title',name,ext];
dframe =  fileinfo.NumFrames;
filelegende = get(handles.edit_legendfile,'string');
imgall = imread(filelegende,'tif');     
[dy,dx,c] = size(imgall ); % attention il faut que le .tif soit plus petit que le film
aviobj = avifile(fileout,'fps',1)
try
for i =1:dframe
    mov = aviread(filein,i); %initial  
    [y,x,c]=size(mov.cdata);
    img=imageresize(imgall,y-100,[]); 
    [dy,dx,c] = size(img);
    posy = (y - dy)/2;
    posx = x;
    mov.cdata = cat(2,mov.cdata, ones(y,dx,3)*255);
    posy = (y - dy)/2;
    posx = x;
    mov.cdata(posy+1:posy+dy, posx+1:posx+dx,:) =img;
    aviobj = addframe(aviobj,mov);  %ajouter la frame au film
end
    aviobj = close(aviobj);
catch    
    'error ColorBar'
    aviobj = close(aviobj);
end



% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


