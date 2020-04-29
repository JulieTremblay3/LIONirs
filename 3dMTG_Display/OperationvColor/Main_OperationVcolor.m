function varargout = Main_OperationVcolor(varargin)
% MAIN_OPERATIONVCOLOR M-file for Main_OperationVcolor.fig
%      MAIN_OPERATIONVCOLOR, by itself, creates a new MAIN_OPERATIONVCOLOR or raises the existing
%      singleton*.
%
%      H = MAIN_OPERATIONVCOLOR returns the handle to a new MAIN_OPERATIONVCOLOR or the handle to
%      the existing singleton*.
%
%      MAIN_OPERATIONVCOLOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_OPERATIONVCOLOR.M with the given input arguments.
%
%      MAIN_OPERATIONVCOLOR('Property','Value',...) creates a new MAIN_OPERATIONVCOLOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_OperationVcolor_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_OperationVcolor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main_OperationVcolor

% Last Modified by GUIDE v2.5 04-Jul-2012 15:32:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Main_OperationVcolor_OpeningFcn, ...
    'gui_OutputFcn',  @Main_OperationVcolor_OutputFcn, ...
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


% --- Executes just before Main_OperationVcolor is made visible.
function Main_OperationVcolor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main_OperationVcolor (see VARARGIN)

% Choose default command line output for Main_OperationVcolor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Main_OperationVcolor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Main_OperationVcolor_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_browse.
function btn_browse_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'UseNativeSystemDialogs',false)
[name,path]=uigetfile('.vcolor','multiselect','on')
if path~=0
    if ~iscell(name)
        name = {name}
    end
    set(handles.edit_Path,'string',path)
    set(handles.list_name,'string',name)
end
% --- Executes on selection change in list_name.
function list_name_Callback(hObject, eventdata, handles)
% hObject    handle to list_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_name contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_name


% --- Executes during object creation, after setting all properties.
function list_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Path as text
%        str2double(get(hObject,'String')) returns contents of edit_Path as a double


% --- Executes during object creation, after setting all properties.
function edit_Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_MEAN.
function btn_MEAN_Callback(hObject, eventdata, handles)
% hObject    handle to btn_MEAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path=get(handles.edit_Path,'string');
name = get(handles.list_name,'string');
all = [];
for iname=1:numel(name)
    load([path,name{iname}],'-mat');
    all = [all,vColor];
end
vColor = mean(all')';
[nameout,pathout]= uiputfile([path,'Avg.vColor'])
if nameout==0
    return
end
save([pathout,nameout],'vColor','-mat')
msgbox([' Avg done, see file : ', pathout,nameout])




% --- Executes on button press in btn_BrowseSRX.
function btn_BrowseSRX_Callback(hObject, eventdata, handles)
% hObject    handle to btn_BrowseSRX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[name,path]=  uigetfile('.srx');
if name==0
    return
end
set(handles.edit_vertex,'string',[path,name])

function edit_vertex_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vertex as text
%        str2double(get(hObject,'String')) returns contents of edit_vertex as a double


% --- Executes during object creation, after setting all properties.
function edit_vertex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_LI.
function btn_LI_Callback(hObject, eventdata, handles,type)
% hObject    handle to btn_LI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%type = 1Laterality mask
%type = 2 Laterality mask
type = 2;
path=get(handles.edit_Path,'string');
name = get(handles.list_name,'string');
filevertex = get(handles.edit_vertex,'string');
[VertexBuffer faces_matrix, NV, NT] = Read_iMagic_SRX( filevertex );
if type == 1
    voxelheadcenter = str2num(get(handles.edit_voxelheadcenter,'string'));
    idLeft= find(VertexBuffer(:,1)>voxelheadcenter/2);
    idRight = find(VertexBuffer(:,1)<voxelheadcenter/2);
    idall = 1:size(VertexBuffer,1);
else
    load(get(handles.edit_MaskL,'string'),'-mat');
    idLeft = find(vColor>0);
    load(get(handles.edit_MaskR,'string'),'-mat');
    idRight = find(vColor >0);
    idall = [idLeft;idRight];
end
 
warning off
figure 
hold on
plot([-10,10],[0,0],'k')
for iname=1:numel(name)
    load([path,name{iname}],'-mat');     
    threshold = min(vColor(idall)):0.1:max(vColor(idall));
    max_tr = max(vColor(idall));
    min_tr = min(vColor(idall));
    
    vleft = vColor;
    vleft(idRight)=0;
    vright= vColor;
    vright(idLeft) =0;
    clear IL
    for i=1:numel(threshold)
        if threshold(i)> 0
            IL(i) = (numel(find(vleft>threshold(i)))-numel(find(vright>threshold(i))))/(numel(find(vleft>threshold(i)))+numel(find(vright>threshold(i))));
            thresholdRatio(i) = threshold(i)/max_tr;
        else
            IL(i) = (numel(find(vleft<threshold(i)))-numel(find(vright<threshold(i))))/(numel(find(vleft<threshold(i)))+numel(find(vright<threshold(i))));
            thresholdRatio(i) = threshold(i)/min_tr;
        end
    end
    
    
    h=plot(threshold,IL,'linewidth',3);
    sujet = name{iname};
    set(h,'displayname',sujet);
    title(sujet);
    set(gca,'fontsize',14);
    ylim([-1.2,1.2]);
    max(vColor)
    xlabel(['T value , max = ', num2str(max(vColor(idall))), '  min = ', num2str(min(vColor(idall)))]);
    ylabel(['LI ']);

    i = find( threshold<=threshold(end)/2);
    set(handles.edit_LImax,'string', num2str(IL(i(end))));
    set(handles.edit_Tmax,'string', num2str(threshold(i(end))));
    i = find( threshold<=threshold(1)/2);
    set(handles.edit_LImin,'string', num2str(IL(i(end))));
    set(handles.edit_Tmin,'string', num2str(threshold(i(end))));

    saveas(h,[path,sujet(1:end-6),'.jpg'],'jpg');
end
warning on



function edit_maskL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maskL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maskL as text
%        str2double(get(hObject,'String')) returns contents of edit_maskL as a double


% --- Executes during object creation, after setting all properties.
function edit_maskL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maskL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_MaskL.
function btn_MaskL_Callback(hObject, eventdata, handles)
% hObject    handle to btn_MaskL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name,path] = uigetfile('.vcolor')
if name==0
    return
end
set(handles.edit_MaskL,'string',[path,name])


function edit_maskR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maskR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maskR as text
%        str2double(get(hObject,'String')) returns contents of edit_maskR as a double


% --- Executes during object creation, after setting all properties.
function edit_maskR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maskR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_maskR.
function btn_maskR_Callback(hObject, eventdata, handles)
% hObject    handle to btn_maskR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[name,path] = uigetfile('.vcolor')
if name==0
    return
end
set(handles.edit_MaskR,'string',[path,name])

function edit_MaskL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaskL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaskL as text
%        str2double(get(hObject,'String')) returns contents of edit_MaskL as a double


% --- Executes during object creation, after setting all properties.
function edit_MaskL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaskL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MaskR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaskR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaskR as text
%        str2double(get(hObject,'String')) returns contents of edit_MaskR as a double


% --- Executes during object creation, after setting all properties.
function edit_MaskR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaskR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2



function edit_voxelheadcenter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_voxelheadcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_voxelheadcenter as text
%        str2double(get(hObject,'String')) returns contents of edit_voxelheadcenter as a double


% --- Executes during object creation, after setting all properties.
function edit_voxelheadcenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_voxelheadcenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_Mask.
function btn_Mask_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edit_tMAX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tMAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tMAX as text
%        str2double(get(hObject,'String')) returns contents of edit_tMAX as a double


% --- Executes during object creation, after setting all properties.
function edit_tMAX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tMAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LImax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LImax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LImax as text
%        str2double(get(hObject,'String')) returns contents of edit_LImax as a double


% --- Executes during object creation, after setting all properties.
function edit_LImax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LImax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Tmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Tmin as text
%        str2double(get(hObject,'String')) returns contents of edit_Tmin as a double


% --- Executes during object creation, after setting all properties.
function edit_Tmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LImin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LImin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LImin as text
%        str2double(get(hObject,'String')) returns contents of edit_LImin as a double


% --- Executes during object creation, after setting all properties.
function edit_LImin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LImin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Tmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Tmax as text
%        str2double(get(hObject,'String')) returns contents of edit_Tmax as a double


% --- Executes during object creation, after setting all properties.
function edit_Tmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_vertexsrx2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vertexsrx2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vertexsrx2 as text
%        str2double(get(hObject,'String')) returns contents of edit_vertexsrx2 as a double


% --- Executes during object creation, after setting all properties.
function edit_vertexsrx2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vertexsrx2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Browse_SRX2.
function btn_Browse_SRX2_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Browse_SRX2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[name,path]=  uigetfile('.srx');
if name==0
    return
end
set(handles.edit_vertexsrx2,'string',[path,name])

% --- Executes on button press in btn_apply.
function btn_apply_Callback(hObject, eventdata, handles)
% hObject    handle to btn_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filevertex = get(handles.edit_vertex,'string');
[VertexBuffer faces_matrix, NV, NT] = Read_iMagic_SRX( filevertex );
filevertex = get(handles.edit_vertexsrx2,'string');
[VertexBuffer2 faces_matrix2, NV2, NT2] = Read_iMagic_SRX( filevertex );
path=get(handles.edit_Path,'string');
name = get(handles.list_name,'string');
if get(handles.popup_projectionmode,'value')==1
    %Replacer au centre        
    VertexBuffertmp(:,1) = VertexBuffer(:,1)-max(VertexBuffer(:,1))/2;
    VertexBuffertmp(:,2) = VertexBuffer(:,2)-max(VertexBuffer(:,2))/2+30;
    VertexBuffertmp(:,3) = VertexBuffer(:,3)-max(VertexBuffer(:,3))/2;
    %Replacer au centre        
    VertexBuffer2(:,1) = VertexBuffer2(:,1)-max(VertexBuffer(:,1))/2; %gauche droite
    VertexBuffer2(:,2) = VertexBuffer2(:,2)-max(VertexBuffer(:,2))/2+30;
    VertexBuffer2(:,3) = VertexBuffer2(:,3)-max(VertexBuffer(:,3))/2; %avant arrière
    [Vertex_sph(:,1),Vertex_sph(:,2),Vertex_sph(:,3)] = cart2sph(VertexBuffertmp(:,1),VertexBuffertmp(:,2),VertexBuffertmp(:,3));
    [Vertex_sph2(:,1),Vertex_sph2(:,2),Vertex_sph2(:,3)] = cart2sph(VertexBuffer2(:,1),VertexBuffer2(:,2),VertexBuffer2(:,3));
    middle = max(VertexBuffer(:,1))/2;
    vColor2 = zeros(size(VertexBuffer2,1),1);
end

for iname = 1:numel(name)
    load([path,name{iname}],'-mat');
    ind = find(vColor);        
    vColor2(:)=0;
    if get(handles.popup_projectionmode,'value')==2 %direct projection % faster good on the side
        delta =5;
        for i = 1:numel(ind)
            if VertexBuffer(ind(i),1)>middle
                xmin = middle;
                xmax= VertexBuffer(ind(i),1)+delta;
            else
                xmin= VertexBuffer(ind(i),1)-delta;
                xmax = middle;
            end
            ymin = VertexBuffer(ind(i),2)-delta;
            ymax= VertexBuffer(ind(i),2)+delta;
            zmin = VertexBuffer(ind(i),3)-delta;
            zmax= VertexBuffer(ind(i),3)+delta;
            %     x = find(VertexBuffer2(:,1)> xmin & VertexBuffer2(:,1)< xmax);
            %     y = find(VertexBuffer2(:,2)> ymin & VertexBuffer2(:,2)< ymax);
            %     z = find(VertexBuffer2(:,3)> zmin & VertexBuffer2(:,3)< zmax);
            %     figure
            %     hold on
            %     plot(x,'x')
            %     plot(y,'rv')
            %       plot(y,'k+')
            indcolor = find(VertexBuffer2(:,1)> xmin & VertexBuffer2(:,1)< xmax&VertexBuffer2(:,2)> ymin & VertexBuffer2(:,2)< ymax & VertexBuffer2(:,3)> zmin & VertexBuffer2(:,3)< zmax) ;
            vColor2(indcolor) = vColor(ind(i));
        end
    elseif get(handles.popup_projectionmode,'value')==1
        clear VertexBuffertmp 
        clear VertexBuffer
        clear VertexBuffer2
        for i = 1:numel(ind)
            i/numel(ind) 
            theta = Vertex_sph(ind(i),1); %centre de l'activation topo
            phi = Vertex_sph(ind(i),2);                 
            indcolor =  find(theta+0.05 > Vertex_sph2(:,1)&...
                    theta-0.05 < Vertex_sph2(:,1)&...
                    phi+0.05 > Vertex_sph2(:,2)&...
                    phi-0.05 < Vertex_sph2(:,2));                 
            vColor2(indcolor) = vColor(ind(i));
        end  
    vColor = vColor2;
    save([path,'cortex',name{iname}],'vColor','-mat')
    max(vColor)
    elseif get(handles.popup_projectionmode,'value')==3 %circulat projection % faster good on the side
        theta = 5:10:45; %echelon de l'angle
        dx = sin(theta*pi/180)*20;
        dy = cos(theta*pi/180)*20;
        buffer =0;
        tic
        for i = 1:numel(ind)
            x = VertexBuffer(ind(i),1); %centre de l'activation topo
            y = VertexBuffer(ind(i),2);
            z = VertexBuffer(ind(i),3);
            if buffer == 30
                i/numel(ind)
                buffer = 0;
            end
            buffer = buffer+1;            
            for j=1:numel(dx) %axe des x
                indcolor =  find(x+dx(j) > VertexBuffer2(:,1)&...
                    x-dx(j) < VertexBuffer2(:,1)&...
                    y+dy(j) > VertexBuffer2(:,2)&...
                    y-dy(j) < VertexBuffer2(:,2)&...
                    z+dy(j) > VertexBuffer2(:,3)&...
                    z-dy(j) < VertexBuffer2(:,3));
                   vColor2(indcolor) = vColor(ind(i));
            end

            for j=1:numel(dx) %axe des y (small box x,z
                indcolor =  find(x+dy(j) > VertexBuffer2(:,1)&...
                    x-dy(j) < VertexBuffer2(:,1)&...
                    y+dx(j) > VertexBuffer2(:,2)&...
                    y-dx(j) < VertexBuffer2(:,2)&...
                    z+dy(j) > VertexBuffer2(:,3)&...
                    z-dy(j) < VertexBuffer2(:,3));
                    vColor2(indcolor) = vColor(ind(i));
            end
            for j=1:numel(dx) %axe des z (small box x,y)
                indcolor =  find(x+dy(j) > VertexBuffer2(:,1)&...
                    x-dy(j) < VertexBuffer2(:,1)&...
                    y+dy(j) > VertexBuffer2(:,2)&...
                    y-dy(j) < VertexBuffer2(:,2)&...
                    z+dx(j) > VertexBuffer2(:,3)&...
                    z-dx(j) < VertexBuffer2(:,3));
                    vColor2(indcolor) = vColor(ind(i));
            end
        end
        toc
    end
    vColor = vColor2;
    save([path,'cortex',name{iname}],'vColor','-mat')
    max(vColor)
end


% --- Executes on button press in btn_LItime.
function btn_LItime_Callback(hObject, eventdata, handles)
% hObject    handle to btn_LItime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
type = 0 %Laterality mask
%type = 2 Laterality mask
path=get(handles.edit_Path,'string');
name = get(handles.list_name,'string');


if 1==get(handles.radio_hemi,'value')
    filevertex = get(handles.edit_vertex,'string');
[VertexBuffer faces_matrix, NV, NT] = Read_iMagic_SRX( filevertex );
    type = 0
end
if 1==get(handles.radio_mask,'value')
    type=1
end
if type == 0
%     voxelheadcenter = str2num(get(handles.edit_voxelheadcenter,'string'));
%     idLeft= find(VertexBuffer(:,1)>181/2);
%     idRight = find(VertexBuffer(:,1)<181/2);
%     load('F:\Fluence\Review\Brodman_Skin\900fve04CAT-10.0dConcHbR.vcolor','-mat')
    minx = min(VertexBuffer(:,1));
    maxx = max(VertexBuffer(:,1));
    center = minx + (maxx-minx)/2;
    idLeft = find(VertexBuffer(:,1)>center);
    idRight = find(VertexBuffer(:,1)<center);    
else
    load(get(handles.edit_MaskL,'string'),'-mat');
    idLeft = find(vColor>0|vColor <0);
    idLeft = find(vColor>0|vColor <0);
    load(get(handles.edit_MaskR,'string'),'-mat');
    idRight = find(vColor >0|vColor <0);
end
Thresholdpositif = 1;
warning off
hfigure=figure
h1 = subplot(3,1,1);
hold on
h2 = subplot(3,1,2);
hold on
h3 = subplot(3,1,3);
hold on
plot([-10,10],[0,0],'k')
ILmaxall = [];
ILmaxleft = [];
ILmaxright = [];
ILall = [];
ILextleft = [];
ILextright = [];
tall = [];
thresholdall = [];
labelthreshold = get(handles.edit_threshold,'string')
for iname=1:numel(name)
    load([path,name{iname}],'-mat');
     vColormask = zeros(size(vColor));
     vColormask(idLeft) = vColor(idLeft);
     vColormask(idRight) = vColor(idRight);
    if str2num(labelthreshold) > 0
        threshold = abs(max(vColormask))*str2num(labelthreshold);
    else
        threshold = -abs(min(vColormask))*str2num(labelthreshold);
    end
    vleft = vColor;
  %old bad
%     vleft(idRight)=0;
%     vright(idLeft) =0;
   %vleft
   vleft=zeros(size(vColor));
   vleft(idLeft)= vColor(idLeft);
   
   vright =zeros(size(vColor));
   vright(idRight)= vColor(idRight);
  
    clear IL
    for i=1:numel(threshold)
        if threshold(i)>= 0
            IL(i) = (numel(find(vleft>threshold(i)))-numel(find(vright>threshold(i))))/(numel(find(vleft>threshold(i)))+numel(find(vright>threshold(i))));
            extLeft = numel(find(vleft>threshold(i)))
            extRight = numel(find(vright>threshold(i)))
        else
            IL(i) = (numel(find(vleft<threshold(i)))-numel(find(vright<threshold(i))))/(numel(find(vleft<threshold(i)))+numel(find(vright<threshold(i))));
            extLeft = NaN
            extRight = NaN
         end
    end
    x = name{iname};
    t = str2num(x(end-17:end-10));
    savetopo(['left',x],vleft,0)  
    savetopo(['right',x],vright,0) 
   
    tall = [tall,t];
    axes(h1)
    h=plot(t,IL,'linewidth',3,'marker','v');
    sujet = name{iname};
    set(h,'displayname',sujet);
    title(sujet);
    set(gca,'fontsize',14);
    ylim([-1.2,1.2]);
    max(vColor);
    ILall = [ILall;IL]; 
    ILextleft = [ILextleft ;extLeft];
    ILextright = [ILextright ;extRight];
    
    axes(h2)
    ILmax = (max(vleft)-max(vright))/(max(vleft)+max(vright))

    plot(t,ILmax,'marker','v','linewidth',3);
    ILmaxall = [ILmaxall;ILmax];
    if threshold(i)>= 0
    ILmaxleft = [ILmaxleft;max(vleft)];
    ILmaxright = [ILmaxright;max(vright)];
    else
        
    end
    
    axes(h3)
    thresholdall = [thresholdall,threshold];
    plot(t,threshold,'marker','v','linewidth',3);

    %         i = find( threshold<=threshold(end)/2);
    %         set(handles.edit_LImax,'string', num2str(IL(i(end))));
    %          set(handles.edit_Tmax,'string', num2str(threshold(i(end))));
    %          i = find( threshold<=threshold(1)/2);
    %         set(handles.edit_LImin,'string', num2str(IL(i(end))));
    %          set(handles.edit_Tmin,'string', num2str(threshold(i(end))));
end


axes(h1);
xlim([1 iname]);
plot([1 iname],[0,0]);
xlabel(['Time']);
ylabel(['LI ']);

axes(h2);
xlim([1 iname]);
ylim([-1.2,1.2]);
plot([1 iname],[0,0]);
xlabel(['Time']);
xlabel(['LI Max left vs max right']);


axes(h3);
xlim([1 iname]);
plot([1 iname],[0,0]);
xlabel(['Time']);
xlabel(['Treshold max/2']);
name = fliplr(strtok(fliplr(path),'\'))
if type==0
    path = [get(handles.edit_pathout,'string'),'Hemi\']
else
    path = [get(handles.edit_pathout,'string'),'Mask\']
end
saveas(hfigure,[path,'LI',name,'.aviHbOTime   ',labelthreshold,'.jpg'],'jpg');
save([path,'LI',name,'.aviHbOTime   ',labelthreshold,'.mat'],'tall','ILmaxall','ILmaxleft','ILmaxright','ILall','ILextleft','ILextright','thresholdall','-mat')
warning on



function edit_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_threshold as a double


% --- Executes during object creation, after setting all properties.
function edit_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popup_projectionmode.
function popup_projectionmode_Callback(hObject, eventdata, handles)
% hObject    handle to popup_projectionmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_projectionmode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_projectionmode


% --- Executes during object creation, after setting all properties.
function popup_projectionmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_projectionmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_pathout_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathout as text
%        str2double(get(hObject,'String')) returns contents of edit_pathout as a double


% --- Executes during object creation, after setting all properties.
function edit_pathout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


