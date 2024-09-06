function varargout = SelectTimeTopo(varargin)
% SELECTTIMETOPO M-file for SelectTimeTopo.fig
%      SELECTTIMETOPO, by itself, creates a new SELECTTIMETOPO or raises the existing
%      singleton*.
%
%      H = SELECTTIMETOPO returns the handle to a new SELECTTIMETOPO or the handle to
%      the existing singleton*.
%
%      SELECTTIMETOPO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTTIMETOPO.M with the given input arguments.
%
%      SELECTTIMETOPO('Property','Value',...) creates a new SELECTTIMETOPO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectTimeTopo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectTimeTopo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectTimeTopo

% Last Modified by GUIDE v2.5 11-Jul-2011 14:18:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectTimeTopo_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectTimeTopo_OutputFcn, ...
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


% --- Executes just before SelectTimeTopo is made visible.
function SelectTimeTopo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectTimeTopo (see VARARGIN)

% Choose default command line output for SelectTimeTopo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectTimeTopo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SelectTimeTopo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time as text
%        str2double(get(hObject,'String')) returns contents of edit_time as a double


% --- Executes during object creation, after setting all properties.
function edit_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_angleface_Callback(hObject, eventdata, handles)
% hObject    handle to edit_angleface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_angleface as text
%        str2double(get(hObject,'String')) returns contents of edit_angleface as a double


% --- Executes during object creation, after setting all properties.
function edit_angleface_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_angleface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_anglecote_Callback(hObject, eventdata, handles)
% hObject    handle to edit_anglecote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_anglecote as text
%        str2double(get(hObject,'String')) returns contents of edit_anglecote as a double


% --- Executes during object creation, after setting all properties.
function edit_anglecote_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_anglecote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_name as text
%        str2double(get(hObject,'String')) returns contents of edit_name as a double


% --- Executes during object creation, after setting all properties.
function edit_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_conc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_conc as text
%        str2double(get(hObject,'String')) returns contents of edit_conc as a double


% --- Executes during object creation, after setting all properties.
function edit_conc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_conc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_display.
function btn_display_Callback(hObject, eventdata, handles)
% hObject    handle to btn_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
time = str2num(get(handles.edit_time,'string'));

outputname = 'epi127_Sz15_Bloc1_610_cTopo2'
angle_cote = [70,-70,70,-70];
%angle_cote = [120,-120,120,-120]; %103
angle_cote = str2num(get(handles.edit_angleface,'string')); %137
angle_face = str2num(get(handles.edit_anglecote,'string'));
conc = str2num(get(handles.edit_conc,'string'));
% plot([time],-3*ones(numel(time),1),'marker','^','markersize',6,'color','k','linewidth',2,'linestyle','none','markerfacecolor','k')
% % % %Ligne EEG crise
% plot([0 4],[0,0],'linewidth',4,'color','k')
%  mark = [3.731541
% 7.242654
% 8.262975
% 9.563568
% 10.564997
% 11.662373
% 12.708943
% 13.958732
% 15.542135
% 16.852889]
%  plot([mark],zeros(numel(mark),1),'linewidth',4,'linestyle','none','marker','+','markersize',12,'color','k')


TotalV = [];
tic
indytext = 0;
p = mfilename('fullpath');
hfigure = figure;
hfigurestd = figure;
[pathstr, name, ext, versn]=fileparts(p);
for iplan = 1:numel(angle_cote)
    TotalH = [];
    for itime=1:numel(time) 
        figure(hfigure)    
        topo_Dang(time(itime),angle_cote(iplan),angle_face(iplan),conc(iplan),hfigure);
        imgall = imread([pathstr,'\temp.tif'],'tif');   
        indy = find(sum(sum(imgall,3)~=765,1));
        indx = find(sum(sum(imgall,3)~=765,2));
%         if iplan==1|iplan==3 %gauche epi114
%             indx = (indx(1)-20):(indx(end)-60);
%             indy = (indy(1)-20):(indy(end)-40); 
%         else %droite
%             indx = (indx(1)-20):(indx(end)-60);
%             indy = (indy(1)+40):(indy(end)+20);
%         end
        if iplan==1|iplan==3 %gauche epi114
            indx = (indx(1)-20):(indx(end)-10);
            indy = (indy(1)-20):(indy(end)+20); 
        else %droite
            indx = (indx(1)-20):(indx(end)-10);
            indy = (indy(1)-20):(indy(end)+20);
        end
        clear img
        img = imgall(indx,indy,:);     
        if iplan == 1 %Temps marqués sur la première colonne     
            figure(hfigurestd)
            cla
            set(gca,'visible','off')
            text(0,0,[num2str(time(itime)),' s'],'FontSize',30,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','arial')
            saveas(hfigurestd,[pathstr,'\temp.tif'],'tif');                   
            clear imgall titre
            imgall = imread([pathstr,'\temp.tif'],'tif'); 
            if indytext == 0 %S'Assurer d'avoir toujours le même découpage en hauteur en initialisant avec le premier
                indytext = find(sum(sum(imgall,3)~=765,2));              
            end
            indx = find(sum(sum(imgall,3)~=765,1));            
            titre = imgall(indytext,indx,:);            
            [dy,dx,c] = size(titre);
            posx = 50;
            posy = 1;
            Margeup = ones(dy,size(img,2),3)*255;
            Margeup(posy:posy+dy-1, posx:posx+dx-1,:)=titre;       
            img = [Margeup;img];
        end
          TotalH = [TotalH,img]; %concatenation horizontal             
    end  
    if size(TotalH,2)< size(TotalV,2)
        addy = size(TotalV,2)-size(TotalH,2);
        TotalH = [TotalH,ones(size(TotalH,1),addy,3)*255];
    elseif size(TotalH,2)> size(TotalV,2)& ~isempty(TotalV)
        TotalH(:,size(TotalV,2)+1:end,:)=[];
    end
    TotalV = [TotalV;TotalH];
end
%HbO

figure(hfigurestd)
cla
posx = 1;
posy = 1;
set(gca,'visible','off')
text(0.5,0.5,['HbO'],'FontSize',50,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(hfigurestd,[pathstr,'\temp.tif'],'tif');   
imgall = imread([pathstr,'\temp.tif'],'tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = (indx(1)-20):(indx(end)+20);
indy = (indy(1)-20):(indy(end)+20);
titre = imgall(indx,indy,:); 
[dy,dx,c] = size(titre);
%ajouter une marge
MargeLeft = ones(size(TotalV,1),dx,3)*255;
MargeLeft(posy:posy+dy-1, posx:posx+dx-1,:)=titre;

%HbR
cla
posx = 1;
[y,x,c]=size(TotalV);
posy = y/2;
set(gca,'visible','off')
text(0.5,0.5,['HbR'],'FontSize',50,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(hfigurestd,[pathstr,'\temp.tif'],'tif'); 
imgall = imread([pathstr,'\temp.tif'],'tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = (indx(1)-20):(indx(end)+20);
indy = (indy(1)-20):(indy(end)+20);
titre = imgall(indx,indy,:); 
[dy,dx,c] = size(titre);
MargeLeft(posy:posy+dy-1, posx:posx+dx-1,:)=titre;
TotalV = [MargeLeft,TotalV];

%color bar
cla
guiHelmet = getappdata(0,'guiHelmet');
handles = guihandles(guiHelmet);
set(gca,'visible','off')
Colormapfig = get(handles.IO_HelmetMTG,'colormap');
set(hfigurestd,'colormap',Colormapfig);
prop = get(handles.axes_Mtg2,'CLim');
set(gca,'CLim',prop);
clear prop
hcolorbar = colorbar;
set(hcolorbar,'fontsize',20)
saveas(hfigurestd,[pathstr,'\temp.tif'],'tif');
imgall = imread([pathstr,'\temp.tif'],'tif');
imgall = imresize(imgall,1.5);
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = ((indx(1)):(indx(end)));
indy = ((indy(1)):(indy(end)));
titre = imgall(indx,indy,:); 
titre = [titre;ones(30,size(titre,2),3)*255];
titre = [ones(30,size(titre,2),3)*255;titre];
titre = [titre,ones(size(titre,1),30,3)*255];
titre = [ones(size(titre,1),30,3)*255,titre];
[dy,dx,c] = size(titre);
[y,x,c]=size(TotalV);
posx = 1;
posy = (y - dy)/2;
MargeRight = ones(y,dx,3)*255;
%MargeRight(posy+1:posy+dy,posx+1:posx+dx,:)=titre;
TotalV = [TotalV,MargeRight];
toc
close(hfigurestd)
close(hfigure)
t1 = sprintf('%02.0f',time(1))
t2 = sprintf('%02.0f',time(end))
[nameout,pathout]=uiputfile([outputname,'t',t1,'to',t2,'.tif']);
imwrite(TotalV,[pathout,nameout],'tif');