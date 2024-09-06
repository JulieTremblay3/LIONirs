%Multiple display
function Multidisplay_IOhelmet(source,eventdata,h_figure,mode)
%Produce a video with the topographic map
%Do you want to save video and vcolor file
[file,pathavi]=uiputfile('.avi')
if file~=0
    mov = avifile([pathavi,file],'fps',1); 
end    
 
%Prepare the display
h_figure = gcf;
handlesfigure =  guihandles(h_figure);
% if get(handlesfigure.radio_avgvcolor,'value') %choose the path of the avg file to open
%         pathavg=uigetdir();
% end
NameFile = get(handlesfigure.edit_NameFile,'string')
angle = get(handlesfigure.editviewangle,'string');
numwindows  = str2num(get(handlesfigure.editwindows,'string'));
a_cote = [];
a_face = [];
rem=angle;
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
Haxes = findobj('userdata','Movies');
if ~isempty(Haxes)
    for i = 1:numel(Haxes)
        henlever = Haxes(i);
        set(henlever ,'visible','off');
        ch = get(henlever ,'children');
        delete(ch);
        delete(henlever);    
    end
end


%Use the data directly from Homer
guiHOMER = getappdata(0,'guiHOMER');
guiHelmet = getappdata(0,'guiHelmet');
if isempty(guiHelmet)
    msgbox('Sorry, you should open the helmet .prj to enable this option')
    return
end
start = str2num(get(handlesfigure.edittimestart,'string'));
postTime= str2num(get(handlesfigure.edittimestop,'string')) ;
if  postTime < 0 
    postTime = 0;
end

%% GRAPH 1
deltay = str2num(get(handlesfigure.editspace,'string'));       %ESPACEMENT Horizontal plus grand = plus rapporché
xaxiszoom = str2num(get(handlesfigure.edit_xaxiszoom,'string'));   %zoom echel de l'axe de chaque graphe
yaxiszoom = str2num(get(handlesfigure.edit_yaxiszoom,'string'));   %zoom echel de l'axe de chaque graphe
zaxiszoom = str2num(get(handlesfigure.edit_zaxiszoom,'string'));   %zoom echel de l'axe de chaque graphe
step_time = str2num(get(handlesfigure.edit_steptime,'string'));
speed = str2num(get(handlesfigure.editspeed,'string'));
xzoom = 1.2;
haut = 0.45;
large = 1/(numwindows);
ydistlabel = 0.2;
col = 1/numwindows - deltay;
ini_col = 0.05;
dist = 0.3;
haxes = [];
ini_row = 0.05;         % image du bas
row = 0.45;
hfiguredisplay=figure
set(hfiguredisplay,'unit', 'normalized','position',[0,0,1,1],'color',[1,1,1])
hback = axes
set(hback ,'units','normalize','position',[ini_col,ini_row,large*numwindows-deltay,haut*2])
set(hback,'xtick',[],'ytick',[],'ztick',[])
resetview_axes(hback);
xlim = get(hback,'clim')
delete(hback)
hback = axes
set(hback ,'units','normalize','position',[ini_col,ini_row,0.8,haut*2]);
set(hback,'xtick',[],'ytick',[],'ztick',[], 'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
caxis(xlim)
colorbar

%Orientation of the topo
for i=0:numwindows-1    
    hsubaxes = axes;%HbO
    haxes=[haxes,hsubaxes];
    angle_cote = a_cote(i+1);
    angle_face = a_face(i+1);
    cp_x = cos(angle_cote*pi/180) * cos(angle_face*pi/180)*dist;
    cp_y = sin(angle_cote*pi/180)* cos(angle_face*pi/180)*dist;
    cp_z = sin(angle_face*pi/180)*dist;  
    set(hsubaxes,'CameraTarget',[0,0,0]);
    set(hsubaxes,'CameraPosition',[cp_x,cp_y,cp_z]);
    set(hsubaxes,'CameraUpVector',[0,0,1]);
    set(hsubaxes,'UserData','Movies');
    set(hsubaxes,'units','normalize','position',[ini_col+col*i,ini_row,large,haut]);    
    set(hsubaxes,'xtick',[],'ytick',[],'ztick',[])
    axis([-1*xaxiszoom,1*xaxiszoom, -1*yaxiszoom,1*yaxiszoom,-1*zaxiszoom,1*zaxiszoom])
    hsubaxes = axes;%HbR
    haxes=[haxes,hsubaxes]; 
    set(hsubaxes,'CameraTarget',[0,0,0]);
    set(hsubaxes,'CameraPosition',[cp_x,cp_y,cp_z]);
    set(hsubaxes,'CameraUpVector',[0,0,1]);
    set( hsubaxes,'units','normalize','position',[ini_col+col*i,row,large,haut]);
    set(hsubaxes,'UserData','Movies');
    axis([-1*xaxiszoom,1*xaxiszoom, -1*yaxiszoom,1*yaxiszoom,-1*zaxiszoom,1*zaxiszoom])
    set(hsubaxes,'xtick',[],'ytick',[],'ztick',[])
end


handles = guihandles(guiHelmet);
PMI = get(guiHOMER,'UserData'); 
global currentsub
cf = PMI{currentsub}.currentFile; 
            
 %Mapping start here
try 
    if get(handlesfigure.radio_avgmat,'value')|get(handlesfigure.radio_avgvcolor,'value') % ouverture des fichiers .mat pour la topo 
        setappdata(0,'UseNativeSystemDialogs',false)
        if get(handlesfigure.radio_avgmat,'value')
          [file, path]=  uigetfile('.mat','multiselect','on');
        elseif get(handlesfigure.radio_avgvcolor,'value')
          [file, path]=  uigetfile('.vcolor','multiselect','on');  
        end
        PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');   
        if get(handles.radio_viewskin,'value')
            type=0
        else
            type=1 
        end
        for ifile=1:numel(file)
            if get(handlesfigure.radio_avgmat,'value')
                typeactivation=1% Peut-importe
                load([path,file{ifile}],'-mat')
                d1 = A;
                hHbO=uicontrol('Style','text','string',[file{ifile}],'units','normalize',...
                    'Position',[0.01,0.9-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
                PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type);
                oMRI = get_MRI_Data(PrjStruct);
                vColor = get_SkinVcolor(oMRI);
                save([pathavi,file{ifile},'.vcolor'], 'vColor','-mat');    
                setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);         
            elseif get(handlesfigure.radio_avgvcolor,'value')
                oMRI = get_MRI_Data(PrjStruct);  
                hHbO=uicontrol('Style','text','string',[file{ifile}],'units','normalize',...
                    'Position',[0.01,0.9-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
                vColor=opentopo([path,file{ifile}]);
                if get(handles.radio_viewskin,'value')
                    [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
                    vColortmp = vColor(1:size(VertexBuffer,1));
                    oMRI = set_SkinVcolor(oMRI,vColortmp);
                    PrjStruct = set_MRI_Data(PrjStruct, oMRI);
                elseif get(handles.radio_viewcortex,'value')
                    [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI ); 
                    vColortmp = vColor(1:size(VertexBuffer,1));
                     oMRI = set_CortexLowResVcolor(oMRI,vColortmp);
                     PrjStruct = set_MRI_Data(PrjStruct, oMRI);
                else
                    disp('No display')
                    return
                end
                 setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct); 
                %HbO display
                for i = 2:2:numel(haxes)
                    resetview_axes(haxes(i));  
                end                
            end
    %HbO display
            for i = 2:2:numel(haxes)
                resetview_axes(haxes(i));  
            end 
            F = getframe(hfiguredisplay);
            mov = addframe(mov,F);
            uiwait(gcf,speed)
            delete(hHbO)
        end        
    end
    
if ~get(handlesfigure.radio_avgmat,'value')
    
    if ~isfield(PMI{currentsub}.data(cf).HRF,'tHRF')
     msgbox('Please average data before');
     return
end
fin = numel(PMI{currentsub}.data(cf).HRF.tHRF);      

htitle = uicontrol('Style','text','string',['T value'],'units','normalize',...
                    'Position',[0.01,0.9,0.05 0.06],'FontSize',12,'BackgroundColor',[1,1,1],'HorizontalAlignment','center');              
        
if  ~get(handlesfigure.radio_avgmat,'value')& ~get(handlesfigure.radio_avgvcolor,'value')
for find_time = start:step_time:(start+postTime)
     x = find_time  
     echantillon_time = find(find_time > PMI{currentsub}.data(cf).t);   
     echantillon_time = find(find_time > PMI{currentsub}.data(cf).HRF.tHRF);      
     echantillon_time = echantillon_time(end);
        %labeltime = sprintf('%4.2f',x)
     labeltime = fixdecimal2string(x,4,2)
     hHbO=uicontrol('Style','text','string',['HbO'],'units','normalize',...
                    'Position',[0.01,0.9-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
     hTime=uicontrol('Style','text','string',['Time   ', labeltime,'(s)'],'units','normalize',...
                    'Position',[0.01,0.85,0.1 0.1],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
                
     hHbR = uicontrol('Style','text','string',['HbR'],'units','normalize',...
                    'Position',[0.01,0.45-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');                
    if get(handles.radio_stop,'value')==1
        return
    end
    PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');   
    if get(handles.radio_viewskin,'value')
        type=0
    else
        type=1
    end

    % POUR HBO    
    if ~get(handlesfigure.radio_avgvcolor,'value')%use local hmr value
        typeactivation=get(handles.popupmenu_display,'value');
        d1=selectimagetype(echantillon_time,typeactivation,1);      
        indna = isnan(d1);
        d1(indna)=0;    
        PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type);
            if file~=0
                oMRI = get_MRI_Data(PrjStruct);
                vColor = get_SkinVcolor(oMRI);          
                 nametopo = [pathavi,file, 'HbOTime   ', labeltime,'(s)','.vcolor'];
                 savetopo(nametopo,vColor,1);
            end
        setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct); 
        for i = 2:2:numel(haxes)
            resetview_axes(haxes(i));  
        end  
        %POUR HBR    
        d1=selectimagetype(echantillon_time,typeactivation,2);
        indna = isnan(d1);
        d1(indna)=0; 
        PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type);   
        oMRI = get_MRI_Data(PrjStruct);
        vColor = get_SkinVcolor(oMRI);
        nametopo = [pathavi, file,'HbRTime   ', labeltime,'(s)','.vcolor'];
        savetopo(nametopo,vColor,1);
        setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct); 
        for i = 1:2:numel(haxes)
            resetview_axes(haxes(i));  
        end    
        if isfield(handlesfigure,'linetime')
            set(handlesfigure.linetime,'XData',[x,x],'YData',[-10e-6,10e-6]);
        else
            h=plot([x,x],[-10e-6,10e-6]);        
        end   
        %HbT
        d1=selectimagetype(echantillon_time,typeactivation,3);
        indna = isnan(d1);
        d1(indna)=0; 
        PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type);   
        oMRI = get_MRI_Data(PrjStruct);
        vColor = get_SkinVcolor(oMRI);
        nametopo = [pathavi, file,'HbTTime   ', labeltime,'(s)','.vcolor'];
        savetopo(nametopo,vColor,1);
      
        
    elseif 0 %use special vcolorlist avec les temps non fonctionnel actuellement 
        oMRI = get_MRI_Data(PrjStruct);
        %chaque seconde un fichier   
             hHbO=uicontrol('Style','text','string',[NameFile],'units','normalize',...
                    'Position',[0.01,0.9-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
            load([pathavg, ,filesep,NameFile, labeltime,'(s).vColor'],'-mat');
        if 1
            oMRI = set_SkinVcolor(oMRI,vColor);
            PrjStruct = set_MRI_Data(PrjStruct, oMRI);
        else
             oMRI = set_CortexLowResVcolor(oMRI,vColor);
        	 PrjStruct = set_MRI_Data(PrjStruct, oMRI);
        end
         setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct); 
        %HbO display
        for i = 2:2:numel(haxes)
            resetview_axes(haxes(i));  
        end
%         %HbR display        
%         load([pathavg, '\HbRAvg', num2str(x),'(s).vColor'],'-mat');
%         oMRI = set_SkinVcolor(oMRI,vColor);
%         PrjStruct = set_MRI_Data(PrjStruct, oMRI);
%         setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);     
%           for i = 1:2:numel(haxes)
%             resetview_axes(haxes(i));  
%           end        
    end
    
      if file~=0
        F = getframe(hfiguredisplay);
        mov = addframe(mov,F);
      end       
    uiwait(gcf,speed)
    delete(hHbO)
    delete(hHbR)
    delete(hTime)
end
end
delete(htitle)
delete(hback)
end
    if file~=0
        mov = close(mov);
    end
catch
   disp('ERROR Multidisplay_IOhelmet')
   mov = close(mov);
end