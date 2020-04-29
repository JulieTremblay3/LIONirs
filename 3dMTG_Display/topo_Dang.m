function topo_Dang(find_time,angle_cote,angle_face,dconc,h)
guiHOMER = getappdata(0,'guiHOMER');
guiHelmet = getappdata(0,'guiHelmet');
handles = guihandles(guiHelmet);
hsubaxes = gca;
cla
if isempty(guiHelmet)
    msgbox('Sorry, you should open the helmet .prj to enable this option')
    return
end
iscolorbar = 1;
istitle = 0;
pixelx = 1000;
pixely = 1000;
xmax = 0.13;
ymax = 0.13;
zmax = 0.13;
size_screen=get(0,'ScreenSize');
%set(h,'unit','pixel','position',[30,30,pixelx,pixely])
%Orientation of the topo
dist = 0.3;
set(gca,'visible','off')
cp_x = cos(angle_cote*pi/180) * cos(angle_face*pi/180)*dist;
cp_y = sin(angle_cote*pi/180)* cos(angle_face*pi/180)*dist;
cp_z = sin(angle_face*pi/180)*dist;
set(hsubaxes,'unit','normalize')
set(hsubaxes,'CameraTarget',[0,0,0]);
set(hsubaxes,'CameraPosition',[cp_x,cp_y,cp_z]);
set(hsubaxes,'CameraUpVector',[0,0,1]);
set(hsubaxes,'UserData','Movies');
set(hsubaxes,'xtick',[],'ytick',[],'ztick',[])
set(hsubaxes,'Position',[0,0,1,1]);
set(hsubaxes,'OuterPosition',[0,0,1,1])
set(hsubaxes,'Xlim',[-xmax,xmax]);
set(hsubaxes,'Ylim',[-ymax,ymax]);
set(hsubaxes,'Zlim',[-zmax/3, zmax]);
set(hsubaxes,'xtick',[])
set(hsubaxes,'ytick',[])
set(hsubaxes,'ztick',[])
Colormapfig = get(handles.IO_HelmetMTG,'colormap');
set(h,'colormap',Colormapfig);
prop = get(handles.axes_Mtg2,'CLim');
set(hsubaxes,'CLim',prop);
clear prop

PMI = get(guiHOMER,'UserData');
global currentsub
cf = PMI{currentsub}.currentFile;
echantillon_time = find(find_time > PMI{currentsub}.data(cf).HRF.tHRF);
echantillon_time = echantillon_time(end);
PrjStruct = getappdata(handles.IO_HelmetMTG,'PrjStruct');
if get(handles.radio_viewskin,'value')
    type=0;
else
    type=1;
end
typeactivation=get(handles.popupmenu_display,'value');
d1=selectimagetype(echantillon_time,typeactivation,dconc);
indna = isnan(d1);
d1(indna)=0;
PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type);
setappdata(handles.IO_HelmetMTG,'PrjStruct',PrjStruct);
resetview_axes(hsubaxes)
p = mfilename('fullpath'); 
[pathstr, name, ext]=fileparts(p);
saveas(h,[pathstr,filesep,'temp.tif'],'tif');

