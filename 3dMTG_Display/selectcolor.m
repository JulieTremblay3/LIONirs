
function selectcolor(hObject, eventdata, handles,i)
%Description : Permet de conserver la nouvelle couleur attribue au canal. 
newcolor = uisetcolor([0,0,0]);
set(hObject,'BackgroundColor',newcolor);
% if get(handles.radio_guiSPMnirsHSJ,'value')==1
%     guiHOMER = getappdata(0,'guiHOMER');
% elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
% elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
%     guiHOMER = getappdata(0,'gui_SPMvideo');
% end

PMI = get(guiHOMER,'UserData'); 
global currentsub;

try
    cf = PMI{currentsub}.currentFile;
catch
    cf=1
end
% 

PMI{currentsub}.color(i,:)= newcolor;
p = mfilename('fullpath');
color = PMI{currentsub}.color;
save(p,'color');
set(guiHOMER,'UserData',PMI);
%updataallview()