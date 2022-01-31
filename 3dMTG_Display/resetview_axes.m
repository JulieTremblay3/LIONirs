function resetview_axes(haxes)    
    guiHelmet = getappdata(0,'guiHelmet');
    PrjStruct = getappdata(guiHelmet,'PrjStruct');  
    handles = guihandles(guiHelmet);
    DispHelm = get_Helmet(PrjStruct);
    sMtg = get_Mtg(DispHelm);
    vHoles = get_vHoles(DispHelm);
    try
    if get(handles.radio_guiSPMnirsHSJ,'value')==1
       guiHOMER = getappdata(0,'guiHOMER');
    elseif get(handles.radio_guiSPMnirsHSJ,'value')==2
       guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
    elseif get(handles.radio_guiSPMnirsHSJ,'value')==3
       guiHOMER = getappdata(0,'gui_SPMvideo');
    end
    catch
        1
    end
    PMI =get(guiHOMER,'UserData');      
    DispParameter.axes1 = haxes;
    DispParameter.viewhelmet = get(handles.radio_holes,'value');
    DispParameter.reset = 1;  
    DispParameter.channel = get(handles.radio_Channel,'value');
    DispParameter.viewfront = get(handles.radio_viewfront,'value');
    DispParameter.dist = get(handles.radio_viewdistance,'value');
    DispParameter.fiducies = get(handles.radio_viewfiducies,'value');
    DispParameter.viewskin = get(handles.radio_viewskin,'value');
    DispParameter.viewcortex = get(handles.radio_viewcortex,'value');
    DispParameter.viewatlas = get(handles.radio_viewatlas,'value');
    DispParameter.viewcortexintersection = get(handles.radio_normal_cortex,'value');
    DispParameter.cthresh = str2num(get(handles.edit_cminthreshold,'string'));
    DispParameter.cmin = str2num(get(handles.edit_climmin,'string'));
    DispParameter.cmax = str2num(get(handles.edit_climmax,'string'));
    DispParameter.MRIviewlight = get(handles.radio_light,'value');
    DispParameter.MRIviewtransparent = get(handles.radio_transparent,'value');
    DispParameter.viewMeasListAct = get(handles.radio_ViewMeasListAct,'value');
    DispParameter.enumaltas = str2num(get(handles.edit_atlasview,'string'));
    DispParameter.viewcortexandatlas  = get(handles.radio_cortexandatlas,'value'); 
    DispParameter.HideNotCover = get(handles.radio_show_cover,'value');
    DispParameter.LineWidth = 4;
    DispParameter.D1label = 0;
    IO_displayHelmet(PrjStruct,PMI,DispParameter);
    guidata(handles.IO_HelmetMTG, handles); 