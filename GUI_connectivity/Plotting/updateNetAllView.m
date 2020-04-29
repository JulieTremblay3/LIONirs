function updateNetAllView(GUI_LookMatrice_handles)
 
        %GUI_LookMat = getappdata(0,'GUI_LookMat');
       % GUI_LookMatrice_handles = guihandles(GUI_LookMat); 
        %DATA = get(GUI_LookMatrice_handles.GUI_LookMat,'UserData');
        GUI_LookMatrice_handles = plot_axesAij(GUI_LookMatrice_handles);
        GUI_LookMatrice_handles = plot_histogrammeAij(GUI_LookMatrice_handles);
        if get(GUI_LookMatrice_handles.popupmenu_linkoption,'value')==1 
            if get(GUI_LookMatrice_handles.popupmenu_view,'value')==1 %channel mode
            GUI_LookMatrice_handles = plot_axes_linkConnectogramm(GUI_LookMatrice_handles);
            end   
        elseif get(GUI_LookMatrice_handles.popupmenu_linkoption,'value')==2
            1
            GUI_LookMatrice_handles = plot_axes_linkMAP(GUI_LookMatrice_handles);
        elseif    get(GUI_LookMatrice_handles.popupmenu_linkoption,'value')==3
            
        end
        
% %          plot_axes_linkConnectogramm
%         guiHelmet = getappdata(0,'GUI_ViewNetwork');       
%         if ~isempty(guiHelmet)
%              guiHelmet_handles = guihandles( guiHelmet)
%             plot_axesviewNET(GUI_LookMatrice_handles,guiHelmet_handles);
%         end
%         guiHelmet = getappdata(0,'guiHelmet');
%         if ~isempty(guiHelmet)
%             Helmethandles = guihandles(guiHelmet); 
%             fhresetview = getappdata(guiHelmet,'fhresetview');
%             fhresetview(Helmethandles);
%         end