function update_ERLCOHView(GUI_LookMatrice_handles)
if get(GUI_LookMatrice_handles.radiobutton1,'value')
        GUI_LookMatrice_handles = plot_axes_MATRICE(GUI_LookMatrice_handles);
        GUI_LookMatrice_handles = plot_axes_connectogramEEG(GUI_LookMatrice_handles);
        GUI_LookMatrice_handles = plot_axes_timefreq(GUI_LookMatrice_handles);
end
