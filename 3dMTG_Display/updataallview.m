function updataallview()
        global currentsub;
        guiHOMER = getappdata(0,'guiHOMER');
        HOMERhandles = guihandles(guiHOMER); 
        PMI = get(guiHOMER,'UserData'); 
        cf = PMI{currentsub}.currentFile; 
        if get(HOMERhandles.HOMer_Tab,'value')
            plotAxes_SDG( HOMERhandles, 0 );
            plotAxes3(HOMERhandles, 0 );
            plotAxes1( HOMERhandles, 0, 0, 0 );
        end
        if get(HOMERhandles.LISA_tab,'value')
            plotAxes_LisaA(HOMERhandles,0);
            plotAxes_LisaB(HOMERhandles,0);          
            plotAxes_LisaAmoinsB(HOMERhandles,0);
        end
        if get(HOMERhandles.MMARGE_Tab,'value')
            plot_MMARGE_axes_avg(HOMERhandles,0);
            plot_MMARGE_axes_avgHistogramme(HOMERhandles,0);
            plot_MMARGE_axes_zonepeak(HOMERhandles,0);
        end
        if ~isempty(getappdata(0,'guiGML'))
            guiGML = getappdata(0,'guiGML');
            GMLhandles = guihandles(guiGML); 
            plot_axes_donnees(GMLhandles,0);
            plot_axes_corr(GMLhandles,0);
            plot_axes_yfrh(GMLhandles,0);
        end
        guiHelmet = getappdata(0,'guiHelmet');
        if ~isempty(guiHelmet)
            Helmethandles = guihandles(guiHelmet); 
            fhresetview = getappdata(guiHelmet,'fhresetview');
            fhresetview(Helmethandles);
        end