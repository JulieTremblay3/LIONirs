function handles = plotAxes_d1(varargin)
% handles.edit_noisemarker %MARKER SIZE
% handles.edit_stepvalue
%Data in PMI %guiHOMER = getappdata(0,'gui_SPMnirsHSJ') & PMI = get(guiHOMER,'UserData');
%handles.NIRS.Dt.fir.pp

parameter_linewidth = 2;
handles=varargin{1};
newfigure = 0;
if numel(varargin)>1
    newfigure = varargin{2};
end
if newfigure
    figure
else
%     figure(handles.figure1)
%     gca;
    axes(handles.display_axes);
end
cla

noisemarkersize = str2num(get(handles.edit_noisemarker,'string'));
hold on
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
plotLst=PMI{currentsub}.plotLst;
set(handles.edit_plotLst,'string',mat2str(plotLst))
dyoffset = -abs(str2num(get(handles.edit_stepvalue,'string')));
offset = 0;


bad_ch = find(~PMI{currentsub}.data(cf).MeasListAct);
ind_bad = ismember(PMI{currentsub}.plotLst,bad_ch');
if ~isempty(ind_bad)
    plot_bad = PMI{currentsub}.plotLst(ind_bad);
else
    plot_bad = [];
end
%plot_bad = PMI{currentsub}.plotLst_bad;%A faire disparaitre toujours utiliser MeasListAct

void_ch = NaN*ones(size(PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(1))));

%Trouver le champ normlist s'il  existe pour afficher les canaux à
%renormaliser en gras
imodulenorm = [];
for i = 1:numel(handles.NIRS.Dt.fir.pp)
    if numel(handles.NIRS.Dt.fir.pp(i).pre) >= 14
        if strcmp(handles.NIRS.Dt.fir.pp(i).pre(1:14),'Step Detection')
            imodulenorm = i;        %Le dernier is plusieur Step Detection 
        end
    end
end  

idfile = get(handles.popupmenu_file,'value');
if ~isempty(imodulenorm)
    NormList = handles.NIRS.Dt.fir.pp(imodulenorm).normlist{idfile,1};
else
    NormList = [];
end

%Timing option
start = str2num(get(handles.edit_start,'string'));
duration = str2num(get(handles.edit_duration,'string'));
idmodulehold = get(handles.popupmenu_module_hold, 'value');
if strcmp(handles.NIRS.Dt.fir.pp(idmodulehold).pre,'READ_RAW_NIRS')&get(handles.radiobutton_hold,'value')
currentsub = 2;
end
try
    if isempty(start)|start == 0;
        xlim([PMI{currentsub}.data(cf).HRF.tHRF(1),PMI{currentsub}.data(cf).HRF.tHRF(end)]);
        xlabel('time (s)');
        set(handles.edit_start,'string','0');
        idstart = 1 ;
        idstop = numel(PMI{currentsub}.data(cf).HRF.tHRF);
    else
        idstart = find(start>PMI{currentsub}.data(cf).HRF.tHRF);
        if isempty(idstart)
            idstart = 1;
            set(handles.edit_start,'string',num2str(start));
        end
        
        if isempty(duration)
            duration = 20;
            set(handles.edit_duration,'string','20')
        end
        idstop = find((duration+start)>PMI{currentsub}.data(cf).HRF.tHRF);
        
        idstart = idstart(end);
        idstop  = idstop(end);
        if idstart >= idstop
            start = PMI{currentsub}.data(cf).HRF.tHRF(idstop) - duration;
            set(handles.edit_start,'string',num2str(start));
            idstart = find(start>PMI{currentsub}.data(cf).HRF.tHRF);
            idstart = idstart(end);
        end
        xlim([PMI{currentsub}.data(cf).HRF.tHRF(idstart),PMI{currentsub}.data(cf).HRF.tHRF(idstop)]);
    end
catch
    1;
    h=msgbox('Time invalid');
    uiwait(h)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mode 1 NORMAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.popup_view,'value') ==1 %Visualisation des canaux selectionnés uniquements
    for i=1:numel(plotLst)
        try
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, plotLst(i),...
            numel(PMI{currentsub}.color)/3);
        catch
            idxc = 1;
        end
        %Affichage courbe
        srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),2));
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(i));
        if isempty(find(plot_bad == plotLst(i))) %Normal plot (bad channels are not displayed)
            h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idstart:idstop),d1(idstart:idstop),'linewidth',parameter_linewidth );
            %Affichage error bar
            if get(handles.popup_average,'value')==2
                AvgStdErr = PMI{currentsub}.data(cf).HRF.AvgStdErr(:,plotLst(i));
                he = errorbar(PMI{currentsub}.data(cf).HRF.tHRF, d1,AvgStdErr );
                set(he,'color',PMI{currentsub}.color(idxc,:));
                set(he,'Tag',num2str(plotLst(i)))
                set(he,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
                if ~newfigure
                    set(he,'uicontextmenu',handles.remove_bad_ch)
                end
            end
        elseif get(handles.radio_showrejectedchannels,'value') %Show bad channels (bad channels are in dotted line)
            h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idstart:idstop),d1(idstart:idstop),'LineStyle',':','linewidth',parameter_linewidth);
            %Affichage error bar
            if get(handles.popup_average,'value')==2
                AvgStdErr = PMI{currentsub}.data(cf).HRF.AvgStdErr(:,plotLst(i));
                he = errorbar(PMI{currentsub}.data(cf).HRF.tHRF, d1,AvgStdErr );
                set(he,'color',PMI{currentsub}.color(idxc,:));
                set(he,'Tag',num2str(plotLst(i)))
                set(he,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
                if ~newfigure
                    set(he,'uicontextmenu',handles.remove_bad_ch)
                end
            end
        else
            h = plot(void_ch,'linewidth',parameter_linewidth);
        end
        set(h,'color',PMI{currentsub}.color(idxc,:));
        set(h,'Tag',num2str(plotLst(i)));
        set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end
    end    
    for i=1:numel(plotLst) 
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, plotLst(i),...
            numel(PMI{currentsub}.color)/3);
        %Affichage courbe
        srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),2));
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(i));
        
        
        ch1 = plotLst(i);
        if ch1 <= handles.NIRS.Cf.H.C.N/2
            ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
        else
            ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
        end
        list = sort([ch1,ch2]);
        ch1 = list(1);
        if ~isempty(NormList)
            ind_on = find(NormList(:,1)==ch1&NormList(:,5)==1); %Trouver les indices du canal à on
        else
            ind_on = [];
        end
           %Afficher en gras les parties à renormaliser
        if ~isempty(ind_on)
            try
            for i_on=1:numel(ind_on)
               if isempty(find(plot_bad == NormList(ind_on(i_on),1)))
                i_start = NormList(ind_on(i_on),2);
                i_stop = NormList(ind_on(i_on),3);
                t = downsample(PMI{currentsub}.data(cf).HRF.tHRF(i_start:i_stop),60);
                x = downsample(d1(i_start:i_stop),60);
                h3 = plot(t,x,'Tag',num2str(plotLst(i)),'line','none','linewidth',2,'color','k','marker','v');   
                end
            end
            catch
            end
        end

        if isfield(PMI{currentsub}.data(cf).HRF,'noise')
            idnoise = find(PMI{currentsub}.data(cf).HRF.noise(:,plotLst(i)));
           if noisemarkersize>0
            if get(handles.radio_showrejectedchannels,'value') %Show markers for every channels in plotLst
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise),'*y','markersize',noisemarkersize);
                end
            elseif isempty(find(plot_bad == plotLst(i))) %Show markers only for good channels
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise),'*y','markersize',noisemarkersize);
                end
            end
           end
        end
    end
    %Line zero
    plot([PMI{currentsub}.data(cf).HRF.tHRF(1),PMI{currentsub}.data(cf).HRF.tHRF(end)],[0,0],'k')
%HOLD ON NORMAL VIEW 2 different step simultenously
if get(handles.radiobutton_hold,'value')
    currentsub = 2;  
     for i=1:numel(plotLst)         
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, plotLst(i),...
            numel(PMI{currentsub}.color)/3);
        %Affichage courbe
        srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),2));
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(i));
        if isempty(find(plot_bad == plotLst(i))) %Normal plot (bad channels are not displayed)
            h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1,'linewidth',parameter_linewidth );
            %Affichage error bar
            if get(handles.popup_average,'value')==2
                AvgStdErr = PMI{currentsub}.data(cf).HRF.AvgStdErr(:,plotLst(i));
                he = errorbar(PMI{currentsub}.data(cf).HRF.tHRF, d1,AvgStdErr );
                set(he,'color',PMI{currentsub}.color(idxc,:));
                set(he,'Tag',num2str(plotLst(i)))
                set(he,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
                if ~newfigure
                    set(he,'uicontextmenu',handles.remove_bad_ch)
                end
            end
        elseif get(handles.radio_showrejectedchannels,'value') %Show bad channels (bad channels are in dotted line)
            h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1,'LineStyle',':','linewidth',parameter_linewidth);
            %Affichage error bar
            if get(handles.popup_average,'value')==2
                AvgStdErr = PMI{currentsub}.data(cf).HRF.AvgStdErr(:,plotLst(i));
                he = errorbar(PMI{currentsub}.data(cf).HRF.tHRF, d1,AvgStdErr );
                set(he,'color',PMI{currentsub}.color(idxc,:));
                set(he,'Tag',num2str(plotLst(i)))
                set(he,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
                if ~newfigure
                    set(he,'uicontextmenu',handles.remove_bad_ch)
                end
            end
        else
            h = plot(void_ch,'linewidth',parameter_linewidth);
        end
        set(h,'color',PMI{currentsub}.color(idxc,:));
        set(h,'Tag',num2str(plotLst(i)));
        set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det,'overlay']);
        set(h,'linewidth',0.5)
        set(h,'LineStyle','--')

        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end        
    end  
    for i=1:numel(plotLst) 
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, plotLst(i),...
            numel(PMI{currentsub}.color)/3);
        %Affichage courbe
        srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),2));
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(i));
        
        
        ch1 = plotLst(i);
        if ch1 <= handles.NIRS.Cf.H.C.N/2
            ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
        else
            ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
        end
        list = sort([ch1,ch2]);
        ch1 = list(1);
        if ~isempty(NormList)
            ind_on = find(NormList(:,1)==ch1&NormList(:,5)==1); %Trouver les indices du canal à on
        else
            ind_on = [];
        end
           %Afficher en gras les parties à renormaliser
        if ~isempty(ind_on)
            try
            for i_on=1:numel(ind_on)
               if isempty(find(plot_bad == NormList(ind_on(i_on),1)))
                i_start = NormList(ind_on(i_on),2);
                i_stop = NormList(ind_on(i_on),3);
                t = downsample(PMI{currentsub}.data(cf).HRF.tHRF(i_start:i_stop),60);
                x = downsample(d1(i_start:i_stop),60);
                h3 = plot(t,x,'Tag',num2str(plotLst(i)),'line','none','linewidth',2,'color','k','marker','v');   
                end
            end
            catch
            end
        end

        if isfield(PMI{currentsub}.data(cf).HRF,'noise')
            idnoise = find(PMI{currentsub}.data(cf).HRF.noise(:,plotLst(i)));
            if noisemarkersize>0
            if get(handles.radio_showrejectedchannels,'value') %Show markers for every channels in plotLst
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise),'*y','markersize',noisemarkersize);
                end
            elseif isempty(find(plot_bad == plotLst(i))) %Show markers only for good channels
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise),'*y','markersize',noisemarkersize);
                end
            end
            end
        end
    end
    %Line zero
    plot([PMI{currentsub}.data(cf).HRF.tHRF(1),PMI{currentsub}.data(cf).HRF.tHRF(end)],[0,0],'k')
     currentsub = 1;
end

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mode 2 Normal (step) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif get(handles.popup_view,'value')==2
    for i=1:numel(plotLst)
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, plotLst(i),...
            numel(PMI{currentsub}.color)/3);
        
        ch1 = plotLst(i);
        if ch1 <= handles.NIRS.Cf.H.C.N/2
            ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
        else
            ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
        end
        list = sort([ch1,ch2]);
        ch1 = list(1);
        if ~isempty(NormList)
            ind_on = find(NormList(:,1)==ch1&NormList(:,5)==1); %Trouver les indices du canal à on
        else
            ind_on = [];
        end
        
        
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(i))+offset;
        srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),2));
        
        if isempty(find(plot_bad == plotLst(i))) %Show markers for every channels in plotLst
            h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1,'Tag',num2str(plotLst(i)),'linewidth',parameter_linewidth);
        elseif get(handles.radio_showrejectedchannels,'value')%Show bad channels (bad channels are in dotted line)
            h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1,'Tag',num2str(plotLst(i)),'LineStyle',':','linewidth',parameter_linewidth);
        else
            h = plot(void_ch,'linewidth',parameter_linewidth);
        end
        set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end
        set(h,'color',PMI{currentsub}.color(idxc,:));
        %Afficher en gras les parties à renormaliser
        if ~isempty(ind_on)
            for i_on=1:numel(ind_on)
                i_start = NormList(ind_on(i_on),2);
                i_stop = NormList(ind_on(i_on),3);
                t = downsample(PMI{currentsub}.data(cf).HRF.tHRF(i_start:i_stop),60);
                x = downsample(d1(i_start:i_stop),60);
                h3 = plot(t,x,'Tag',num2str(plotLst(i)),'line','none','linewidth',2,'color','k','marker','v');   
            end
        end
        set(h,'Tag',num2str(plotLst(i)))
        set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end
        
        
        if isfield(PMI{currentsub}.data(cf).HRF,'noise')
            idnoise = find(PMI{currentsub}.data(cf).HRF.noise(:,plotLst(i)));
            if get(handles.radio_showrejectedchannels,'value') %Show markers for every channels in plotLst
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise),'*y','markersize',noisemarkersize);
                end
            elseif isempty(find(plot_bad == plotLst(i))) %Show markers only for good channels
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise),'*y','markersize',noisemarkersize);
                end
            end
        end
        set(h,'Tag',num2str(plotLst(i)))
        set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end
        offset = offset + dyoffset;
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mode 3 two Wavelenght %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif get(handles.popup_view,'value')==3 %HbO HbR or 830 690 nm sans le step avec dotted pour hbr
    all  =[];
    for i=1:numel(plotLst)
        ch1 = plotLst(i);
        if ch1 <= handles.NIRS.Cf.H.C.N/2
            ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
        else
            ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
        end
        if isempty(find(all==ch1))|isempty(find(all==ch2))
            all = [all,ch1,ch2];
        end
    end
    [plotLst] = sort([all]);
    set(handles.edit_plotLst,'string',mat2str(plotLst));
    
    for i=1:numel(plotLst)
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, plotLst(i),...
            numel(PMI{currentsub}.color)/3);
        %Affichage courbe
        srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),2));
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(i));
        if PMI{currentsub}.data.MeasListAct(plotLst(i)) %Normal plot (bad channels are not displayed)
            if plotLst(i) <= handles.NIRS.Cf.H.C.N/2
                offset = 0;
                modeHbO = 1;
            else
                offset =0;
                modeHbO = 0;
            end
            if modeHbO
                h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1+offset,'linewidth',parameter_linewidth);
            else
                h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1+offset,'linestyle','--','linewidth',parameter_linewidth);
            end
            set(h,'Tag',num2str(plotLst(i)));
            set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
            if ~newfigure
                set(h,'uicontextmenu',handles.remove_bad_ch)
            end
        %Affichage error bar
        if get(handles.popup_average,'value')==2
                AvgStdErr = PMI{currentsub}.data(cf).HRF.AvgStdErr(:,plotLst(i));
                he = errorbar(PMI{currentsub}.data(cf).HRF.tHRF, d1+offset,AvgStdErr );
                set(he,'color',PMI{currentsub}.color(idxc,:));
                set(he,'Tag',num2str(plotLst(i)))
                set(he,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
                if ~newfigure
                    set(he,'uicontextmenu',handles.remove_bad_ch)
                end
            end
        elseif get(handles.radio_showrejectedchannels,'value') %Show bad channels (bad channels are in dotted line)
            if plotLst(i) <= handles.NIRS.Cf.H.C.N/2
                offset = 0;
            else
                offset = dyoffset;
            end
            h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1+offset,'LineStyle',':','linewidth',parameter_linewidth);
            %Affichage error bar
            if get(handles.popup_average,'value')==2
                AvgStdErr = PMI{currentsub}.data(cf).HRF.AvgStdErr(:,plotLst(i))+offset;
                he = errorbar(PMI{currentsub}.data(cf).HRF.tHRF, d1,AvgStdErr );
                set(he,'color',PMI{currentsub}.color(idxc,:));
                set(he,'Tag',num2str(plotLst(i)))
                set(he,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
                if ~newfigure
                    set(he,'uicontextmenu',handles.remove_bad_ch)
                end
            end
        else
            h = plot(void_ch,'linewidth',parameter_linewidth);
        end
        set(h,'color',PMI{currentsub}.color(idxc,:));
        set(h,'Tag',num2str(plotLst(i)));
        set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end
        
        
        if isfield(PMI{currentsub}.data(cf).HRF,'noise')
            idnoise = find(PMI{currentsub}.data(cf).HRF.noise(:,plotLst(i)));
            if get(handles.radio_showrejectedchannels,'value') %Show markers for every channels in plotLst
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise)+offset,'*y','markersize',noisemarkersize);
                end
            elseif isempty(find(plot_bad == plotLst(i))) %Show markers only for good channels
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise)+offset,'*y','markersize',noisemarkersize);
                end
            end
        end
        set(h,'Tag',num2str(plotLst(i)));
        set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mode 4 two Wavelenght (step) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif get(handles.popup_view,'value')==4 %830 et 690 avec step between
    all  =[];
    for i=1:numel(plotLst)
        ch1 = plotLst(i);
        if ch1 <= handles.NIRS.Cf.H.C.N/2
            ch2 = ch1 + handles.NIRS.Cf.H.C.N/2;
        else
            ch2 = ch1 - handles.NIRS.Cf.H.C.N/2;
        end
        if isempty(find(all==ch1))|isempty(find(all==ch2))
            all = [all,ch1,ch2];
        end
    end
    [plotLst] = sort([all]);
    set(handles.edit_plotLst,'string',mat2str(plotLst));
    
    for i=1:numel(plotLst)
        if plotLst(i) <= handles.NIRS.Cf.H.C.N/2
            offset = 0;
            modeHbO = 1;
        else
            offset =0;
            modeHbO = 0;
        end
        idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList, plotLst(i),...
            numel(PMI{currentsub}.color)/3);
        %Affichage courbe
        srs = SDPairs2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),1));
        det = SDDet2strboxy(PMI{currentsub}.data(cf).MeasList(plotLst(i),2));
        d1 = PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(i));
        if PMI{currentsub}.data.MeasListAct(plotLst(i)) %Normal plot (bad channels are not displayed)
            if plotLst(i) <= handles.NIRS.Cf.H.C.N/2
                offset = 0;
            else
                offset =dyoffset;
            end
            if modeHbO
                h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1+offset,'linewidth',parameter_linewidth);
            else
                h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1+offset,'linestyle','--','linewidth',parameter_linewidth);
            end
            set(h,'Tag',num2str(plotLst(i)));
            set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
            if ~newfigure
                set(h,'uicontextmenu',handles.remove_bad_ch)
            end
            %Affichage error bar
            if get(handles.popup_average,'value')==2
                AvgStdErr = PMI{currentsub}.data(cf).HRF.AvgStdErr(:,plotLst(i));
                he = errorbar(PMI{currentsub}.data(cf).HRF.tHRF, d1+offset,AvgStdErr );
                set(he,'color',PMI{currentsub}.color(idxc,:));
                set(he,'Tag',num2str(plotLst(i)))
                set(he,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
                if ~newfigure
                    set(he,'uicontextmenu',handles.remove_bad_ch)
                end
            end
        elseif get(handles.radio_showrejectedchannels,'value') %Show bad channels (bad channels are in dotted line)
            if plotLst(i) <= handles.NIRS.Cf.H.C.N/2
                offset = 0;
            else
                offset = dyoffset;
            end
            h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1+offset,'LineStyle',':','linewidth',parameter_linewidth);
            %Affichage error bar
            if get(handles.popup_average,'value')==2
                AvgStdErr = PMI{currentsub}.data(cf).HRF.AvgStdErr(:,plotLst(i));
                he = errorbar(PMI{currentsub}.data(cf).HRF.tHRF, d1+offset,AvgStdErr );
                set(he,'color',PMI{currentsub}.color(idxc,:));
                set(he,'Tag',num2str(plotLst(i)))
                set(he,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
                if ~newfigure
                    set(he,'uicontextmenu',handles.remove_bad_ch)
                end
            end
        else
            h = plot(void_ch,'linewidth',parameter_linewidth);
        end
        set(h,'color',PMI{currentsub}.color(idxc,:));
        set(h,'Tag',num2str(plotLst(i)));
        set(h,'Displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end
        
        
        if isfield(PMI{currentsub}.data(cf).HRF,'noise')
            idnoise = find(PMI{currentsub}.data(cf).HRF.noise(:,plotLst(i)));
            if get(handles.radio_showrejectedchannels,'value') %Show markers for every channels in plotLst
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise)+offset,'*y','markersize',noisemarkersize);
                end
            elseif isempty(find(plot_bad == plotLst(i))) %Show markers only for good channels
                if ~isempty(idnoise)
                    h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise)+offset,'*y','markersize',noisemarkersize);
                end
            end
        end
        set(h,'Tag',num2str(plotLst(i)));
        set(h,'displayname',['ch',num2str(plotLst(i)),'_',srs, '_', det]);
        if ~newfigure
            set(h,'uicontextmenu',handles.remove_bad_ch)
        end
    end
    plot([PMI{currentsub}.data(cf).HRF.tHRF(1),PMI{currentsub}.data(cf).HRF.tHRF(end)],[0,0],'k')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mode 5 zone list step %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif get(handles.popup_view,'value')==5 %Visualisation des zones regroupés avec une couleur pour chaque zone
    firsttime = 1;
    if ~isfield(PMI{currentsub},'zone')
        set(handles.popup_view,'value',1)
        set(handles.edit_zonelist,'enable','off')
        handles.newlist=0;
        msgbox('Please create zone before using this mode')
        return
    end
    izonelist = str2num(get(handles.edit_zonelist,'string'));    
    
    if isempty(izonelist)
        izonelist = 1:numel(PMI{currentsub}.zone.plotLst);
    end   
    plotLstall = [];
    for izone=izonelist;
        plotLst = [];
        color = [];
        try
            plotLst =[plotLst;PMI{currentsub}.zone.plotLst{izone}];
        catch
            msgbox('Please enter an integer as the identification number of the zone you want to display')
            return
        end
        for izoneplotLst=1:numel(PMI{currentsub}.zone.plotLst{izone});
            color = [color;PMI{currentsub}.zone.color(izone,:)];
        end
        %Line zero
        
        plot([PMI{currentsub}.data(cf).HRF.tHRF(1),PMI{currentsub}.data(cf).HRF.tHRF(end)],[offset,offset],'k')
        for i=1:numel(plotLst)
           % PMI{1}.color(plotLst(i),:)=color(i,:); %affichage de couleur des canaux = couleur des zones
            d1 = PMI{currentsub}.data(cf).HRF.AvgC(:,plotLst(i))+ offset;
            if isempty(find(plot_bad == plotLst(i))) %Normal plot (bad channels are not displayed)
                h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1,'Tag',num2str(plotLst(i)),'linewidth',parameter_linewidth);
            elseif get(handles.radio_showrejectedchannels,'value')  %Show bad channels (bad channels are in dotted line)
                h = plot(PMI{currentsub}.data(cf).HRF.tHRF,d1,'Tag',num2str(plotLst(i)),'LineStyle',':','linewidth',parameter_linewidth);
            else
                h = plot(void_ch,'linewidth',parameter_linewidth);
            end
            if 0==mod(izone,2)
                set(h,'linestyle','--')
            end
            plotLstall = [plotLstall, plotLst(i)];
            name = ['zone',num2str(izone),'ch',num2str(plotLst(i))];
            set(h,'color',color(i,:),'displayname',[name]);
            if ~newfigure
                set(h,'uicontextmenu',handles.context_Artefact)
            end
            if isfield(PMI{currentsub}.data(cf).HRF,'noise')
                idnoise = find(PMI{currentsub}.data(cf).HRF.noise(:,plotLst(i)));
                if get(handles.radio_showrejectedchannels,'value')  %Show markers for every channels in plotLst
                    if ~isempty(idnoise)
                        h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise),'*y','markersize',noisemarkersize);
                    end
                elseif isempty(find(plot_bad == plotLst(i))) %Show markers only for good channels
                    if ~isempty(idnoise)
                        h = plot(PMI{currentsub}.data(cf).HRF.tHRF(idnoise'),d1(idnoise),'*y','markersize',noisemarkersize);
                    end
                end
            end
        end
        if numel(plotLst)>1
            offset = offset + dyoffset;
        elseif firsttime==1
            firsttime =0;
            offset = offset + dyoffset;
        end
    end
    set(handles.edit_plotLst,'string',mat2str(plotLstall))
    id = zeros(numel(PMI{currentsub}.data.MeasListAct),1);
    id(plotLstall)=1;
    if numel(plotLstall)< 1
        plotLst = find(plotLstall(1:end/2)| plotLstall((end/2+1):end));
    else
        plotLst = plotLstall;
    end
    PMI{currentsub}.plotLst=plotLst;
    PMI{currentsub}.plot = [PMI{currentsub}.data(cf).MeasList(plotLst,1),PMI{currentsub}.data(cf).MeasList(plotLst,2)];
    set(guiHOMER,'UserData',PMI)        
    %Line zero
    plot([PMI{currentsub}.data(cf).HRF.tHRF(1),PMI{currentsub}.data(cf).HRF.tHRF(end)],[0,0],'k')
    
end


%%%%%%%%%%%%%%%%%% GENERAL ALL MODE AJUSTE Y axis%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ajuste Y axis
ylim('auto')
if ~get(handles.radio_autoaxis,'value')
    ymin = str2num(get(handles.edit_min,'string'));
    ymax = str2num(get(handles.edit_max,'string'));
    if ymin < ymax
        set(gca,'YLIM',[ymin,ymax])
    else
        yli = get(gca,'YLim');
        set(handles.edit_min,'string',num2str(yli(1)));
        set(handles.edit_max,'string',num2str(yli(2)));
    end
end

ylabel('');
label = 0;
for i = 1:size(handles.module_list,2)
        if strncmp('ModifyBeerLambertLaw',handles.module_list{1,i},numel('ModifyBeerLambertLaw')) == 1
            if handles.popupmenu_module.Value < i
                ylabel('Light intensity (a.u.)')
            else
                ylabel('dCONC (uM/L)')
            end
            label = 1;
            break            
        end
    end
   

if label == 0
    ylabel('Light intensity (a.u.)')
end


    



%%%%%%%%%%%%%%%%%% GENERAL TIMING OPTION %%%%%%%%%%%%%%%%%%%%%%%%%%
yli = get(gca,'YLim');
posxstart = str2num(get(handles.edit_time_start,'string'));
posxstop = str2num(get(handles.edit_time_stop,'string'));
if ~isempty(posxstart)
    h = plot([posxstart,posxstart],[yli(1), yli(2)],'linewidth',4,'color','b');
end
if ~isempty(posxstop)
    h = plot([posxstop,posxstop],[yli(1), yli(2)],'linewidth',4,'color','r');
end
 
%Plot triggers


if isfield(handles,'selected_trig')
    idfile = get(handles.popupmenu_file,'value');
    idmodulehold = get(handles.popupmenu_module_hold, 'value');
%     if strcmp(handles.NIRS.Dt.fir.pp(idmodulehold).pre,'READ_RAW_NIRS')&get(handles.radiobutton_hold,'value')
%         selected_trig=handles.NIRS.Dt.fir.aux5ini{idfile};
%     else
%         selected_trig = handles.selected_trig;
%     end
    try 
    tmp = [];
    trig_list=get(handles.popupmenu8, 'String');
    for itrig = 1:numel(trig_list)
       trigtype =trig_list{itrig};
       numtrig = str2num(trigtype(1:end-4));
        if  strcmp(trigtype(end-3:end),'  on')
         tmp = [tmp;handles.triggers(ismember(handles.triggers(:,1), numtrig),:)];
        end
        
        %handles.selected_trig(itrig,:) = handles.triggers(ismember(handles.triggers(:,1), str2num(trig_name{trig_ind(itrig)})),:);
    end
   selected_trig = tmp;
    catch
    end
    
    
    if get(handles.radio_triggers,'value')
        
        for i = 1:size(selected_trig,1)
            try
            trig_val = selected_trig(i,1);
            trig_ind = selected_trig(i,2);
            if trig_ind <=numel(PMI{currentsub}.data(cf).HRF.tHRF)
               trig_ind = PMI{currentsub}.data(cf).HRF.tHRF(trig_ind);
               if isfield(handles.NIRS.Dt.fir,['auxtrig',num2str(trig_val)]) 
                   if isfield(eval(['handles.NIRS.Dt.fir.','auxtrig',num2str(trig_val)]),'color')                       
                        colornew = eval(['handles.NIRS.Dt.fir.','auxtrig',num2str(trig_val),'.color']);
                   else
                       colornew = [0,1,0];
                   end
                   if isfield(eval(['handles.NIRS.Dt.fir.','auxtrig',num2str(trig_val)]),'.label')
                        labelnew= eval(['handles.NIRS.Dt.fir.','auxtrig',num2str(trig_val),'.label']);       
                   else 
                       labelnew=  num2str(selected_trig(i,1));
                   end
                   h = plot([trig_ind, trig_ind],[yli(1), yli(2)],'displayname',labelnew, 'color',colornew,'uicontextmenu',handles.Context_removetrigger);   
                   text(trig_ind,yli(2),  labelnew)
               else
                   h = plot([trig_ind, trig_ind],[yli(1), yli(2)], 'g','displayname',num2str(selected_trig(i,1)),'uicontextmenu',handles.Context_removetrigger);
                   text(trig_ind,yli(2), num2str(selected_trig(i,1)));
               end
                 
            end
               catch
            end
        end
    end
else
    %msgbox('NO TRIGGER AVAILABLE')
end
if 0 %get(handles.radio_3dmontageupdate,'value')
    option.conc = 1;
    option.SPMgui = 1;
    option.restorestetting = 0;
    IO_HelmetMTG_Display(handles, option)
end
idmodule = get(handles.popupmenu_module, 'value');
if idmodule>numel(handles.NIRS.Dt.fir.pp)
    set(handles.popupmenu_module,'value',numel(handles.NIRS.Dt.fir.pp));
    idmodule = numel(handles.NIRS.Dt.fir.pp);
end
fileid=get(handles.popupmenu_file,'value');
    nbhalf=numel(PMI{currentsub}.data(cf).MeasListAct)/2;
    ch1 = find(PMI{currentsub}.data(cf).MeasListAct);%(1);
    idfirst = find(ch1<=nbhalf);
    idlast = find(ch1>nbhalf);
    halflist = zeros(nbhalf,2);
    halflist(ch1(idfirst),:)=1;
    halflist(ch1(idlast)-nbhalf,:)=1;  
if ~isempty( strfind(handles.NIRS.Dt.fir.pp(idmodule).pre, 'Epoch averaging')); 
    handles.NIRS.Cf.H.C.okavg(:,:) = 0;
    handles.NIRS.Cf.H.C.okavg(find(halflist(:)),:)=1 ;
else
    handles.NIRS.Cf.H.C.ok(:,fileid)=0;
    handles.NIRS.Cf.H.C.ok(find(halflist(:)),fileid) = 1 ;
end

NIRS = handles.NIRS;
save(handles.NIRSpath{handles.subjectnb,1},'NIRS','-mat');

set(guiHOMER,'UserData',PMI)