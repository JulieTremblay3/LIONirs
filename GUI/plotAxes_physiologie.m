function handles = plotAxes_physiologie(varargin)


parameter_linewidth = 2;
handles=varargin{1};
newfigure = 0;
if numel(varargin)>1
    newfigure = varargin{2};
end
if newfigure
    figure
else
    axes(handles.axes_physiologie)
end
%%%%%%%%%%%%%PLOT AUX HERE%%%%%%%
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

setcolor = lines;
cla;
hold on;
if isfield( PMI{currentsub}.data(cf),'AUX')  
    for ich=1:numel(PMI{currentsub}.data(cf).AUX.view)
        if PMI{currentsub}.data(cf).AUX.view{ich}
              fs = PMI{currentsub}.data(cf).AUX.fs{ich};
            daux = PMI{currentsub}.data(cf).AUX.data{ich};
            Label = PMI{currentsub}.data(cf).AUX.label{ich};
            taux=1/fs:1/fs:(1/fs*numel(daux));
            plot(taux,daux,'linewidth',2,'color',setcolor(ich,:),'displayname',Label )
        end
    end
end

%Timing option
start = str2num(get(handles.edit_start,'string'));
duration = str2num(get(handles.edit_duration,'string'));
    if isempty(start)|start == 0;
        xlim([PMI{currentsub}.data(cf).HRF.tHRF(1),PMI{currentsub}.data(cf).HRF.tHRF(end)]);
        set(handles.edit_start,'string',0)
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
    if ~newfigure
  %  set(gca,'yticklabel',{' '})
    set(gca,'xticklabel',{' '})
    end
 