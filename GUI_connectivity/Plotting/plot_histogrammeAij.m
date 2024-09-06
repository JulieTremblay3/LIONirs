function handles=plot_histogrammeAij(varargin)
%use the handles of the gui to greate histogram of the connectivity data
% varargin1 handles
% varargin2 newfigure
 
handles = varargin{1};
if numel(varargin)>1; newfigure = varargin{2};else; newfigure = 0;end 
if newfigure ==0
    axes(handles.axes_Zone);cla;hold on
else
    figure;hold on
end

DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value')
cmax=str2num(get(handles.edit_cmax,'string'))
cmin=str2num(get(handles.edit_cmin,'string'))
MAT = DATA{id}.MAT;
 

%Ordoner en ordre de zone
if get(handles.popupmenu_viewhistogramme,'value')==1 %pearson mean (connectivity intensity ) 
    %and std on selected zone combination
    nameall = [];
    for id=1:numel(DATA) %trouver les indices de moyennes de groupes.
        nameall = [nameall, {DATA{id}.name}];
        groupeall(id)=DATA{id}.GR ;     
    end
    for igroupe = 1:max(groupeall)
        idG(igroupe,1) = strmatch(['AVG groupe ',num2str(igroupe)],nameall,'exact');   
    end
    %Selected zone
    listok = get(handles.listbox_selectedzone,'string');
    idlabelall  = [];
    for izone = 1:numel(DATA{id}.zone.plotLst)
        idlabelall = [idlabelall, {DATA{idG(1,1)}.zone.label{izone}}];
    end
    for idzone = 1:numel(listok)
        try
            idok = strmatch(upper(listok(idzone)),upper(idlabelall),'exact');   
            if ~isempty(idok)
                    idGzone(idzone) =idok
            else
                msgbox(['CONNECTOGRAMME ', listok(idzone),' not found'])
            end
        catch
            1
        end
    end
    id = 1;
    if isempty(listok)
        return
    end
    for i=1:numel(idGzone)
        j = 1;        
        while j<i | j==i
        idGzonelist(id,1) = idGzone(i);
        idGzonelist(id,2) = idGzone(j);
        labelzone{id} =  [DATA{idG(1,1)}.zone.label{idGzonelist(i)},'/' ,DATA{idG(1,1)}.zone.label{idGzonelist(j)}]
        j = j + 1;
        id = id+1; 
        end
    end
    meanbar = zeros(size(idGzonelist,1),max(groupeall))
    stdbar = zeros(size(idGzonelist,1),max(groupeall))
     for i = 1:size(idGzonelist,1)
         for j = 1:max(groupeall)
            meanbar(i,j)=mean(DATA{idG(j)}.MATtrial(idGzonelist(i,1),idGzonelist(i,2),: ));
            stdbar(i,j)=std(DATA{idG(j)}.MATtrial(idGzonelist(i,1),idGzonelist(i,2),: ));            
         end
     end    
    	hbar = bar(meanbar);
%     for j = 1:numel(hbar)
%         xd = get(hbar(j),'xdata')
%         errorbar(xd ,meanbar(:,j),stdbar(:,j),'x')        
%     end
    
    set(gca,'XTick',[1:size(idGzonelist,1)])
    set(gca,'XTickLabel',labelzone)
end    
if get(handles.popupmenu_viewhistogramme,'value')==2 %subject pearson mean (connectivity intensity ) 
    idsubject = get(handles.popup_listsujet, 'value')
    listok = get(handles.listbox_selectedzone,'string');
    idlabelall  = [];
    for izone = 1:numel(DATA{idsubject}.zone.plotLst)
        idlabelall = [idlabelall, {DATA{idsubject}.zone.label{izone}}];
    end
    for idzone = 1:numel(listok)
        idGzone(idzone) = strmatch(listok(idzone),idlabelall,'exact');   
    end
    id = 1;
    if isempty(listok)
        return
    end
    for i=1:numel(idGzone)
        j = 1;        
        while j<i | j==i
        idGzonelist(id,1) = idGzone(i);
        idGzonelist(id,2) = idGzone(j);
        labelzone{id} =  [DATA{id}.zone.label{idGzonelist(i)},'/' ,DATA{id}.zone.label{idGzonelist(j)}]
        j = j + 1;
        id = id+1; 
        end
    end
    meanbar = zeros(size(idGzonelist,1),1)
     for i = 1:size(idGzonelist,1)
          ch1 = DATA{idsubject}.zone.chMAT{idGzonelist(i,1)}
          ch2 = DATA{idsubject}.zone.chMAT{idGzonelist(i,2)}
          
          val = DATA{idsubject}.MAT(ch1,ch2 )
          
          idbad = find(val==0)
          val(idbad)=nan;
          meanbar(i,1)=nanmean(val(:));
     end    
    	hbar = bar(meanbar);
    
    set(gca,'XTick',[1:size(idGzonelist,1)])
    set(gca,'XTickLabel',labelzone)

end




end


