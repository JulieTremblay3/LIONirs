function handles=plot_axesAij(varargin)
%AFFICHAGE DES MATRICES DE CONNECTIVITE PAR ZONE DANS L'INTERFACE
% varargin1 handles
% varargin2 newfigure
 
handles = varargin{1};
if numel(varargin)>1; newfigure = varargin{2};else; newfigure = 0;end 
if newfigure ==0
    axes(handles.axes_AJD);cla;hold on
else
    figure;hold on
    colormap(jet)
end
axis ij
DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value');
cmax = str2num(get(handles.edit_cmax,'string'));
cmin = str2num(get(handles.edit_cmin,'string'));
MAT = DATA{id}.MAT;
%Ordoner en ordre de zone
if get(handles.radio_fisher,'value')
    MAT =1/2*(log((1+MAT )./(1-MAT )));
end
List = strvcat(DATA{id}.ZoneList);
ML = DATA{id}.zone.ml;
idzone = [];
idlist = [];
idlabelall = [];
%POUR IDENTIFICATION DE TOUTE LES ZONES
for izone = 1:numel(DATA{id}.zone.plotLst)
    chzone = DATA{id}.zone.plotLst{izone};
    for ichzone = 1:numel(chzone)
        ich = chzone(ichzone);
        switch  DATA{id}.System
            case 'ISS Imagent'
                strDet = SDDet2strboxy_ISS(ML(ich,2));
                strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                idch = strmatch([strDet, ' ',strSrs ],List,'exact');
            case 'NIRx'
                strDet = SDDet2strboxy(ML(ich,2));
                strSrs = SDPairs2strboxy(ML(ich,1));
                idch = strmatch([strDet, ' ',strSrs ],List,'exact');
            otherwise
                strDet = SDDet2strboxy_ISS(ML(ich,2));
                strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                idch = strmatch([strDet, ' ',strSrs ],List,'exact');
        end
        idlist = [idlist, idch];
        if ichzone==1
            idzone =[idzone, izone];
        else
            idzone =[idzone, 0];
        end
    end
    idlabelall = [idlabelall, {[DATA{id}.zone.label{izone}]}];
end
set(handles.listbox_zone,'string',idlabelall);
%FIN

%POUR L'AFFICHAGE DES ZONES SELECTIONNÉES CANAL PAR CANAL
if get(handles.popupmenu_view,'value')==1%view zone
    listok = get(handles.listbox_selectedzone,'string');
    idlist = [];
    idlabel=[];
    idzone =[];
    izonebelong = [];
    for ilistzone = 1:numel(listok)
        for izone = 1:numel(DATA{id}.zone.plotLst)
            chzone = DATA{id}.zone.plotLst{izone};
            labelzone = DATA{id}.zone.label{izone};
            x = strmatch({labelzone} , {listok{ilistzone}}, 'exact');
            if ~isempty(x)
%                 for ichzone = 1:numel(chzone)
%                     ich = chzone(ichzone);
%                     if strcmp(DATA{id}.System,'ISS')
%                         strDet = SDDet2strboxy_ISS(ML(ich,2));
%                         strSrs = SDPairs2strboxy_ISS(ML(ich,1));
%                         idch = strmatch([strDet, ' ',strSrs ],List,'exact');
%                     end
                    idch=DATA{id}.zone.chMAT{izone};
                    idlist = [idlist, idch];
%                     if ichzone==1
%                         idzone =[idzone, izone];
%                     else
                        idzone =[idzone,izone, zeros(1,numel(idch)-1)];
                        izonebelong  = [ izonebelong ,ones(1,numel(idch)).*izone];
%                     end
%             end
            idlabel = [idlabel, {[DATA{id}.zone.label{izone}, sprintf('_%03.0f',ilistzone)]}];
            end
        end
    end
    idline = [find(idzone)-0.5,numel(idzone)+0.5];
   
    %AFFICHAGE THRESHOLD 
    if ~isempty(idlist)
        tr = str2num(get(handles.edit_threshold,'string'));
        if isempty(tr)
            tr = 0;
        end
        idtr = find(abs(MAT)<tr);
        MAT(idtr)=0;
        h = imagesc(MAT(idlist,idlist));
        nbchbadbyzone = [];
        for izone=1:max(izonebelong(:))
            ch = find(izonebelong==izone);
            nbchbad= numel(find(sum(0==(MAT(idlist,idlist(ch))))> (0.70*numel(idlist)))) ;  
            if isempty(nbchbad)
                nbchbadbyzone = 0;
            end
            nbchbadbyzone = [nbchbadbyzone,nbchbad, numel(idlist(ch))];
        end 
        set(handles.editchnb,'string',num2str(nbchbadbyzone));
        
        if newfigure ==0
            
            for iexbutton = 1:numel(idlist)
                hplot=plot(-1,iexbutton,'>');             
                ich = idlist(iexbutton);
                %strDet = SDDet2strboxy_ISS(ML(ich,2));
                %strSrs = SDPairs2strboxy_ISS(ML(ich,1));                    
                %set(hplot,'Displayname',[num2str(ich),'_', strDet, ' ',strSrs ]);
                set(hplot,'Displayname',[num2str(ich),' ',DATA{id}.ZoneList{idlist(iexbutton)}]); %take from the zone list information 
                set(hplot,'uicontextmenu',handles.context_mat);
                
            end
        else
           title(DATA{id}.name)
        end
%         for id = 1:(numel(idline)-1)
%             plot([idline(id), idline(id)],[idline(id+1), idline(id)],'k','linewidth',4)
%             plot([idline(id+1), idline(id)],[idline(id), idline(id)],'k','linewidth',4)
%             plot([idline(id), idline(id+1)],[idline(id+1), idline(id+1)],'k','linewidth',4)
%             plot([idline(id+1), idline(id+1)],[idline(id+1), idline(id)],'k','linewidth',4)
%         end
    end
    
    caxis([cmin cmax]);
    set(gca,'xtick', find(idzone));
    set(gca,'xticklabel', idlabel);
    set(gca,'ytick', find(idzone));
    set(gca,'yticklabel', idlabel);
    colorbar
    %TEMP
    % subMAT.mat = MAT(idlist,idlist);
% subMAT.pos = DATA{id}.zone.pos(idlist,1:3)
% Mat = threshold_absolute(subMAT.mat , 0.5);
% [X,Y,Z] = adjacency_plot_und(Mat,subMAT.pos)
% figure;plot3(X,Y,Z); 
    if get(handles.radio_individualplot,'value')
        MATtrial = DATA{id}.MATtrial;
        lar =  ceil(size(MATtrial,3)/4);
        figure;
        for id=1:size(MATtrial,3)
            subplot(4,lar,id);
            tr = str2num(get(handles.edit_threshold,'string'));
            idtr = find(abs(MATtrial)<tr);
           % MATtrial(idtr)=nan;
            MATtrial(idtr)=0;
            h = imagesc(MATtrial( idlist, idlist,id));
            set(h,'uicontextmenu',handles.context_mat);
            caxis([cmin cmax]);                
            set(gca,'xtick', find(idzone));
            set(gca,'xticklabel', idlabel);
            set(gca,'ytick', find(idzone));
            set(gca,'yticklabel', idlabel);
            %title(DATA{id}.name)
            title(['Bloc nb:',num2str(id)])
            axis xy
        end
        
    end
    %POUR L'AFFICHAGE DES ZONES SELECTIONNÉES MOYENNES DES ZONES
elseif  get(handles.popupmenu_view,'value')==2%view avg zone
    listok = get(handles.listbox_selectedzone,'string');
    MATAVG = zeros(numel(listok));
    idlist = [];
    idlabel=[];
    idzone =[];
    for adji = 1:numel(listok)
        for adjj = 1:numel(listok)            
            labelzone = listok{adji};
            x = strmatch({labelzone} ,idlabelall, 'exact');            
            labelzone = listok{adjj};
            y = strmatch({labelzone} ,idlabelall, 'exact');
            if isempty(x)|isempty(y)
                msgbox('problem zone in subject');
            end
            chzone = DATA{id}.zone.plotLst{x};
            idlisti = [];
            for ichzone = 1:numel(chzone);
                ich = chzone(ichzone);
                if strcmp(DATA{id}.System,'ISS')
                    strDet = SDDet2strboxy_ISS(ML(ich,2));
                    strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                    idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                end
                idlisti = [idlisti, idch];
            end
            
            chzone = DATA{id}.zone.plotLst{y};
            idlistj = [];
            for ichzone = 1:numel(chzone)
                ich = chzone(ichzone);
                 if strcmp(DATA{id}.System,'ISS')
                    strDet = SDDet2strboxy_ISS(ML(ich,2));
                    strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                    idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                 end
                idlistj = [idlistj, idch];
            end
            if isempty(idlisti)|isempty(idlistj)
                MATAVG(adji,adjj)=nan;
            else
                temp = MAT(idlisti, idlistj);
                MATAVG(adji,adjj)= nanmean(temp(:));
            end
            
        end
    end
    tr = str2num(get(handles.edit_threshold,'string'));
    idtr = find(abs(MATAVG)<tr);
    MATAVG(idtr)=0;
    h = imagesc(MATAVG);
    if newfigure ==0
        set(h,'uicontextmenu',handles.context_mat);
    else
         title(DATA{id}.name)
    end
    idzone = 1:numel(listok);
    idlabel = listok;
    set(gca,'xtick', find(idzone));
    set(gca,'xticklabel', idlabel);
    set(gca,'ytick', find(idzone));
    set(gca,'yticklabel', idlabel);
    caxis([cmin cmax]);
    colorbar
end %POUR AFFICHAGE DES ZONES SELECTIONNÉES

