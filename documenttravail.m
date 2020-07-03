

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
             idch=DATA{id}.zone.chMAT{izone};
             idlist = [idlist, idch];
             idzone =[idzone,izone, zeros(1,numel(idch)-1)];
             izonebelong  = [ izonebelong ,ones(1,numel(idch)).*izone];
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
        MAT2=MAT(idlist,idlist);
        cmax=MAT2(1,1);
        cmin=MAT2(1,1);
        dim=size(MAT2);
    end
    
    %POUR L'AFFICHAGE DES ZONES SELECTIONNÉES MOYENNES DES ZONES
elseif get(handles.popupmenu_view,'value')==2%view avg zone
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
                msgbox('problem zone in subject')
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
                MATAVG(adji,adjj)=nan
            else
                temp = MAT(idlisti, idlistj)
                MATAVG(adji,adjj)= nanmean(temp(:))
            end
            
        end
    end
    tr = str2num(get(handles.edit_threshold,'string'));
    idtr = find(abs(MATAVG)<tr);
    MATAVG(idtr)=0;
    cmax = MATAVG(1,1);
    cmin = MATAVG(1,1);
    dim = size(MATAVG);
    MAT2 = MATAVG;
end 
for i=1:dim(1)
    for j=1:dim(2)
        if MAT2(i,j)< cmin
            cmin= MAT2(i,j);
        elseif MAT2(i,j)> cmax
            cmax = MAT2(i,j);
        end
    end
end

cmin = sprintf('%0.2g',cmin);
cmax = sprintf('%0.2g',cmax);
set(handles.edit_cmin,'string',[cmin]);
set(handles.edit_cmax,'string',[cmax]);
cmin
cmax
guidata(handles.GUI_LookMat, handles);
updateNetAllView(handles);