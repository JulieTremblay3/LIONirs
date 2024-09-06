function  [ listHBOch, listHBRch]= findzoneinlistgood(NIRS, zone, listgood,label)
%Description: Find a specific zone label in the zone definition. 
%As we are definding a listgood channel description for component find the
%index in the listgood channel description corresponding to the zone
%definition. 
ML_new= [NIRS.Cf.H.C.id(2:3,:)',...
        ones(numel(NIRS.Cf.H.C.wl(:)),1),...
        NIRS.Cf.H.C.wl(:)];
ML_listgood= [NIRS.Cf.H.C.id(2:3,listgood)',...
        ones(numel(NIRS.Cf.H.C.wl(listgood)),1),...
        NIRS.Cf.H.C.wl(listgood)'];
%Check if ML subject and ML zone are identical
 if sum( sum(zone.ml == ML_new)) == numel(ML_new(:))
    %use the zone directly 
    %find the zone label 
    for izone = 1:numel(zone.label) 
        if strcmp(zone.label{izone}, label)
        zone.label{izone}
        listHBOch = []
        listHBRch = []
        plotLst = zone.plotLst{izone}
        for iplotLst = 1:numel(plotLst)
            id = find(ML_new(plotLst(iplotLst),1)==ML_listgood(:,1) & ML_new(plotLst(iplotLst),2)==ML_listgood(:,2) & ML_listgood(:,4)==1)
            listHBOch = [listHBOch,id]
            id = find(ML_new(plotLst(iplotLst),1)==ML_listgood(:,1) & ML_new(plotLst(iplotLst),2)==ML_listgood(:,2) & ML_listgood(:,4)==2)
            listHBRch = [listHBRch,id]
        end 
        end
    end
 end
    
    
    