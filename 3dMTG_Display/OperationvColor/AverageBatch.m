function AverageBatch(path,namesubject,time,pathout) 
for itime = time
    try
    allHbO = [];
    for iname=1:numel(namesubject);
        load([path,namesubject{iname,1},'.aviHbOTime   ', num2str(itime),'(s).vcolor'],'-mat');
        allHbO = [allHbO,vColor];
    end
    vColor = mean(allHbO')';
    [nameout]= [pathout,'HbOAvg',num2str(itime),'(s).vColor']
    save([nameout],'vColor','-mat')
    catch 
        'HbO not averaged'
    end
    try 
    allHbR = [];
    for iname=1:numel(namesubject);
        load([path,namesubject{iname,1},'.aviHbRTime   ', num2str(itime),'(s).vcolor'],'-mat');
        allHbR = [allHbR,vColor];
    end
    vColor = mean(allHbR')';
    [nameout]= [pathout,'HbRAvg',num2str(itime),'(s).vColor']
    save([nameout],'vColor','-mat')
    catch
        'HbR not averaged'
    end
    try 
    allHbT = [];
    for iname=1:numel(namesubject);
        load([path,namesubject{iname,1},'.aviHbTTime   ', num2str(itime),'(s).vcolor'],'-mat');
        allHbT = [allHbT,vColor];
    end
    vColor = mean(allHbT')';
    [nameout]= [pathout,'HbTAvg',num2str(itime),'(s).vColor']
    save([nameout],'vColor','-mat')
    catch 
        'HbT not averaged'
    end
end


