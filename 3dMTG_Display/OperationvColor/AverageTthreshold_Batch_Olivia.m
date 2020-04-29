
path = 'I:\IO\Tonotopie\video\Video_SonLeft\'
pathout = 'I:\IO\Tonotopie\video\Test\'
thresholdHbO = 1.5;
thresholdHbR = -1.5;
thresholdHbT = 1.5;
Timingtocheck = 1:30 % ouvrir les fichier contenant les temps suivent 
FlagTimeWindow = 1; %0 faire chaque temps, 1 faire une fenêtre avec tout les temps

namesubject = {'Tono02JT_viewBack';
    'Tono03JP_viewBack'
    'Tono05GB_viewBack'
    'Tono07ED_viewBack'
    'Tono09CTL_viewBack'
    'Tono10MC_viewBack'
    'Tono13IVR_viewBack'
    'Tono15AS_viewBack'
    'Tono16NN_viewBack'
    'Tono17HJB_viewBack';}

clear labeltime

if FlagTimeWindow ==1
    listwindowtime = 1;
else
    listwindowtime = 1:numel(Timingtocheck);
end
    
for iwindowtime = listwindowtime
    %try
    if FlagTimeWindow ==1
        labeltime = cell(numel(Timingtocheck),1)
        for i = 1:numel(Timingtocheck)
            labeltime{i}=sprintf('%03.0f',Timingtocheck(i))
        end
    else
        labeltime{1}=sprintf('%03.0f',Timingtocheck(iwindowtime))
    end
    data = []
    for iname=1:numel(namesubject);
        idok = [];
        allHbO = [];
        for itime=1:numel(labeltime)
            load([path,namesubject{iname,1},'.aviHbOTime   ', labeltime{itime},'(s).vcolor'],'-mat');
            allHbO = [allHbO,vColor];
            idok = [idok;find(vColor>thresholdHbO)];
            if ~isfield(data,'subjectsum')
                data.subjectsum =zeros(size(vColor),numel(namesubject))
            end
        end
          %  h=figure
            if FlagTimeWindow ==1
                allHbO = max(allHbO');                     
            end
            ind = find(allHbO);
            [n,xout] =hist(allHbO(ind),-10:10);           
            nallHbO(iname,:)=n;
            
         if ~isempty(idok)
                data.subjectsum(idok,iname) = 1;
         end      
    end
    vColor =  sum(data.subjectsum,2)/numel(namesubject);
    if FlagTimeWindow==1
        [nameout]= [pathout,'HbOT',num2str(thresholdHbO),'val',labeltime{1},'to',labeltime{end},'(s).vColor']
    else
        [nameout]= [pathout,'HbOT',num2str(thresholdHbO),'val',labeltime{1},'(s).vColor']
    end
    save([nameout],'vColor','-mat')
    %catch 
    %    'HbO not averaged'
   % end
   %HbR
    data = []
    for iname=1:numel(namesubject);
        idok = [];
        allHbR = [];
        for itime=1:numel(labeltime)
            load([path,namesubject{iname,1},'.aviHbRTime   ', labeltime{itime},'(s).vcolor'],'-mat');
            idok = [idok;find(vColor<thresholdHbR)];
             allHbR = [allHbR,vColor];
            if ~isfield(data,'subjectsum')
                data.subjectsum =zeros(size(vColor),numel(namesubject))
            end
        end
         if ~isempty(idok)
                data.subjectsum(idok,iname) = 1;
         end  
            if FlagTimeWindow ==1
                allHbR = min(allHbR');                     
            end
            ind = find(allHbR);
            [n,xout] =hist(allHbR(ind),-10:10);           
            nallHbR(iname,:)=n;
    end
            vColor =  sum(data.subjectsum,2)/numel(namesubject);
            
         
    if FlagTimeWindow==1
        [nameout]= [pathout,'HbRT',num2str(thresholdHbR),'val',labeltime{1},'to',labeltime{end},'(s).vColor']
    else
        [nameout]= [pathout,'HbRT',num2str(thresholdHbR),'val',labeltime{1},'(s).vColor']
    end
    save([nameout],'vColor','-mat')
    %HbT
    data = []
    for iname=1:numel(namesubject);
        idok = [];
        allHbT = [];
        for itime=1:numel(labeltime)
            load([path,namesubject{iname,1},'.aviHbTTime   ', labeltime{itime},'(s).vcolor'],'-mat');
            idok = [idok;find(vColor>thresholdHbT)];
            allHbT = [allHbT,vColor];
            if ~isfield(data,'subjectsum')
                data.subjectsum =zeros(size(vColor),numel(namesubject))
            end
        end
        
         if ~isempty(idok)
                data.subjectsum(idok,iname) = 1;
         end    
            if FlagTimeWindow ==1
                allHbT = max(allHbT');                     
            end
            ind = find(allHbT);
            [n,xout] =hist(allHbT(ind),-10:10);           
            nallHbT(iname,:)=n;
    end
    vColor =  sum(data.subjectsum,2)/numel(namesubject);
    if FlagTimeWindow==1
        [nameout]= [pathout,'HbTT',num2str(thresholdHbT),'val',labeltime{1},'to',labeltime{end},'(s).vColor']
    else
        [nameout]= [pathout,'HbTT',num2str(thresholdHbT),'val',labeltime{1},'(s).vColor']
    end
    save([nameout],'vColor','-mat')
    
   
end

figure
subplot(3,1,1)
imagesc(xout,1:numel(namesubject),nallHbO)
xlabel('T value')
ylabel('Subject')
title('HbO')
subplot(3,1,2)
imagesc(xout,1:numel(namesubject),nallHbR)
xlabel('T value')
ylabel('Subject')
title('HbR')
subplot(3,1,3)
imagesc(xout,1:numel(namesubject),nallHbT)
xlabel('T value')
ylabel('Subject')
if FlagTimeWindow==1
    title(['HbT', labeltime{1},'to',labeltime{end}])
else
    title(['HbT',labeltime{end}]);
end
%     try 
%     allHbR = [];
%     for iname=1:numel(namesubject);
%         load([path,namesubject{iname,1},'.aviHbRTime   ', labeltime,'(s).vcolor'],'-mat');
%         allHbR = [allHbR,vColor];
%     end
%     vColor = mean(allHbR')';
%     [nameout]= [pathout,'HbRAvg',labeltime,'(s).vColor']
%     save([nameout],'vColor','-mat')
%     catch
%         'HbR not averaged'
%     end
%     try 
%     allHbT = [];
%     for iname=1:numel(namesubject);
%         load([path,namesubject{iname,1},'.aviHbRTime   ', labeltime,'(s).vcolor'],'-mat');
%         allHbT = [allHbT,vColor];
%     end
%     vColor = mean(allHbT')';
%     [nameout]= [pathout,'HbTAvg',labeltime,'(s).vColor']
%     save([nameout],'vColor','-mat')
%     catch 
%         'HbT not averaged'
%     end
% 
% 
% 
