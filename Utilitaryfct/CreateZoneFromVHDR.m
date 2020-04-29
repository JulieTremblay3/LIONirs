%Create zone from vhdr eeg marker
function zone = CreateZoneFromVHDR(file)
%File vhdr to open and createzone
% [name,path]=uigetfile('.vhdr')
% file = [path,name];
[path,name,ext]=fileparts(file); 
fid=fopen(file)
Chnb = [];
EEGlist = [];
ListtoGroupe = [];
plotid = []; 
while ~feof(fid)
    line = fgetl(fid);
    if numel(line)>2
        if line(1:2)=='Ch';
        [rem,tok]= strtok(line,'=');
        Chnb = [Chnb,str2num(rem(3:end))];
        [rem,tok]=strtok(tok,'_');
        EEGlist = [EEGlist,{rem(2:end)}];
        if isempty(find(strcmp(ListtoGroupe,{rem(2:end)})));
            ListtoGroupe = [ListtoGroupe,{rem(2:end)}];
        end
        [rem,tok] = strtok(tok,'_');
        num_src = Strboxy2SDPairs(rem);
        [rem,tok] = strtok(tok,',');
        num_det = StrBoxy2SDDet(rem(2:end));
        plotid = [plotid;num_src,num_det];
       
        end
    end
end
fclose(fid);

zone = [];
ListtoGroupe = sort(ListtoGroupe);
for i=1:numel(ListtoGroupe)
    list = [];
    for ich = 1:numel(Chnb)        
        if strcmp(ListtoGroupe{i},EEGlist{ich})
            list = [list,ich];
        end
    end
    zone.color(i,:) = [0,0,0];
    zone.label{i}= ListtoGroupe{i};
    labelzone= ListtoGroupe{i};

    if strcmp(labelzone(1),'C')
        if strcmp(labelzone(end),'z')
            
        elseif mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [0/255,0/255,255/255]; %Left
        elseif ~mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [255/255,153/255,0/255]; %Right
        end
    elseif strcmp(labelzone(1),'F')|strcmp(labelzone(1),'A')
        if strcmp(labelzone(end),'z')
            
        elseif mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [112/255,16/255,128/255]; %Left
        elseif ~mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [255/255,0/255,0/255]; %Right
        end
    elseif strcmp(labelzone(1),'T')|strcmp(labelzone(1:2),'FT')
        if strcmp(labelzone(end),'z')
            
        elseif mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [0/255,255/255,255/255]; %Left
        elseif ~mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [204/255,0/255,0/255]; %Right
        end
    elseif strcmp(labelzone(1),'P')
         if strcmp(labelzone(end),'z')
            
        elseif mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [11/255,132/255,199/255];%Left
        elseif ~mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [255/255,153/255,200/255]; %Right
        end
    elseif strcmp(labelzone(1),'O')
         if strcmp(labelzone(end),'z')
            
        elseif mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [0/255,127/255,0/255];%Left
         elseif ~mod(str2num(labelzone(end)),2)
            zone.color(i,:) = [255/255,102/255,0/255]; %Right
        end        
    end

    zone.plotLst{i}=list;
    zone.plot{i} = plotid(list,:);
    
end
save([path,name(1:end-5),'EEG.zone'],'zone','-mat');        
