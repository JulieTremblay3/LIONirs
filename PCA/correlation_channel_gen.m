function zone = correlation_channel_gen(varargin)
%arg1 data (data time x ch)
%arg2 measlistact (noise)
%arg3 listch
%arg4 coefficient correlation
%arg5 wavelenght 
d1 = varargin{1};
measlistact = varargin{2};

if numel(varargin) > 2 
    MeasList= varargin{3};
else
    MeasList = ones(numel(measlistact),4);
end
if numel(varargin) > 3
    coefcorr = varargin{4};
else
    coefcorr = 0.8;
end

if numel(varargin)>4
    usewv = varargin{5}; %First wavelengt half or second wavelenght half
else
    usewv = 1;
end
%d1 time x nb ch
nbch = size(d1,2)/2;
matcorr = zeros(nbch,nbch);

if usewv == 1
    chok = find(measlistact(1:end/2));
    for i = 1:numel(chok)
        ch = chok(i);
        d1ok = d1(:,ch);
        for j = i:numel(chok)
            ch2 = chok(j);
            d2ok = d1(:,ch2);
             matcorr(ch,ch2) = corr(d1ok,d2ok); 
             matcorr(ch2,ch) = matcorr(ch,ch2);
        end
    end
elseif usewv == 2
    chok = find(measlistact(end/2+1:end))+nbch;
    for i = 1:numel(chok)
    ch = chok(i);
    d1ok = d1(:,ch);
    for j = i:numel(chok)
        ch2 = chok(j);
        d2ok = d1(:,ch2);
         matcorr(ch-nbch,ch2-nbch) = corr(d1ok,d2ok); 
         matcorr(ch2-nbch,ch-nbch) = matcorr(ch-nbch,ch2-nbch);
    end
    end
    chok = find(measlistact(end/2+1:end));
    
elseif usewv == 3

    chokHbO = find(measlistact(1:end/2));
    chokHbR = find(measlistact(end/2+1:end))+nbch;
    chok = chokHbO ;
    for i = 1:numel(chok)
        ch = chokHbO(i);
        chHbO =  chokHbO(i);
        chHbR =  chokHbR(i);
        %d1ok = (abs(d1(:,chHbO)) + abs(d1(:,chHbR)))/2;
        for j = i:numel(chok)
            ch2 = chokHbO(j);
            chHbO =  chokHbO(j);
            chHbR =  chokHbR(j);
            d2ok = (abs(d1(:,chHbO)) + abs(d1(:,chHbR)))/2;
             matcorr(ch,ch2) = corr(d1ok,d2ok); 
             matcorr(ch2,ch) = matcorr(ch,ch2);

        end
        maxamp(ch) = max([abs(d1ok(:));abs(d1ok(:))]);
    end
   
    chok = find(measlistact(end/2+1:end));


elseif usewv ==4
      chok = find(measlistact(1:end/2));
    for i = 1:numel(chok)
        ch = chok(i);
        d1ok = d1(:,ch);
        for j = i:numel(chok)
            ch2 = chok(j);
            d2ok = d1(:,ch2);
             matcorr(ch,ch2) = abs(corr(d1ok,d2ok)); 
             matcorr(ch2,ch) = matcorr(ch,ch2);
        end
    end


end
%figure;imagesc(matcorr)
%1
%load('C:\data\default\01_default_pre\mat_default20trial.mat')
%matcorr = matall

listch = ones(nbch ,1);
together830 = [];
togethernb = [];
meanmaxamp  = [];
for i = 1:numel(chok)
    ch = chok(i);
    if listch(ch)
        ind = find(matcorr(ch,:)>coefcorr);
        indok  = find(listch(ind));
        listch(ind) = 0;                
        together830 = [together830,{ind(indok)}];     
        togethernb = [togethernb,numel(ind(indok))];
        if usewv == 3
            if numel(indok)==1
                meanmaxamp = [meanmaxamp, 0];
            else
                meanmaxamp = [meanmaxamp, mean(maxamp(ind(indok)))*numel(indok)];
            end
        end
    end            
end

colorlist = [1 0 0;
    0,0,1
    1,0,1
    85/255,0,1
colorcube(numel(together830))];
if  usewv == 1 | usewv == 2 | usewv == 4 
    [nb,ind]= sort(togethernb,'descend');
    together830 = together830(ind);
    numel(together830);
end 
if  usewv == 3  
    [nb,ind]= sort(meanmaxamp,'descend');
    together830 = together830(ind);
    numel(together830);
end 

zone = []; 
zone.color =[];
idzone = 0;
for nb_zone = 1:numel(together830);
    if ~isempty(together830{nb_zone})    
    if 0 %numel(together830{nb_zone})>1
        zone.plot{nb_zone} = 0; 
        zone.plotLst{nb_zone} = together830{nb_zone}';
        zone.label{nb_zone} = num2str(nb_zone);   
        zone.color= [zone.color;[0,0,0]];
    elseif 0 %830 seulement
        zone.plot{nb_zone} = 0; %[DOT{currentsub}.data(1).MeasList(together830{nb_zone},1),DOT{currentsub}.data(1).MeasList(together830{nb_zone},2)];
        zone.plotLst{nb_zone} = together830{nb_zone}';
        zone.label{nb_zone} = num2str(nb_zone);
        zone.color= [zone.color; colorlist(nb_zone,:)] ;   
    else
        if numel(together830{nb_zone})
        %830
        idzone = idzone + 1;
        zone.plot{idzone} = [MeasList(together830{nb_zone},1),MeasList(together830{nb_zone},2)];
        zone.plotLst{idzone} = [together830{nb_zone}'];
        zone.label{idzone} = [num2str(nb_zone),'_wl1'];
        zone.color= [zone.color; colorlist(nb_zone,:)] ;  
        %690
        idzone = idzone + 1;
        zone.plot{idzone} = [MeasList(together830{nb_zone},1),MeasList(together830{nb_zone},2)];
        zone.plotLst{idzone} = [together830{nb_zone}'+nbch];
        zone.label{idzone} = [num2str(nb_zone),'_wl2'];
        zone.color= [zone.color; colorlist(nb_zone,:)] ; 
        end
    end
    end
end

