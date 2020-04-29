setappdata(0,'UseNativeSystemDialogs',false)
[name,path]=uigetfile('.vcolor','multiselect','on')


if ~iscell(name)
    name = {name}
end
all = [];
allnbpos = [];
allvMax = [];
vmaxall = [];
for iname=1:numel(name)
    load([path,name{iname}],'-mat');
   all = [all,vColor];
   [vmax,vind]= max(abs(vColor)) ;
   vmaxall = [vmaxall,vmax]
   vColor(:)=0;
   vColor(vind)=1;
   allvMax = [allvMax,vColor];
end
vmaxall'
figure
plot([ vmaxall])

vColor = sum(allvMax')';
save([path,'posallvMax.vcolor'],'vColor','-mat');
% vColor = sum(allnbneg')';
% save([path,'negallnb.vcolor'],'vColor','-mat');

% max(VertexBuffer(:,1))
% % max(VertexBuffer(:,2))
% % max(VertexBuffer(:,3))
% save([path,'all.vcolor'],'vColor','-mat');
% 
% vColor = mean(all')'./std(all')';
% save([path,'Tall.vcolor'],'vColor','-mat');
