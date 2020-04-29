setappdata(0,'UseNativeSystemDialogs',false)
[name,path]=uigetfile('.vcolor','multiselect','on')
if ~iscell(name)
    name = {name}
end
all = [];
allnbpos = [];
allnbneg = [];
for iname=1:numel(name)
    load([path,name{iname}],'-mat');
   all = [all,vColor];
   vpos= [abs(vColor) > 7 ] ;
   vneg= -[vColor < -3 ] ;
  allnbpos = [allnbpos, vpos];
  allnbneg = [allnbneg, vneg];
end

vColor = sum(allnbpos')';
save([path,'s1_40absallnb7.vcolor'],'vColor','-mat');
vColor = sum(allnbneg')';
save([path,'negallnb.vcolor'],'vColor','-mat');

% % max(VertexBuffer(:,1))
% % max(VertexBuffer(:,2))
% % max(VertexBuffer(:,3))
% vColor = mean(all')'
% save([path,'1_40all.vcolor'],'vColor','-mat');
% 
% vColor = mean(all')'./std(all')';
% save([path,'Tall.vcolor'],'vColor','-mat');
