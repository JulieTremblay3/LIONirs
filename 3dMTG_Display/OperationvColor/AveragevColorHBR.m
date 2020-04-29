[name,path]=uigetfile('.vcolor','multiselect','on')
h=figure
hold on
if ~iscell(name)
    name = {name}
end
all = [];
idLeft= find(VertexBuffer(:,1)>181/2);
idRight = find(VertexBuffer(:,1)<181/2);
for iname=1:numel(name)
    load([path,name{iname}],'-mat');
        plot([-10,0],[0,0],'k')

    vleft = vColor;
    vleft(idRight)=0;
    vright= vColor;
    vright(idLeft) =0;
  %  v= [vColor > 2 ] ;
  % all = [all,v]
   %  all = [all,vColor]
   clear IL
   clear  threshold
   threshold = -2:-1:-6
    for i=1:numel(threshold)
        IL(i) = (numel(find(vleft<threshold(i)))-numel(find(vright<threshold(i))))/(numel(find(vleft<threshold(i)))+numel(find(vright<threshold(i))))
    end
   
    h=plot(threshold,IL,'linewidth',3)
    sujet = name{iname}
    set(h,'displayname',sujet(1:end-4))
        ylim([-1.2,1.2])
%     title(sujet)
%     saveas(h,[path,sujet(1:end-4),'.jpg'],'jpg')
    save('all.vcolor','vColor','-mat')   
end

%vColor = sum(all')';
%vColor = mean(all')';
% max(VertexBuffer(:,1))
% max(VertexBuffer(:,2))
% max(VertexBuffer(:,3))

save('all.vcolor','vColor','-mat')