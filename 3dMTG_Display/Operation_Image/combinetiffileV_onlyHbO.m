function [pathout,nameout] = combinetiffileV(path,file,cmin,cmax,cthresh)
Total = [];
for i = 1:numel(file)
    [X] = imread([path,file{i}],'tif'); 
    Total = [Total;X];
end
figure
imagesc(Total)
[nbxT,nby,z]= size(Total)

addmarge = ones(nbxT,840-420,3,'uint8')*255;
h=figure
text(0.5,0.5,['HbO'],'FontSize',50,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');     
nb = imgall(330:420,610:840,:); % attention y,x sur la figure pour couper le chiffre
[nbx,nby,x] = size(nb);
addmarge(1:nbx,1:nby,:)= nb;


% h=figure
% text(0.5,0.5,['HbR'],'FontSize',50,'HorizontalAlignment','left','VerticalAlignment','bottom')
% saveas(h,'temp','tif');
% close(h)    
% imgall = imread('temp','tif');     
% nb = imgall(330:420,610:840,:); % attention y,x sur la figure pour couper le chiffre
% [nbx,nby,x] = size(nb);
% addmarge(round(nbxT/2)+(1:nbx),(1:nby),:)= nb;



h=figure
haxe = axes
set(haxe,'xtick',[])
map = jet(50);
nb_color = size(map,1);
delta_map = (cmax-cmin)/nb_color;
nb_map_min = floor(abs(cmin+cthresh)/delta_map);  
nb_map_max = ceil(abs(cmin-cthresh)/delta_map);
nb_map_list = nb_map_min:nb_map_max;
nb_map = size(nb_map_list,2);
for i_map = 1:nb_map;
map(nb_map_list(i_map),:) = [255/256 255/256 179/256];
end
set(h,'colormap',map)
set(haxe,'clim',[cmin,cmax])
hcol = colorbar
set(hcol,'fontsize',20)
pos = get(hcol,'position')
set(hcol,'position',[pos(1),pos(2),0.1,pos(4)])
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');   
figure
imagesc(imgall)
nb = imgall(1:850,980:1200,:); % attention y,x sur la figure pour couper le chiffre
[nbx,nby,x] = size(nb);
addmarge((round(nbxT/2)-size(nb,1)/2):(round(nbxT/2)+size(nb,1)/2-1),round(nby/2)+(1:nby),:)= nb;

Total = [Total,addmarge];
figure
imagesc(nb)
[nameout,pathout]=uiputfile('all.tif');
imwrite(Total,[pathout,nameout],'tif');

