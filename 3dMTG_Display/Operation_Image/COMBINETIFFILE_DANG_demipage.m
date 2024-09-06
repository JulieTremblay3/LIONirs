%IPSICONTRO
% [file,path]= uigetfile('*.tif')
% [img] = imread([path,file],'tif'); 
% Montage=imageresize(img,400,[]);
% imwrite(Montage,[path,file, 'all.tif'],'tif');

%NOTE le ipsi contro doit être sauvegarder manuellement


path = 'E:\Data\Epilepsie\Frontal(Reanalyse)\Article_janv2012\Figure\epi103\'
%Montage
sujet = 'epi103_Sz01_Bloc1_205'
delta = 0; 




%IPSICONTRA
file = [sujet,'_ipsicontra.fig']
open([path,file])
set(gcf,'unit','normalized','position',[0,0,1,0.5])
% set(gca,'unit','normalized','position',[0,0,1,0.5]) 
lignes =get(gca,'children');

xaxis = get(gca,'xlim')
set(gca,'xlim',xaxis+delta)
for i=1:numel(lignes)
    xdata=get(lignes(i),'xdata');  
    xnew = (xdata+delta);
    set(lignes(i),'xdata',xnew);
end 
    
%saveas(gcf,[path,file(1:end-3),'tif'],'tif') 
%A faire manuellement  


[imgall] = imread([path,file(1:end-4),'_2.jpg'],'jpg'); 

indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
img = imgall(indx,indy ,:); 
IPSI=imageresize(img,[],1150);
spacel = ones(size(IPSI,1),100,3)*255; 
spacer = ones(size(IPSI,1),150,3)*255; %
IPSI = [spacel,IPSI,spacer];


%TOPO
file = [sujet,'_Topo.tif']
[img] = imread([path,file],'tif'); 
TOPO=imageresize(img,[],1400);



spaceB = 255*ones(50,size(TOPO,2),3);
h = figure    
set(gca,'visible','off')
text(0,0,['B'],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
titre = imgall(indx,indy,:); 
[dy,dx,c] = size(titre);
spaceB(1:dy,(1:dx)+40,:)=titre;

spaceC = 255*ones(50,size(TOPO,2),3);
h = figure    
set(gca,'visible','off')
text(0,0,['C'],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
titre = imgall(indx,indy ,:); 
[dy,dx,c] = size(titre);
spaceC(1:dy,(1:dx)+40,:)=titre;

spaceD = 255*ones(50,size(TOPO,2),3);
h = figure    
set(gca,'visible','off')
text(0,0,['D'],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
titre = imgall(indx,indy,:); 
[dy,dx,c] = size(titre);
spaceD(1:dy,(1:dx)+40,:)=titre;


  

Total = [spaceB;IPSI;spaceC;TOPO];

%MONTAGE ajouter à gauche
file = [sujet(1:6),'_Montage.TIF']
[img] = imread([path,file],'tif'); 
Montage = imageresize(img,300,[]);




h = figure  
spaceA = 255*ones(50,size(Montage,2),3);
set(gca,'visible','off')
text(0,0,['A'],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
titre = imgall(indx,indy,:); 
[dy,dx,c] = size(titre);
spaceA(1:dy,(1:dx)+40,:)=titre;

Montage = [spaceA;Montage];
imwrite(Montage,[path,sujet, 'alldemi.tif'],'tif');

spacebottom = ones(size(Total,1)-size(Montage,1),size(Montage,2),3)*255; %
Montage = [Montage;spacebottom];
Total = [Montage,Total];

imwrite(Total,[path,sujet, 'alldemi.tif'],'tif');

close all

