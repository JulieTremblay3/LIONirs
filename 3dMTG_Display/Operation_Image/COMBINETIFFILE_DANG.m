%IPSICONTRO
% [file,path]= uigetfile('*.tif')
% [img] = imread([path,file],'tif'); 
% Montage=imageresize(img,400,[]);
% imwrite(Montage,[path,file, 'all.tif'],'tif');



path = 'D:\Data\Epilepsie\Frontal(Reanalyse)\Article_janv2012\Figure\epi137\'
%Montage
sujet = 'epi137_Sz01_Bloc1_278'
delta = 0; 
%MONTAGE 
file = [sujet(1:6),'_Montage.TIF']
[img] = imread([path,file],'tif'); 
Montage=imageresize(img,300,[]);
spacel = ones(300,100,3)*255; 
spacer = ones(300,1300-size(Montage,2),3)*255; %
Montage = [spacel,Montage,spacer];
  

%IPSICONTRA
file = [sujet,'_ipsicontra.fig']
open([path,file])
set(gcf,'unit','normalized','position',[0,0,1,0.5])

lignes =get(gca,'children');
for i=1:numel(lignes)
    xdata=get(lignes(i),'xdata')
    xnew = (xdata+delta)
    set(lignes(i),'xdata',xnew)
end 
    
saveas(gcf,[path,file(1:end-3),'jpg'],'jpg')
[imgall] = imread([path,file(1:end-4),'_2','.jpg'],'jpg'); 

indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
img = imgall(indx,indy ,:); 
IPSI=imageresize(img,[],1150);
spacel = ones(size(IPSI,1),100,3)*255; 
spacer = ones(size(IPSI,1),150,3)*255; %
IPSI = [spacel,IPSI,spacer];


%EEG
file = [sujet,'_EEG.png']
[imgall] = imread([path,file],'png'); 
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);   
img = imgall(indx,indy ,:);  
EEG=imageresize(img,[],1200);
spacel = ones(size(EEG,1),100,3)*255; 
spacer = ones(size(EEG,1),100,3)*255; %
EEG = [spacel,EEG,spacer];
figure;image(imgall)
%TOPO
file = [sujet,'_Topo.tif']
[img] = imread([path,file],'tif'); 
TOPO=imageresize(img,[],1400);

spaceA = 255*ones(50,size(TOPO,2),3);
h = figure    
set(gca,'visible','off')
text(0,0,['A'],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
titre = imgall(indx,indy ,:); 
[dy,dx,c] = size(titre);
spaceA(1:dy,(1:dx)+40,:)=titre;

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



Total = [spaceA;Montage;spaceB;IPSI;spaceC;EEG;spaceD;TOPO];
imwrite(Total,[path,sujet, 'all.tif'],'tif');

close all

