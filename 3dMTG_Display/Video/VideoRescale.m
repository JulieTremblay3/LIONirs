function VideoRescale(filein,fileout,dx,dy,dframe,time)
aviobj = avifile(fileout,'fps',1) %fichier à creer
try
for i=dframe
    mov = aviread(filein,i); %initial à couper  
    mov2.cdata = mov.cdata(dx,dy,:); 
%     mov2.cdata(70:160,1130:1176,:)=255;   
    mov2.colormap = [];        
    %ADD time en haut
    h = figure    
    labelt=sprintf('%01.2f',time(i))
    text(0.5,0.5,['Time : ', labelt,'s'],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
    saveas(h,'temp','tif');
    close(h)    
    imgall = imread('temp','tif');  
    % figure; imagesc(imgall)
    nb = imgall(330:450,600:900,:); 
    [nbx,nby,x] = size(nb);
    mov2.cdata (1:nbx,1:nby,:) = nb;   
    aviobj = addframe(aviobj,mov2);  %ajouter la frame au film
end
aviobj = close(aviobj);
catch
    'error Rescale'
    aviobj = close(aviobj);
end