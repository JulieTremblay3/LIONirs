function VideoColorBar(filename,newfilename,filelegende,dframe)
 
imgall = imread(filelegende,'tif');     
[dy,dx,c] = size(imgall ); % attention il faut que le .tif soit plus petit que le film
aviobj = avifile(newfilename,'fps',1)
try
for i =dframe
    mov = aviread(filename,i); %initial  
    [y,x,c]=size(mov.cdata);
    img=imageresize(imgall,y-100,[]); 
    [dy,dx,c] = size(img);
    posy = (y - dy)/2;
    posx = x;
    mov.cdata = cat(2,mov.cdata, ones(y,dx,3)*255);
    posy = (y - dy)/2;
    posx = x;
    mov.cdata(posy+1:posy+dy, posx+1:posx+dx,:) =img;
    aviobj = addframe(aviobj,mov);  %ajouter la frame au film
end
    aviobj = close(aviobj);
catch    
    'error ColorBar'
    aviobj = close(aviobj);
end
