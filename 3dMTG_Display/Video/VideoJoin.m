[name,path] = uigetfile('.avi')
file1 = [path,name];
[name,path] = uigetfile('.avi')
file2 = [path,name];
[name,path] = uiputfile('.avi')
fileout = [path,name];

aviobj = avifile(fileout,'fps',1) %fichier à creer
for i =1:25
    mov = aviread(file1,i); %initial à couper
    mov2.cdata = mov.cdata; 
    mov2.colormap = [];       
    mov = aviread(file2,i); %initial à couper
    mov2.cdata = cat(1,mov2.cdata, mov.cdata);
    %ADD time en haut
    h = figure    
    close(h)     
    aviobj = addframe(aviobj,mov2);  %ajouter la frame au film
end
aviobj = close(aviobj);
