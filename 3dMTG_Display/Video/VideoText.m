function VideoText(filename,newfilename,texttoadd,posx,posy,dframe)
%filename = 'Avg_MotScale.avi';
% texttoadd = ' Word reading';
% newfilename = 'Avg_MotScale_titre2.avi';

% filename = 'Avg_NonMotScale.avi';
% texttoadd = ' Non-Word reading';
% newfilename = 'Avg_NonMotScale_titre.avi';
% % 
% filename = 'videomot_nonmot_souspersubject.avi';
% texttoadd = 'Word vs. Non-Word reading';
% newfilename = 'videomot_nonmot_souspersubject_titre.avi';

h = figure    
set(gca,'visible','off')
text(0,0,[texttoadd],'FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(h,'temp','tif');
close(h)    
imgall = imread('temp','tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = indx(1):indx(end);
indy = indy(1):indy(end);
titre = imgall(indx,indy,:); 
[dy,dx,c] = size(titre);
try
aviobj = avifile(newfilename,'fps',1)
for i =dframe
    mov = aviread(filename,i); %initial 
    mov.colormap = [];        
    mov.cdata(posy:posy+dy-1, posx:posx+dx-1,:) =titre;
    aviobj = addframe(aviobj,mov);  %ajouter la frame au film
end
aviobj = close(aviobj);
catch
    'error text'
    aviobj = close(aviobj);
end