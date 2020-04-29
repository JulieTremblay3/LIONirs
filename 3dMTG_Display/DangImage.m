
angle_face = [40,40,40,40]
angle_cote = [120,-120,120,-120]
conc = [1,1,2,2] 
time = [1,2,3,4,5,6]
outputname = 'epi103_Sz02_Bloc1_321'
% path = 'D:\Data\Epilesie\Frontal(Reanalyse)\epi109OC\'
% angle_face = [10,10,10,10]
% angle_cote = [0,0,10,10]
% conc = [1,2,1,2]
% outputname = 'epi102_Sz01_Bloc6_401'
% time = [2,5,7,10,12,14]
 
% outputname = 'epi141_Sz20_Bloc2_780_Topo'
%  angle_cote = [180,90,-90,180,90,-90]; %137
%  angle_face = [35,35,35,35,35,35];
%  conc = [1,1,1,2,2,2];
 
%  %DIMA TEST AVEC VCOLOR FILE
%  outputname = 'MOT'
%  angle_cote = [90,-90,180]; %137
%  angle_face = [35,35,35];


%  conc = [1,1,1];
%  
%  timefile = {'HbTT2.179val000to4(s).vColor',
%      'HbTT2.179val005to8(s).vColor',
%      'HbTT2.179val009to12(s).vColor',
%      'HbTT2.179val003to16(s).vColor',
%      'HbTT2.179val017to20(s).vColor'}
 

%  angle_cote = [0,90,-90,0,90,-90]

% angle_cote = [180,180];
% angle_face = [40,40];
% conc = [1,2]
% angle_cote = [90,-90,90,-90]
% angle_face = [30,30,30,30]
% conc = [1,1,2,2]
% % %Marker de tempsiom
% time = [2,5,6,10,20,50];

% time = [2,5,10,15,20,25]
% %time = [8,10,16,25]
% time = [70,74,78]
%time = [20,40,60,80]
% time = [1,3,5,7,15];
% time = [6,10,22,26]
% time = [2,5]
%time= [3,10,15,20]
plot([time],-2*ones(numel(time),1),'marker','+','markersize',6,'color','k','linewidth',2,'linestyle','none','markerfacecolor','k')
% % % %Ligne EEG crise
% plot([0 4],[0,0],'linewidth',4,'color','k')
%  mark = [3.731541
% 7.242654
% 8.262975
% 9.563568
% 10.564997
% 11.662373
% 12.708943
% 13.958732
% 15.542135
% 16.852889]
%  plot([mark],zeros(numel(mark),1),'linewidth',4,'linestyle','none','marker','+','markersize',12,'color','k')


TotalV = [];
tic
indytext = 0;
p = mfilename('fullpath');
hfigure = figure;
hfigurestd = figure;
[pathstr, name, ext, versn]=fileparts(p);
for iplan = 1:numel(angle_cote)
    TotalH = [];
    cadragemaxx = 500; %attention mettre un nombre pair
    cadragemaxy = 500;
    couperbas = 130
    for itime=1:numel(time) 
        figure(hfigure)    
        topo_Dang(time(itime),angle_cote(iplan),angle_face(iplan),conc(iplan),hfigure);
        imgall = imread([pathstr,'\temp.tif'],'tif');  
        
        indy = find(sum(sum(imgall,3)~=765,1));
        indx = find(sum(sum(imgall,3)~=765,2));
        dx = indx(end)-indx(1);
        dy = indy(end)-indy(1);
        if dx <= cadragemaxx                      
            if ~mod(dx,2)%even 
               addx=(cadragemaxx-dx)/2;
               indx = (indx(1)-addx+1):(indx(end)+addx);
            else %odd
               addx=(cadragemaxx-dx+1)/2;
               indx = (indx(1)-addx+1):(indx(end)+addx-1);
            end            
        end
        if dy <= cadragemaxy                      
            if ~mod(dy,2) %odd
               addy=(cadragemaxy-dy)/2;
               indy = (indy(1)-addy+1):(indy(end)+addy);
            else %even
               addy=(cadragemaxy-dy+1)/2;
               indy = (indy(1)-addy+1):(indy(end)+addy-1);
            end            
        end       
        
        %Couper le cadrage du bas
        indx = indx(1):(indx(end)-couperbas);
        
%         if iplan==1|iplan==3 %gauche epi114
%             indx = (indx(1)-20):(indx(end)-60);
%             indy = (indy(1)-20):(indy(end)-40); 
%         else %droite
%             indx = (indx(1)-20):(indx(end)-60);
%             indy = (indy(1)+40):(indy(end)+20);
%         end
%         if iplan==1|iplan==3 %gauche epi114
%             indx = (indx(1)-20):(indx(end)-10);
%             indy = (indy(1)-20):(indy(end)+20); 
%         else %droite
%             indx = (indx(1)-20):(indx(end)-10);
%             indy = (indy(1)-20):(indy(end)+20);
%         end

        clear img
        img = imgall(indx,indy,:);    
        if iplan == 1 %Temps marqués sur la première colonne     
            figure(hfigurestd)
            cla
            set(gca,'visible','off')
            text(0,0,[num2str(time(itime)),' s'],'FontSize',30,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','arial')
            saveas(hfigurestd,[pathstr,'\temp.tif'],'tif');                   
            clear imgall titre
            imgall = imread([pathstr,'\temp.tif'],'tif'); 
            if indytext == 0 %S'Assurer d'avoir toujours le même découpage en hauteur en initialisant avec le premier
                indytext = find(sum(sum(imgall,3)~=765,2));              
            end
            indx = find(sum(sum(imgall,3)~=765,1));            
            titre = imgall(indytext,indx,:);            
            [dy,dx,c] = size(titre);
            posx = 50;
            posy = 1;
            Margeup = ones(dy,size(img,2),3)*255;
            Margeup(posy:posy+dy-1, posx:posx+dx-1,:)=titre;       
            img = [Margeup;img];
        end
          TotalH = [TotalH,img]; %concatenation horizontal             
    end  
    if size(TotalH,2)< size(TotalV,2)
        addy = size(TotalV,2)-size(TotalH,2);
        TotalH = [TotalH,ones(size(TotalH,1),addy,3)*255];
    elseif size(TotalH,2)> size(TotalV,2)& ~isempty(TotalV)
        TotalH(:,size(TotalV,2)+1:end,:)=[];
    end
    TotalV = [TotalV;TotalH];
end
%HbO

figure(hfigurestd)
cla
posx = 1;
posy = 1;
set(gca,'visible','off')
text(0.5,0.5,['HbO'],'FontSize',50,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(hfigurestd,[pathstr,'\temp.tif'],'tif');   
imgall = imread([pathstr,'\temp.tif'],'tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = (indx(1)-20):(indx(end)+20);
indy = (indy(1)-20):(indy(end)+20);
titre = imgall(indx,indy,:); 
[dy,dx,c] = size(titre);
%ajouter une marge
MargeLeft = ones(size(TotalV,1),dx,3)*255;
MargeLeft(posy:posy+dy-1, posx:posx+dx-1,:)=titre;

%HbR
cla
posx = 1;
[y,x,c]=size(TotalV);
posy = y/2;
set(gca,'visible','off')
text(0.5,0.5,['HbR'],'FontSize',50,'HorizontalAlignment','left','VerticalAlignment','bottom')
saveas(hfigurestd,[pathstr,'\temp.tif'],'tif'); 
imgall = imread([pathstr,'\temp.tif'],'tif');     
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = (indx(1)-20):(indx(end)+20);
indy = (indy(1)-20):(indy(end)+20);
titre = imgall(indx,indy,:); 
[dy,dx,c] = size(titre);
MargeLeft(posy:posy+dy-1, posx:posx+dx-1,:)=titre;
TotalV = [MargeLeft,TotalV];

%color bar
cla
guiHelmet = getappdata(0,'guiHelmet');
handles = guihandles(guiHelmet);
set(gca,'visible','off')
Colormapfig = get(handles.IO_HelmetMTG,'colormap');
set(hfigurestd,'colormap',Colormapfig);
prop = get(handles.axes_Mtg2,'CLim');
set(gca,'CLim',prop);
clear prop
hcolorbar = colorbar;
set(hcolorbar,'fontsize',20)
saveas(hfigurestd,[pathstr,'\temp.tif'],'tif');
imgall = imread([pathstr,'\temp.tif'],'tif');
imgall = imresize(imgall,1.5);
indy = find(sum(sum(imgall,3)~=765,1));
indx = find(sum(sum(imgall,3)~=765,2));
indx = ((indx(1)):(indx(end)));
indy = ((indy(1)):(indy(end)));
titre = imgall(indx,indy,:); 
titre = [titre;ones(30,size(titre,2),3)*255];
titre = [ones(30,size(titre,2),3)*255;titre];
titre = [titre,ones(size(titre,1),30,3)*255];
titre = [ones(size(titre,1),30,3)*255,titre];
[dy,dx,c] = size(titre);
[y,x,c]=size(TotalV);
posx = 1;
posy = (y - dy)/2;
MargeRight = ones(y,dx,3)*255;
MargeRight(posy+1:posy+dy,posx+1:posx+dx,:)=titre;
figure;image(MargeRight/255)
TotalV = [TotalV,MargeRight];
toc
close(hfigurestd)
close(hfigure)
t1 = sprintf('%02.0f',time(1))
t2 = sprintf('%02.0f',time(end))
[nameout,pathout]=uiputfile([outputname,'t',t1,'to',t2,'.tif']);
imwrite(TotalV,[pathout,nameout],'tif');

figure
image(titre)