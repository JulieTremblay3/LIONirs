path = {'F:\Fluence\Review\MaxT_time\G1\vcolor\'
     'F:\Fluence\Review\MaxT_time\G2\vcolor\'
     'F:\Fluence\Review\MaxT_time\G3\vcolor\'
	'F:\Fluence\Review\MaxT_time\G4\vcolor\'}
file_MaskL ={ 'F:\Fluence\MTG_MRI\NORS_NIRS\L_Broca.vcolor',
    'F:\Fluence\MTG_MRI\NORS_NIRS\L_Wernicke.vcolor'
    'F:\Fluence\MTG_MRI\NORS_NIRS\Hemi_left.vcolor'}
file_MaskR ={'F:\Fluence\MTG_MRI\NORS_NIRS\R_Broca.vcolor',
    'F:\Fluence\MTG_MRI\NORS_NIRS\R_Wernicke.vcolor',
    'F:\Fluence\MTG_MRI\NORS_NIRS\Hemi_right.vcolor'}
ratiothreshold = ' 0.80 '
ILallzone = cell(3,1);
for izone = 1:3   
    ILall = [];
for ipath=1:numel(path)
file = dir(path{ipath})
    for i = 1:numel(file)
    [pathstr, name, ext] = fileparts(file(i).name)
    if strcmp(ext,'.vcolor')  
            IL =  []
            load(file_MaskL{izone},'-mat');
            idLeft = find(vColor>0);
            load(file_MaskR{izone},'-mat');
            idRight = find(vColor >0);
            load([path{ipath}, name, ext],'-mat');
            vleft = vColor;
            vleft(idRight)=0;
            vright= vColor;
            vright(idLeft) =0;
            threshold = max(vColor)*str2num(ratiothreshold);
            for itresh=1:numel(threshold)
                if threshold(itresh)> 0
                    IL(itresh) = (numel(find(vleft>threshold(itresh)))-numel(find(vright>threshold(itresh))))/(numel(find(vleft>threshold(itresh)))+numel(find(vright>threshold(itresh))));
                else
                    IL(itresh) = (numel(find(vleft<threshold(itresh)))-numel(find(vright<threshold(itresh))))/(numel(find(vleft<threshold(itresh)))+numel(find(vright<threshold(itresh))));
                end
            end
            ILall = [ILall;IL]
    end            
    end
end
ILallzone{izone}= ILall;
end

colorset = jet(200);
for i=1:80;colorset(i,:)= [0;0;0.5];end
for i=81:120;colorset(i,:)= [0;0.5;0];end
for i=121:200;colorset(i,:)= [0.5;0;0];end

figure
colormap(colorset)
subplot(1,3,1)
imagesc(-ILallzone{1})
title(['Broca',ratiothreshold])
xlabel('Time')
ylabel('Subject')
caxis([-1,1])
subplot(1,3,2)
imagesc(-ILallzone{2})
title(['Wernicke',ratiothreshold])
xlabel('Time')
ylabel('Subject')
caxis([-1,1])
subplot(1,3,3)
imagesc(-ILallzone{3})
title(['Hemi',ratiothreshold])
xlabel('Time')
ylabel('Subject')
caxis([-1,1])



figure
x = ILallzone{1};
subplot(1,3,1)
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
title(['Broca',ratiothreshold])
subplot(1,3,2)
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
x = ILallzone{2};
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
title(['Wernicke',ratiothreshold])
subplot(1,3,3)
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
x = ILallzone{3};
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
title(['Hemi',ratiothreshold])

%Enfant
figure
x = ILallzone{1};
x = x(1:10)
subplot(4,3,1)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
title(['G1 Broca',ratiothreshold])
subplot(4,3,2)
x = ILallzone{2};
x = x(1:10)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
title(['G1 Wernicke' ,ratiothreshold])
subplot(4,3,3)
x = ILallzone{3};
x = x(1:10)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
title(['G1 Hemi',ratiothreshold])
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))

x = ILallzone{1};
x = x(11:20)
subplot(4,3,4)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2));
title('G2 Broca')
subplot(4,3,5)
x = ILallzone{2};
x = x(11:20)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2));
title('G2 Wernicke')
subplot(4,3,6)
x = ILallzone{3};
x = x(11:20)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
title('G2 Hemi')

x = ILallzone{1};
x = x(21:30)
subplot(4,3,7)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
title('G3 Broca')
subplot(4,3,8)
x = ILallzone{2};
x = x(21:30)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
title('G3 Wernicke')
subplot(4,3,9)
x = ILallzone{3};
x = x(21:30)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
caxis([1,2])
title('G3 Hemi')

x = ILallzone{1};
x = x(31:40)
subplot(4,3,10)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
title('G4 Broca')
subplot(4,3,11)
x = ILallzone{2};
x = x(31:40)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
title('G4 Wernicke')
subplot(4,3,12)
x = ILallzone{3};
x = x(31:40)
pie([sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2)]);
correctcaxis(sum(x>0.2),sum(x<=0.2&x>=-0.2),sum(x<-0.2))
title('G4 Hemi')
