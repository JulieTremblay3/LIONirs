%Fluence Lateralité index par région
file_MaskL ={ 'F:\Fluence\MTG_MRI\NORS_NIRS\L_Broca.vcolor',
    'F:\Fluence\MTG_MRI\NORS_NIRS\L_Wernicke.vcolor'
    'F:\Fluence\MTG_MRI\NORS_NIRS\Hemi_left.vcolor'}
file_MaskR ={'F:\Fluence\MTG_MRI\NORS_NIRS\R_Broca.vcolor',
    'F:\Fluence\MTG_MRI\NORS_NIRS\R_Wernicke.vcolor',
    'F:\Fluence\MTG_MRI\NORS_NIRS\Hemi_right.vcolor'}
pathout = 'F:\Fluence\Review\LI_Time_Broca\'

path = 'F:\Fluence\Review\Video\'
name = {'01_fvl313'
    '02_FVE04'              %4
    '03_FVE30'               %4
    '04_FVE17'            %5
    '05_FVE34'               %5
    '06_FVE23'               %5
    '07_FVE28'               %6
    '08_FVE06'              %6
    '09_FVE14'            %6
    '10_FVE03'              %6
    '11_FVE01'              %7
    '12_FVE24'                 %7
    '13_FVE02'                %8
    '14_FVE33'                 %8
    '15_FVE26'                 %8
    '16_FVE07'                %9
    '17_FVE35'                 %9
    '18_FVE05'              %10
    '19_FVE32'                 %10
    '20_FVE36'                 %10
    '21_FVE22'                 %11
    '22_FVE15'               %11
    '23_FVE16'               %13
    '24_FVE13'               %13
    '25_FVE31'                  %14
    '26_FVE18'              %14
    '27_FVE09'                 %15
    '28_FVE11'               %15
    '29_FVE12'              %15
    '30_FVE08'                 %16
    '31_FVE27'                  %16
    '32_FVE10'               %A
    '33_FVE19'              %A
    '34_FVE25'                  %A
    '35_FVE21'                  %A
    '36_FVE20'                  %A
    '37_FVE29'                  %A
    '38_pFVE310'                 %A
    '39_FVE311'                 %A
    '40_FVE312'}                 %A
filevertex = 'F:\Fluence\MTG_MRI\NORS_NIRS\skinlow.SRX';
allIL_time = cell(3,1)
Thresholdpositif = ones(40,1); %Plus petit ou plus grand que le threshold
% Thresholdpositif(2)=-1
% Thresholdpositif(6)=-1
% Thresholdpositif(39)=-1
% Thresholdpositif(40)=-1

t = round([38.256 %tmax
11
27.1201
7.408
24.8161
23.5361
26.3521
5.36
26.8641
32.2401
12.2721
6.64
33.0081
5.872
28.1441
23.7921
14.8321
23.5361
29.1681
22.5121
18.1601
10.7361
16.1121
28.1441
13.5521
11.7601
13.5521
15.0881
20.9761
18.6721
7.408
21.4881
17.3921
15.6001
29.4241
28.9121
24.0481
11.888
26.36
28.92])


labelthreshold = '0.95' % pourcentage du T max
%t = 1:40;
for izone = 1:3
    load(file_MaskL{izone},'-mat');
    idLeft = find(vColor>0);
    load(file_MaskR{izone},'-mat');
    idRight = find(vColor >0);
    for iname = 1:numel(name)
        hfigure = figure;
        id_time = 1;
        for itime = t(iname)
            try
            load([path,name{iname},'.aviHbOTime   ',num2str(itime),'(s).vcolor'],'-mat');

            if Thresholdpositif(iname) == 1
                threshold = max(vColor)*str2num(labelthreshold);
            else
                threshold = min(vColor)*str2num(labelthreshold);
            end
            vleft = vColor;
            vleft(idRight)=0;
            vright= vColor;
            vright(idLeft) =0;
            clear IL
            for i=1:numel(threshold)
                if max(vColor)< 1.812461
                    IL(i) = NaN
                else
                    if threshold(i)> 0
                        IL(i) = (numel(find(vleft>threshold(i)))-numel(find(vright>threshold(i))))/(numel(find(vleft>threshold(i)))+numel(find(vright>threshold(i))));
                    else
                        IL(i) = (numel(find(vleft<threshold(i)))-numel(find(vright<threshold(i))))/(numel(find(vleft<threshold(i)))+numel(find(vright<threshold(i))));
                    end
                end
            end
            IL_time(id_time) = IL;
            threshold_time(id_time) = threshold;
            id_time = id_time+1;
            catch
                IL_time(id_time) = 0;               
                threshold_time(id_time) =  0;
                id_time = id_time+1;
            end
        end
        subplot(2,1,1)
        hold on
        plot(t,IL_time)
        ylim([-1,1])
        xlabel('Time (s)')
        ylabel('IL')
        title(['IL,',labelthreshold ,' ', name{iname},'HbO'])
        allIL_time{izone}  = [allIL_time{izone} ;IL_time]
        plot([t(1),t(end)],[0,0],'k')
        subplot(2,1,2)
        plot(t,threshold_time)
        xlabel('Time (s)')
        ylabel('T value')
        saveas(hfigure,[pathout,num2str(izone),labelthreshold ,name{iname},'.jpg'],'jpg')
        close(hfigure)
        
    end
end 
% figure
% hold on
% plot(mean(allIL_time(1:10,:),1),'r')
% plot(mean(allIL_time(11:20,:),1),'b')
% plot(mean(allIL_time(21:30,:),1),'g')
% plot(mean(allIL_time(31:40,:),1),'k')
% figure
% imagesc(allIL_time)
% caxis([-1,1])
% imagesc(mean(allIL_time'))
% caxis([-1,1])
colorjet = colormap(jet(101));

figure
subplot(1,3,1)
imagesc(-allIL_time{1})
title('Broca')
xlabel('Time')
ylabel('Subject')
subplot(1,3,2)
imagesc(-allIL_time{2})
title('Wernicke')
xlabel('Time')
ylabel('Subject')
subplot(1,3,3)
imagesc(-allIL_time{3})
title('Hemi')
xlabel('Time')
ylabel('Subject')