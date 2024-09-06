% Utiliser le fichier de setting.m qui contient les blocs à ouvrir
% en format .nir (de spm nirs) raw on fait une normalisation diviser par la
% on ajoute un filtre passe bas et on normalise par la mediane du bloc
% ou en format .nirs de homer
% pathname = 'C:\data\LFC01_SM\pc00-00\';
% filename = [{'LFC101s1-DC.nirs'}];%
% ecrire les temps à combiner de même durée dans la variable
% dt = 5 (+ ou - 5 seconde autour du spike
% %timingfile{ifilename} = [44]
% timingfile{1} = [timingfile{1}-dt,timingfile{1}+dt];
% le but est d,utiliser par la suite le programme eegpermutation
% pour enlever les blocs outliers et de faire les permutations pour voir si
% la courbe du moyennage donne une composante significative
%******************** Changer ici
filesetting = 'C:\data\LFC01_SM\Trig_Spike_pc00_interieurbloclecture.m'
%filesetting = 'C:\data\LFC01_SM\Trig_Spike_small_nir.m'
fileele = 'C:\data\LFC01_SM\allele.ele'
fileout = 'spike2.dat'
pathout = 'C:\data\LFC01_SM\Analyse01_pc00interieurbloclecture\'
 
baseline = 1:138;   %-10 à 3 pour la comparaison T test
baseline = 138:197; %-3 à 0
baseline = 1:197;   %-10 à 0
clear baseline
baseline = 1:138;   %-10 à 3
peaklist =138:5:392; %en indice les temps pour sauvegarder une topo




%********************
lpf =0.5; % filtre passe bas
run(filesetting)% Contient  - filename fichier à ouvrir en .nir ou .nirs)
%- timingfile (une cellule par nom de fichier et les temps à regarder
fixsize  =[];
data = [];

for ifilename = 1:numel(filename)
    file = [pathname,filename{ifilename}];
    [pathstr, name, ext] = fileparts(file)
    datatime =timingfile{ifilename}; %trig
    
    if strcmp(ext,'.nir') % nir file spm
        filevhdr = [file(1:end-3),'vhdr']
        info = read_vhdr(filevhdr)
        d = fopen_NIR(file,info.NumberOfChannels);
        samplingintervals = info.SamplingInterval/1000000;
        pAstart = 1; % round(Astart/samplingintervals);
        pAstop =size(d,2);%round(Astop/samplingintervals);
        t = pAstart*samplingintervals:samplingintervals:pAstop*samplingintervals;
    else
        load([pathname,filename{1}],'-mat') %PC00-00 contains Normalised unfiltered pulse corrected data and corresponding averaged data
        d = d';
        samplingintervals  =1/19.5312;
    end
    
    
    %petit filtre pour nir spm
    if strcmp(ext,'.nir')
        AdvOptions.FiltOptions.FilterType=1;  %Default Butterworth
        AdvOptions.FiltOptions.FilterOrder=3;  %Default 3rd order
        fs = 1/(t(2)-t(1));
        [fb,fa] = MakeFilter(AdvOptions.FiltOptions.FilterType,AdvOptions.FiltOptions.FilterOrder,fs,lpf,'low');
        d= filtfilt(fb,fa,d');
        d = d';
    end
    
    for istim=1:size(datatime,1)
        %Segment par
        t1=datatime(istim,1); %temps debut en seconde
        t2=datatime(istim,2);%temps debut en seconde
        t1i = find(t<=t1);
        t1i=t1i(end);
        t2i = find(t<=t2);
        t2i=t2i(end);
        indt = t1i:t2i;
        
        if isempty(fixsize)
            if mod(numel(indt),2)
                indt = indt(1):(indt(1)+numel(indt));
            end
            fixsize = numel(indt);
            tfinal = -fixsize/2*samplingintervals:samplingintervals:samplingintervals*(fixsize/2-1);
        else
            indt = indt(1):(indt(1)+fixsize-1);
        end
        clear intensnorm
        if  strcmp(ext,'.nir')
            meanIntens = median(d(:,indt),2);
            tmp = d';
            tmp = tmp(indt,:)./(ones(numel(indt),1)*meanIntens'); %intensite normalise non filtré ; %intensite normalise non filtré
            intensnorm = tmp';
            [m,p]=size(intensnorm);
            tmp = reshape(intensnorm,m,1,p);
        else
            intensnorm = d(:,indt);
            [m,p]=size(intensnorm);
            tmp = reshape(intensnorm,m,1,p);
        end
        data = cat(2,data,tmp);
        %electrode, trial, temps
    end
end

figure
hold on 

indgood = find(std(mean(abs(data(:,:,:)),2),1,3)< 0.01);
A = zeros(size(data,1));
A(indgood)=1
save([pathout,'measlistact.mat'],'A','-mat')

colorlst = lines(numel(indgood))
%%%%%%TEST JULIE retranche le baseline à chaque essais
% ibaseline = zeros(size(data));
% for iele = 1:size(data,1)
%     for itrial = 1:size(data,2)
%         data(iele,itrial,:) = data(iele,itrial,:) - mean(data(iele,itrial,baseline),3);
%         ibaseline(iele,itrial,:) = mean(data(iele,itrial,baseline),3);
%     end
% end
%%%%%%%%%%%%%
%figure;plot(squeeze(d(1,:)'));


%AFFICHAGE dOD moyennée
figure
hold on
offset = 0;
for i = 1:numel(indgood)
    if 1
        d1 = squeeze(mean(data(indgood(i),:,:),2))'+offset;
    else
        d1 = squeeze(data(indgood(i),1,:))'+offset;
    end
    plot(tfinal,d1,'displayname',num2str(i),'color',colorlst(i,:))
    if ~mod(i,30)
        offset = offset + 0.2;
    end
end
ch1 = data(1,:,:);

%AFFICHAGE T et sauvegarde de chaque topo en .mat
figure
hold on
offset = 0;
doffset = 10;
Tval = zeros(size(data,1),1);

idlist = 1;
for peak = peaklist
    for i = 1:numel(indgood)
        d1 = (squeeze(mean(data(indgood(i),:,:),2))'-mean(squeeze(mean(data(indgood(i),:,baseline),2))'))./(squeeze(std(data(indgood(i),:,:),1,2))'./sqrt(size(data,2)));
        if idlist == 1
            plot(tfinal,d1,'displayname',num2str(i),'color',colorlst(i,:))
        end
        if ~mod(i,30)
            offset = offset + doffset;
        end
        Tval(indgood(i)) = d1(peak);
    end
    
    % figure
    % plot(tfinal)
    %Save les topos !
     
    labelpeak = sprintf('%01.2f',tfinal(peak));
    labelidlist = sprintf('%03.0f',idlist)
    label_bstart = sprintf('%01.2f',tfinal(baseline(1)));
    label_bstop =sprintf('%01.2f',tfinal(baseline(end)));
    title(['Baseline Time : ',label_bstart,'to',label_bstop])
    A = Tval(1:end/2);
    save([pathout,'HbO_peak',labelidlist,'_',labelpeak,'_',label_bstart,'to',label_bstop,'.mat'],'A','-mat');
    A = Tval(end/2+1:end);
    save([pathout,'HbR_peak',labelidlist,'_',labelpeak,'_',label_bstart,'to',label_bstop,'.mat'],'A','-mat');
    
    idlist = 1 + idlist
end
time  = tfinal(peaklist);
save([pathout,'time.mat'],'time','-mat')
%Cross Correlation !
figure
hold on
offset = 0;
doffset = 10;
crossvalmax=zeros(size(data,1),1);
paradigm = zeros(size(data,2),size(data,3));
timeparadigm = 197:197+5
paradigm(:,timeparadigm)=1;
paradigm = paradigm';
for i = 1:numel(indgood)
    d1 = squeeze(data(indgood(i),:,:))';
    c = xcorr(d1(:),paradigm(:));
    %plot(c(:),'displayname',num2str(i),'color',colorlst(i,:));
    if ~mod(i,30)
        offset = offset + doffset;
    end
    crossvalmax(indgood(i)) = max(abs(c));
end

label_bstart = sprintf('%01.2f',tfinal(timeparadigm(1)));
label_bstop =sprintf('%01.2f',tfinal(timeparadigm(end)));
A =  crossvalmax(1:end/2);
save([pathout,'CrosscorrelationHbO',label_bstart,'to',label_bstop,'.mat'],'A','-mat');
A =  crossvalmax(end/2+1:end);
save([pathout,'CrosscorrelationHbR',label_bstart,'to',label_bstop,'.mat'],'A','-mat');



%
% A =1-FUniv;
% figure;plot(A)
% save('Test.mat','A','-mat')
% write_dat(data-1,tfinal,fileele,fileout)

