%Exemple avec les fichier .nirs (donnée non traité de homer)
%PCA 
filesetting = 'G:\Artefact_Nirsspm\Tono03_JP\checknoise_tono3_sondroit.m'
run(filesetting)

load([pathname,filename{1}],'-mat')
%initialise output 
intensnormall=zeros(size(d));
noise = zeros(size(d));
datatime =timingfile{1};

for istim=1:1
%Segment par 
t1=datatime(istim,1); %temps debut en seconde
t2=datatime(istim,2);%temps debut en seconde
t1i = find(t<=t1);
t1i=t1i(end);
t2i = find(t<=t2);
t2i=t2i(end);
indt = t1i:t2i;
figure;plot(t,d)
%filtre 
AdvOptions.FiltOptions.FilterType=1;  %Default Butterworth
AdvOptions.FiltOptions.FilterOrder=3;  %Default 3rd order
lpf = 0.2;
fs = 1/(t(2)-t(1)); 
[fb,fa] = MakeFilter(AdvOptions.FiltOptions.FilterType,AdvOptions.FiltOptions.FilterOrder,fs,lpf,'low');
d = filtfilt(fb,fa,d);    
meanIntens = median(d(indt,:));
%Normalisation 
intensnorm = d(indt,:)./(ones(numel(indt),1)*meanIntens); %intensite normalise non filtré 
measlistact = ones(size(intensnorm,2),1);
indbad =find(mean(d)<100);
measlistact(indbad) = 0; 
%Zonem correlation inter channel 0.8 
zone = correlation_channel_gen(intensnorm,measlistact); %ok
save('corr.zone','zone','-mat')
%pour l'exemple on enlève les zones à 1 canal
for i=1:numel(zone.plotLst)
    if numel(zone.plotLst{i})==1
        measlistact( zone.plotLst{i})=0;
    end
end
%PCA
A = measlistact
save('measlistact.mat','A','-mat')
%measlistact good or bad bloc if is necessary
[h,dataPCA,measlistact]  = Peak_Marking_PCA(intensnorm,zone,t(indt),measlistact );    %d1ok segment of normalise intensity                            
%dataPCA with the correction apply if used
    %measlistact bloc mark as noise
    indtoremove = find(measlistact==0);
    noise(indt,indtoremove,1)=4;
    intensnormall(indt,:) = dataPCA;
end
   save('testPCA','intensnormall','noise','-mat')

