filesetting = 'C:\data\LFC01_SM\Trig_Spike_small.m'
%Attent des données .nir outpout de la lecture de SPM
run(filesetting)
file = [pathname,filename{1}];
filevhdr = [file(1:end-3),'vhdr']
info = read_vhdr(filevhdr)
d = fopen_NIR(file,info.NumberOfChannels);
samplingintervals = info.SamplingInterval/1000000;
pAstart = 1; % round(Astart/samplingintervals);
pAstop =size(d,2);%round(Astop/samplingintervals);  
t = pAstart*samplingintervals:samplingintervals:pAstop*samplingintervals;
datatime =timingfile{1};
d = d';
figure
plot(t,d)
title('no filter')
AdvOptions.FiltOptions.FilterType=1;  %Default Butterworth
AdvOptions.FiltOptions.FilterOrder=3;  %Default 3rd order
lpf =3;
fs = 1/(t(2)-t(1)); 
[fb,fa] = MakeFilter(AdvOptions.FiltOptions.FilterType,AdvOptions.FiltOptions.FilterOrder,fs,lpf,'low');
d = filtfilt(fb,fa,d);  

figure
plot(t,d)
title('filter') 
for istim=1:size(datatime,1)
%Segment par 
t1=datatime(istim,1); %temps debut en seconde
t2=datatime(istim,2);%temps debut en seconde
t1i = find(t<=t1);
t1i=t1i(end);
t2i = find(t<=t2);
t2i=t2i(end);
indt = t1i:t2i;
meanIntens = median(d(indt,:)); 
intensnorm = d(indt,:)./(ones(numel(indt),1)*meanIntens); %intensite normalise non filtré 
measlistact = ones(size(intensnorm,2),1);
indbad =find(mean(abs(intensnorm-1))> 0.03);
measlistact(indbad) = 0;
plotlst = find(measlistact);
figure
plot(t(indt),intensnorm(:,plotlst))
zone = correlation_channel_gen(intensnorm,measlistact); %ok
save('corr.zone','zone','-mat')
[h,dataPCA,measlistact]  = Peak_Marking_PCA(intensnorm,zone,t(indt),measlistact );    %d1ok segment of normalise intensity    
figure
plot(t(indt),dataPCA(:,plotlst))
end 