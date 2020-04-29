function [SD,varargout] = nirs_boxy_convert2(SD,fName,fileOut,Idx_file,preprocess)
%This function reads the BOXY txt file format and extracts the data. 
%The varargin input is for fast preprocessing (preprocess, noisy,
%h_amp, amp_max,cf_max). The varargout returns the list of good channels.
% 
% if SD.save_bin1 %need to delete, as file will be appended
%     if exist(fileOut,'file'), delete(fileOut); end
% end
flagACPH = 1;
fid=fopen(fName,'r');
%hwaitbar = waitbar(0,'Open raw file');
if fid == -1
   %Test valid
   warndlg(['Cannot open file ', fName])
   error('Cannot open file');
end    
fseek(fid, 0, 'eof');
endoffile = ftell(fid);
fseek(fid, 0, 'bof');
%Now find the start/stop of data
warning('off','MATLAB:dispatcher:InexactCaseMatch');
while(1)
   tline=fgets(fid);
   if  ~isempty(strfind(tline,'#DATA'))
       startDataPos=ftell(fid);
       break;
   end
   if ~isempty(strfind(tline,'#ACQ '))
       ACQoptions=readACQ(fid);
   end
end %end while(1)
warning('on','MATLAB:dispatcher:InexactCaseMatch');

fclose(fid);

numSrc = size(SD.PosLasers,1);
ml=[]; 

%To concatenate several files
DataPoints = 1;
SD.Markers{1,1}=[];  
SD.fs = ACQoptions.Fs;

SD.nb_Mux = ACQoptions.MUXChannels;
fid=fopen(fName,'r');
if fid == -1
     msgbox([fName, 'doesn''t exist'])
   return
end
while(1)
   tline=fgetl(fid);
   if  ~isempty(strfind(tline,'#DATA'))
       endDataPos=ftell(fid);
       break;
   end
end
%nb_file = 1; %SD.N_dfiles;
%Create a New Segment marker
SD.Markers{Idx_file,1}.Type='New Segment';
SD.Markers{Idx_file,1}.Description='';
SD.Markers{Idx_file,1}.Position=DataPoints;
SD.Markers{Idx_file,1}.Size=1;
SD.Markers{Idx_file,1}.ChNumber=0; %all channels

%initialise the new data  
%GO to the start position and read out the data
frewind(fid);
fseek(fid,startDataPos,-1);
MeasMap=readMeasMap(fid);
nMeasMap=AdjustMeasMap(ACQoptions,MeasMap);
MeasMap=nMeasMap;
linetest=fgetl(fid); %One more blank line
DetectorLabels = char((1:ACQoptions.MUXChannels)+64);  %"A","B","C"

%load all the data
data = textscan(fid,'%f');   
%shape into channels x time
data = reshape(data{1},size(MeasMap,1),[]);

listdownsample = ['record', 'exmux','time','group','flag'  ];
listmaxdownsample = ['digaux','mark'];
ml=[];
d_dc = []; %DC component
d_ph = []; %PH component
d_ac = []; %AC component
if SD.resample >1
    %introduce other array, in case we are downsampling
    data1 = [];
    for i_ch=1:size(data,1)               
        if strfind(listdownsample,deblank(MeasMap(i_ch,:)))
           data1(i_ch,:)= downsample(data(i_ch,:),SD.resample);                   
        elseif strfind(listmaxdownsample,deblank(MeasMap(i_ch,:)))
            nb_sample = numel(data(i_ch,:));
            rest = mod(nb_sample,SD.resample); 
            if rest ~= 0, rest = SD.resample-rest; end
            nbrow = ceil(nb_sample/SD.resample); %not defined if resample = 1!
            data1(i_ch,:)= max(reshape([data(i_ch,:),zeros(1,rest)], SD.resample,nbrow));
        else
            %Note that decimate is extremely slow - but it does include a
            %low pass filter
            data1(i_ch,:)  = decimate(data(i_ch,:),SD.resample);
        end
    end %end for 
    data = data1;
end %end if SD.resample

for Sidx = 1:numSrc
   Detectors = SD.Pairs{Sidx,1};  
   lasers = SD.PosLasers{Sidx,1};
   for LaserIdx=1:length(lasers)% de 1 à 2
       for Didx=1:length(Detectors)% de 1 à 8 ou 16
           %Check mux mode used %JT
           muxindex = mod(lasers(LaserIdx), SD.nb_Mux);
           if(muxindex==0) %important ex : mod(16,16) = 0 alors que c'est la source 16 qui doit être prise 
               muxindex = SD.nb_Mux;
           end      
           strDC = [DetectorLabels(Detectors(Didx)) '-DC' num2str(muxindex)];
           if flagACPH
           strPh = [DetectorLabels(Detectors(Didx)) '-Ph' num2str(muxindex)];
           strAC = [DetectorLabels(Detectors(Didx)) '-AC' num2str(muxindex)];
           end
           % Match the detector and laser in the measurements
           dcIdx = strmatch(strDC,MeasMap,'exact');
           if flagACPH
                phIdx = strmatch(strPh,MeasMap,'exact');
                acIdx = strmatch(strAC,MeasMap,'exact');
           end
           ml=[ml; lasers(LaserIdx) Detectors(Didx) 1 0];                   
           d_dc=[d_dc; data(dcIdx,:)];  %DC data 
           if flagACPH
                d_ph=[d_ph; data(phIdx,:)];  %phase data
                d_ac=[d_ac; data(acIdx,:)];  %AC data
           end
       end
   end %end for LaserIdx=1:length(lasers)
   Sidx = Sidx + 1; 
end %end for Sidx   
% % Store digital aux 
DigIdx = strmatch('digaux',MeasMap,'exact');
SD.aux5 = sparse(data(DigIdx,:));

id =strmatch('aux-1',MeasMap,'exact');
if ~isempty(id); SD.aux1 = data(id,:);end
id =strmatch('aux-2',MeasMap,'exact');
if ~isempty(id); SD.aux2  = data(id,:);end
id = strmatch('aux-3',MeasMap,'exact');
if ~isempty(id);SD.aux3  = data(id,:);end
id = strmatch('aux-4',MeasMap,'exact'); 
if ~isempty(id);SD.aux4 = data(id,:);end
id =strmatch('digaux',MeasMap,'exact');
if ~isempty(id); SD.aux5i = data(id,:);end

% DigIdx = strmatch('aux-1',MeasMap,'exact');
% Aux4 = sparse(data(DigIdx,:));
% max(data)
% [a,id]=max(data')
% a(55)=0
% figure; plot(a)
% [b,id]=max(a)c
% figure;
% fs = 19.5312  
% t = (1/fs):(1/fs): (numel(Aux3)*1/fs);
% plot(t,Aux3)
% figure;
% fs = 19.5312  
% t = (1/fs):(1/fs): (numel(Aux4)*1/fs);
% plot(t,Aux4)
% figure; plot(t,Aux2)



% figure
% hold on
% h=fdesign.lowpass('Fp,Fst,Ap,Ast',0.3,0.6,1,60);
% d=design(h,'equiripple'); %Lowpass FIR filter
% y=filtfilt(d.Numerator,1,Aux2); %zero-phase filtering
% plot(Aux1,'r','displayname','aux1')
% plot(Aux2,'b','displayname','aux2')
% plot(Aux3,'k','displayname','aux3')
% plot(Aux4,'g','displayname','aux4')
% % plot(data(DigIdx,:))
% SD.aux5 = sparse(data(DigIdx,:)); 
%create measure list 
for i=1:numel(SD.PosLasers)    
    numWaveall(i) = size(SD.PosLasers{i,1},2);
end
numWave =  max(numWaveall);
ml_tmp=ml;
for idx=1:numSrc
    for idx2=1:numWave
        if ~isempty(SD.PosLasers{idx,1})
            ml_tmp(ml(:,1)==SD.PosLasers{idx,1}(idx2),1)=idx;
            ml_tmp(ml(:,1)==SD.PosLasers{idx,1}(idx2),4)=idx2;
        end
    end
end
ml = ml_tmp; 
%now sort the data
[~, sortlst]=sort(ml(:,4));
ml=ml(sortlst,:);
%Ajout des aux dans les données a la suite de 830 et 690
SD.ml=ml;

%sort the data
try
    d_dc = d_dc(sortlst,:); 
catch
    msgbox('Please verify mux, number of detector and number of source')
    return
end
if flagACPH
    d_ph=d_ph(sortlst,:);  %phase data
    d_ac=d_ac(sortlst,:);  %AC data
end


    if  preprocess.STDampenable || preprocess.DCampenable
        Fs = ACQoptions.Fs; %Sampling rate
        [d_dc,ok] = fast_preprocessing(d_dc,preprocess,Fs);
    else
        ok = ones(size(d_dc,1),1);
    end
    varargout{1} = ok;    

SD.sizebloc = size(d_ph,2);



%save binary files 
    fwrite_NIR(fileOut,d_dc);
    if flagACPH
        fwrite_NIR([fileOut(1:end-4),'PH','.nir'],d_ph);
        fwrite_NIR([fileOut(1:end-4),'AC','.nir'],d_ac);     
    end
end             


%***************************************
function ACQoptions=readACQ(fid)
%This function parses the file acquisition settings
ACQoptions=[];
while(1)
   line= fgets(fid);
   if ~isempty(strfind(line,'#DATA'))
       fseek(fid,-20,'cof');
       break;
   end
   if ~isempty(line) && ~all(isspace(line))
       [value1, attribute] = strtok(line);
       attribute(isspace(attribute))=[];
       value1=str2num(value1);
       switch(attribute)
           case 'DetectorChannels'
               ACQoptions.DetChannels=value1;
           case 'ExternalMUXChannels'
               ACQoptions.MUXChannels=value1;
           case 'AuxiliaryChannels'
               ACQoptions.AUXChannels=value1;
           case 'Waveform(CCF)Frequency(Hz)'
               ACQoptions.WaveformFreq=value1;
           case 'WaveformsSkipped'
               ACQoptions.WaveformSkip=value1;
           case 'WaveformsAveraged'
               ACQoptions.WaveformAvg=value1;
           case 'CyclesAveraged'
               ACQoptions.CyclesAvg=value1;
           case 'AcquisitionsperWaveform'
               ACQoptions.ACQperWaveform=value1;
           case 'UpdataRate(Hz)'
               ACQoptions.Fs=value1;
       end
   end
end
end


%***************************************************************
function MeasMap=readMeasMap(fid)
%This returns a mapping of the measurements
line=fgetl(fid);
whitespaces=regexp(line,'\s');
whitespaces=[0 whitespaces];
MeasMap=[];
for idx=1:length(whitespaces)-1
   MeasMap=strvcat(MeasMap,line(whitespaces(idx)+1:whitespaces(idx+1)));
end
end

function MeasMap=AdjustMeasMap(Acq,oldMeasMap)
ndet = Acq.DetChannels;
nmux = Acq.MUXChannels;
nmmap=size(oldMeasMap,1);
MeasMap=[];

for index=1:nmux
   if(index==1) % If first line, take all but only number detectors
       for jindex=1:nmmap
           if(jindex>2 && jindex<= 2+3*ndet)
               MeasMap=strvcat(MeasMap, [strtrim(oldMeasMap(jindex,:)),num2str(index)]);
           else
               MeasMap=strvcat(MeasMap, strtrim(oldMeasMap(jindex,:)));
           end
       end
   else
       for jindex=1:(2+3*ndet) % Only detectors in this case
           MeasMap=strvcat(MeasMap, [strtrim(oldMeasMap(jindex,:)),num2str(index)]);
       end
   end
end
end



