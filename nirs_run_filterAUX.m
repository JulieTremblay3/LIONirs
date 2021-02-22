%___________________________________________________________________
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine
% www.lionlab.umontreal.ca
%___________________________________________________________________
function out = nirs_run_filterAUX(job)
% Filter AUX (EEG) data based on parameters already used for the filtering
% of NIRS data. Create new .dat files and update the NIRS.mat info
% regarding the AUX files.
% INPUT:
%       job: structure with fieldname NIRSmat, containing at least one path+file of a NIRSmat file
%            eg.job.NIRSmat={'C:...\NIRS.mat'};

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
    NIRS = [];
    load(job.NIRSmat{filenb});
    
    
    liststeps={NIRS.Dt.fir.pp.pre};
    for xx = 1:length(liststeps)
        if contains(liststeps{xx},'Filter')
            fstep=xx;
        end
    end
    if ~fstep
        error('No filter applied to the NIRS data %s/ntherefore no filtering parameters', (job.NIRSmat{filenb}))
    end
    
    %parameters
    [lowcut,applylowcut] = str2num(NIRS.Dt.fir.pp(fstep).job.lowcutfreq);
    [highcut ,applyhighcut] = str2num(NIRS.Dt.fir.pp(fstep).job.highcutfreq);
    paddingsym = NIRS.Dt.fir.pp(fstep).job.paddingsymfilter; %symetrie padding on the signal to avoid edge on the filtering
    interpolate = NIRS.Dt.fir.pp(fstep).job.interpolatebadfilter;  %symetrie padding on the signal to avoid edge on the filtering
    filt_ord = NIRS.Dt.fir.pp(fstep).job.filterorder;
    DelPreviousData  = NIRS.Dt.fir.pp(fstep).job.DelPreviousData;
    
    %use last step of preprocessing
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N; %number channels
    fs = NIRS.Cf.dev.fs; %sampling rate
    
    if ~isfield(NIRS.Dt,'AUX')
        error('No auxiliary attached to the NIRS file %s', (job.NIRSmat{filenb}))
    end
    
    for iAUX = 1:numel(NIRS.Dt.AUX)
        newAUX=numel(NIRS.Dt.AUX)+1;
        if contains(NIRS.Dt.AUX(iAUX).label,'AUX') || contains(NIRS.Dt.AUX(iAUX).label,'EEG')
            fprintf('AUX found in NIRS.mat. Creating new AUX .dat file (downsampled, filtered, normalized)...\n')
            if ~isfield(NIRS.Dt.AUX(iAUX).pp(end),'sync_timesec') %check if there were synchronisation -take the last field
                disp('No segmentation have been made -- ensure that aux synchronisation are ok')
                disp('FilterAUX not completed')
            else
                for i=1:length(NIRS.Dt.AUX(iAUX).pp(end).p)
                    nameAUX=NIRS.Dt.AUX(iAUX).pp(end).p{i};
                    try
                        tstart=NIRS.Dt.AUX(iAUX).pp(end).sync_timesec{i};
                    catch
                        tstart = 0;
                    end
                    
                    timeSECaxe = (1/fs):(1/fs):(NIRS.Dt.fir.sizebloc{i}*1/fs);
                    tstop = tstart+timeSECaxe(end);
                    [cPATH,cFILE,cEXT]=fileparts(nameAUX); %current path file and extension
                    [data,infoBV,marker,ind_dur_ch] = fopen_EEG(nameAUX, tstart, tstop);
                    
                    fsAUX =1/(infoBV.SamplingInterval/1000000); %Frequence echantillonage Hz
                    [~,q] = rat(fs/fsAUX,0.0001);
                    AUXupdate.ind_dur_ch=ind_dur_ch((ind_dur_ch(:,1)>tstart*fsAUX) & (ind_dur_ch(:,1)<tstop*fsAUX),:);
                    AUXupdate.marker=marker((ind_dur_ch(:,1)>tstart*fsAUX) & (ind_dur_ch(:,1)<tstop*fsAUX),:);
                    
                    
                    for ich=1:numel(infoBV.name_ele)
                        
                        tmp = data(:,ich);
                        
                        
                        %downsample to the same as fs NIRS
                        tmpr=downsample( tmp , q);
                        
                        
                        %FILTER
                        if paddingsym
                            dtmp = [fliplr(tmpr);tmpr;fliplr(tmpr)];
                            tstartd = size(tmpr,1)+1;
                            tstopd = size(tmpr,1)*2;
                            d = dtmp;
                        else
                            d=tmpr;
                        end
                        id = find(isnan(d));
                        if ~isempty(id)
                            d(id) = 0;
                        end
                        
                        %bandpass
                        if applylowcut && applyhighcut %band pass (low and high)
                            W2 = lowcut*2/fs;
                            W1 = highcut*2/fs;
                            [fb,fa]=butter(filt_ord,[W1, W2]);
                            dfilt = filtfilt(fb,fa,d);
                        elseif applylowcut==true && applyhighcut==false %only low pass
                            Wn = lowcut*2/fs;
                            [fb,fa]=butter(filt_ord,Wn);
                            dfilt = filtfilt(fb,fa,d);
                        elseif applylowcut==false && applyhighcut==true %only high pass
                            Wn = highcut*2/fs;
                            [fb,fa]=butter(filt_ord,Wn,'high');
                            dfilt = filtfilt(fb,fa,d);
                        end
                        
                        if paddingsym
                            dfilt = dfilt(tstartd:tstopd);
                        end
                        
                        %END OF FILTER SECTION
                        
                        %normalize the AUX data (z score)
                        tmpn=(dfilt-mean(dfilt))/std(dfilt);
                        
                        % we cut aux to the data initial size
                        if numel(tmpn)<numel(timeSECaxe)
                            nplus = numel(timeSECaxe)-numel(tmpn);
                            try
                                tmpn = [tmpn; tmpn(end-nplus:end) ];
                            catch
                                try
                                    tmpn = [tmpn, tmpn(end-nplus:end) ];
                                catch
                                    
                                    msgbox(sprintf('Too short AUX\n%s b%d',cFILE,i))
                                end
                            end
                        elseif numel(tmpn)>numel(timeSECaxe)
                            tmpn = tmpn(1:numel(timeSECaxe));
                        end
                        
                        dataNEW(:,ich)=tmpn;
                    end
                    %OVERWRITE INFOS in infoBV
                    infoBV.DataPoints=size(dataNEW,1);
                    infoBV.SamplingInterval=(1/fs)*1000000;
                    
                    AUXupdate.infoBV=infoBV;
                    AUXupdate.data=dataNEW;
                    
                    if ~exist([cPATH '\filAUX\'],'dir')
                        mkdir([cPATH '\filAUX\'])
                    end
                    fileoutAUX=[cPATH '\filAUX\' cFILE 'b' num2str(i) cEXT];
                    fwrite_EEG(fileoutAUX,AUXupdate,1,AUXupdate.infoBV.DataPoints );
                    disp(fileoutAUX)
                    NIRS.Dt.AUX(newAUX).pp.p{i,1}=fileoutAUX;
                    NIRS.Dt.AUX(newAUX).pp.sync_timesec{i,1}=0;
                end
                NIRS.Dt.AUX(newAUX).label = ['fil' NIRS.Dt.AUX(iAUX).label ];
                
            end
        end
    end
    
    fprintf('Update NIRSmat COMPLETED ...%s\n*\n**\n***\n',job.NIRSmat{filenb})
    save(job.NIRSmat{filenb,1},'NIRS');
    
end
out.NIRSmat = job.NIRSmat;
end
