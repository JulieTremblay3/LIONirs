function out = nirs_run_readEEGMarker(job)

load(job.NIRSmat{1},'-mat');
NC = NIRS.Cf.H.C.N;
rDtp = NIRS.Dt.fir.pp(end).p; 
fs = NIRS.Cf.dev.fs;

for ifile =1:numel( NIRS.Dt.EEG.pp(end).p)
    nameEEG = NIRS.Dt.EEG.pp(end).p{ifile};
    if isfield(NIRS.Dt.EEG.pp(end),'sync_timesec')
        offset = NIRS.Dt.EEG.pp(end).sync_timesec{ifile};    
    else
        msgbox('Please first perform segmentation to ensure synchronization between fNIRS and EEG')
        return
    end
    d = fopen_NIR(rDtp{ifile,1},NC);
    timeNIRS = 1/fs:1/fs:1/fs*(size(d,2));
    clear d;
    [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch] = fopen_EEG(nameEEG);
    timeEEG = EEG.infoBV.SamplingInterval/1e6:EEG.infoBV.SamplingInterval/1e6:EEG.infoBV.SamplingInterval/1e6*size(EEG.data,1);
    clear EEG.data
    str = job.e_readEEGMarker_Marker;
    remain = str;
    while ~isempty(remain);
        [token,remain] = strtok(remain,',');
        [markerlabel, newtrig]=strtok(token,':');
        trigmk = str2num(newtrig(2:end));
        mklist = [];
        for imk=2:size(EEG.marker,1)
            if strcmp(upper(EEG.marker{imk,2}),upper({markerlabel}));
                mklist = [mklist,imk];
            end
        end
        %Transfert mklist in the AUX.
        tmarker = timeEEG(EEG.ind_dur_ch(mklist,1))-offset; % sample in EEG samping rate must be convert fNIRS closest point
        for imk= 1:numel(tmarker)
            if tmarker(imk)>0
                id = find(timeNIRS<=tmarker(imk));
                if ~isempty(id)
                    NIRS.Dt.fir(end).aux5{ifile} =[  NIRS.Dt.fir(end).aux5{ifile} ;trigmk,id(end)];
                end
            end
        end
        disp(['MARKER: ', markerlabel,' where transfert in trig ', num2str(trigmk) ,' at fNIRS time :', num2str(  tmarker )  ])
    end
end


save(job.NIRSmat{1},'NIRS');

out.NIRSmat = job.NIRSmat;