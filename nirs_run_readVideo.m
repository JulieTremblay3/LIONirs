function out = nirs_run_readVideo(job)
load(job.NIRSmat{1},'-mat')  
for ifile=1:numel(job.Video_files)
    NIRS.Dt.Video.pp(1).p{ifile} = job.Video_files{ifile};
    disp(['Link Video: ',  job.Video_files{ifile}])
     %when the data are normalised the NIRS data and video offset to be synchronised 
    if job.video_files_sync ==0; %EEG
        NIRS.Dt.Video.syncref = 'EEG';
        NIRS.Dt.Video.pp(1).sync_timesec{ifile}= nan ; 
        disp(['Sync Video using EEG time reference'])
    elseif job.video_files_sync ==1 %NIRS
        NIRS.Dt.Video.syncref = 'NIRS';
        NIRS.Dt.Video.pp(1).sync_timesec{ifile}= 0 ;
        disp(['Sync Video using NIRS time reference '])  
    elseif job.video_files_sync ==2
        NIRS.Dt.Video.syncref = 'AUX';
        disp(['Sync Video using AUX time reference '])
        NIRS.Dt.Video.pp(1).sync_timesec{ifile}= nan ;
    end
    if isfield(job.c_videooffset,'b_videooffset_yes')
        NIRS.Dt.Video.pp(1).offset{ifile} = job.c_videooffset.b_videooffset_yes.i_videooffset(ifile);
        disp(['Additionnal offset with reference: ', num2str(job.c_videooffset.b_videooffset_yes.i_videooffset(ifile)),' secondes']);
    elseif isfield(job.c_videooffset    ,'b_videooffsetwithtrig_yes')
        offset_sec = job.c_videooffset.b_videooffsetwithtrig_yes.i_videooffset(ifile);
        disp(['Additionnal offset with reference: ', num2str(job.c_videooffset.b_videooffsetwithtrig_yes.i_videooffset(ifile)),' secondes']);
        trig = NIRS.Dt.fir.aux5{1};
        trigid= job.c_videooffset.b_videooffsetwithtrig_yes.i_videooffsetwithtrig   
        id = find(trig(:,1)==trigid);   
        if isempty(id)
            disp('Warning: Could not synchronize video file')
            NIRS.Dt.Video.pp(1).offset{ifile} = 0;
        else
        triglag_sec = trig(id(1),2)*  1/NIRS.Cf.dev.fs;
        NIRS.Dt.Video.pp(1).offset{ifile} =   offset_sec-triglag_sec;
        end

    else        
        NIRS.Dt.Video.pp(1).offset{ifile} = 0;
    end
end
save(job.NIRSmat{1},'NIRS'); 
out.NIRSmat = job.NIRSmat;