function out = nirs_run_readAudio(job)

load(job.NIRSmat{1},'-mat') 
for ifile=1:numel(job.Audio_files)
    NIRS.Dt.Audio.pp.p{ifile} = job.Audio_files{ifile};
     %when the data are normalised the NIRS data and video offset to be synchronised 
    if job.Audio_files_sync ==0 %EEG
        NIRS.Dt.Audio.syncref = 'EEG';
        NIRS.Dt.Audio.pp.sync_timesec{ifile}= nan ;    
    elseif job.Audio_files_sync ==1 %NIRS
        NIRS.Dt.Audio.syncref = 'NIRS';
        NIRS.Dt.Audio.pp.sync_timesec{ifile}= 0 ;
    elseif job.Audio_files_sync ==2
        NIRS.Dt.Audio.syncref = 'AUX';
        NIRS.Dt.Audio.pp.sync_timesec{ifile}= nan ;
    end
     if isfield(job.c_Audiooffset,'b_Audiooffset_yes')
        NIRS.Dt.Audio.pp.offset{ifile} = job.c_Audiooffset.b_Audiooffset_yes.i_Audiooffset(ifile);
    else        
        NIRS.Dt.Audio.pp.offset{ifile} = 0
    end
end

save(job.NIRSmat{1},'NIRS'); 

out.NIRSmat = job.NIRSmat;