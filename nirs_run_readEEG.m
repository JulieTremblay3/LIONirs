function out = nirs_run_readEEG(job)

load(job.NIRSmat{1},'-mat'); 
job.EEG_files;
for ifile=1:numel(job.EEG_files)
    NIRS.Dt.EEG.pp.p{ifile} = job.EEG_files{ifile};
    disp(['Link EEG: ', job.EEG_files{ifile}])
end

save(job.NIRSmat{1},'NIRS'); 

out.NIRSmat = job.NIRSmat;