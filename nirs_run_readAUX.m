function out = nirs_run_readAUX(job)
    
 
load(job.NIRSmat{1},'-mat'); 
%remove old.... 
if isfield(NIRS.Dt,'AUX')
    NIRS.Dt = rmfield(NIRS.Dt,'AUX');
end
if 1 ==job.AUX_files_paired_NIRS ;
     for ifile=1:numel(job.AUX_files)
        NIRS.Dt.AUX(1).pp.p{ifile} = job.AUX_files{ifile};
        disp(['Link AUX: ',  job.AUX_files{ifile}])
        [pathstr, name, ext]=fileparts(job.AUX_files{ifile});
        NIRS.Dt.AUX.label = name;
    end
else
    for ifile=1:numel(job.AUX_files)
        NIRS.Dt.AUX(ifile).pp.p{1} = job.AUX_files{ifile};
        disp(['Link AUX: ',  job.AUX_files{ifile}])
        [pathstr, name, ext]=fileparts(job.AUX_files{ifile});
        NIRS.Dt.AUX(ifile).label = name;
    end
end
save(job.NIRSmat{1},'NIRS');



out.NIRSmat = job.NIRSmat;