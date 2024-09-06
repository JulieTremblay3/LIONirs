function out = nirs_run_E_viewNIRS(job)

for filenb=1:size(job.NIRSmat,1)
      file.name = job.NIRSmat{filenb,1};
      job.NIRSmat{filenb,1} = NIRS_GUI(file);
end

out.NIRSmat = job.NIRSmat;