function out = nirs_run_GUI(job)
% Call the GUI for data display.
%
subjectnb = 1; 
plot_sessions_GUI(job.NIRSmat,subjectnb,job);
% for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
%     plot_sessions_GUI(job.NIRSmat{filenb,1});
% end

out.NIRSmat = job.NIRSmat; 

