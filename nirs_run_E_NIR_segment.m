function out = nirs_run_E_NIR_segment(job)

out = [];
 NIRS = [];
 load(job.NIRSmat{1,1});
 NC = NIRS.Cf.H.C.N;   
 fs = NIRS.Cf.dev.fs;
 d = fopen_NIR(job.NIR_FileIn{1},NC);
 t = 1/fs:1/fs:1/fs*size(d,2); 
 id_Start = find(t<job.NIR_START_TIME);
 if isempty(id_Start)
     id_Start = 1;
 end
  id_Stop =  find(t<job.NIR_STOP_TIME);
 
 Data = d(:,id_Start(1):id_Stop(end));
 [pathstr, name, ext, versn] = fileparts(job.NIR_FileIn{1}) ;
 outfile = [pathstr,'\','cut_',num2str(job.NIR_START_TIME),'_to_',num2str(job.NIR_STOP_TIME),'_s',name, ext];
fwrite_NIR(outfile, Data);
fprintf('%s\n',outfile);