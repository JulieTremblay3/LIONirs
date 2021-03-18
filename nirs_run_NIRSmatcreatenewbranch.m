function out = nirs_run_NIRSmatcreatenewbranch(job)

load(job.NIRSmat{1,1});
job.e_NIRSmatdirnewbranch;
[NIRSmatnewdir,name,ext] = fileparts(job.NIRSmat{1,1});



dir2 = [NIRSmatnewdir,filesep,job.e_NIRSmatdirnewbranch];
file2= fullfile(dir2,'NIRS.mat' );
if ~isfolder(dir2)
    mkdir(dir2);
else
    try
        rmdir(dir2, 's');
        mkdir(dir2);
    catch
        disp(['Error unable to remove' , dir2])
        return
    end
end
save(file2,'NIRS');
if job.m_newbranchcomponent
    try
        infileSelectedfactors = fullfile(NIRSmatnewdir,'SelectedFactors.mat');
        outfileSelectedfactors = fullfile(dir2,'SelectedFactors.mat');
        copyfile(infileSelectedfactors,outfileSelectedfactors);
    catch
    end
    try
        infileCorrectionApply = fullfile(NIRSmatnewdir,'CorrectionApply.mat');
        outfileCorrectionApply = fullfile(dir2,'CorrectionApply.mat');
        copyfile(infileCorrectionApply,outfileCorrectionApply);
    catch
    end
    
end

out.NIRSmat = {file2};