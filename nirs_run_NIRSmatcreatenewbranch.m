function out = nirs_run_NIRSmatcreatenewbranch(job)

load(job.NIRSmat{1,1});
[NIRSmatnewdir,name,ext] = fileparts(job.NIRSmat{1,1});

   

dir2 = [NIRSmatnewdir,filesep,job.e_NIRSmatdirnewbranch];
file2= fullfile(dir2,'NIRS.mat' );
disp(['Create new branch: ',file2 ])

if ~isdir(dir2)
    mkdir(dir2);
else
    try
        rmdir(dir2, 's');
        mkdir(dir2);
        
    catch
        disp(['Error unable to write to ' , dir2, ' verify access or open file'])
        return
    end
end
save(file2,'NIRS');
if job.m_newbranchcomponent
    try
        infileSelectedfactors = fullfile(NIRSmatnewdir,'SelectedFactors.mat');
        outfileSelectedfactors = fullfile(dir2,'SelectedFactors.mat');
        copyfile(infileSelectedfactors,outfileSelectedfactors);
        disp(['Copy components: ',outfileSelectedfactors])
    catch
        disp(['No components to copy empty: ',outfileSelectedfactors]);
    end
    try
        infileCorrectionApply = fullfile(NIRSmatnewdir,'Apply.mat');
        outfileCorrectionApply = fullfile(dir2,'CorrectionApply.mat');
        copyfile(infileCorrectionApply,outfileCorrectionApply);
        disp(['Copy corrections: ',outfileCorrectionApply])
    catch
        disp(['No corrections to copy empty: ',outfileCorrectionApply]);     
    end
else
    disp(['No components & corrections in the new branch'])
end

out.NIRSmat = {file2};