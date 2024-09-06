function out = nirs_run_NIRSmatdiradjust(job)
for imat =1:numel(job.NIRSmat)
load(job.NIRSmat{imat,1});
[NIRSmatnewdir,name,ext] = fileparts(job.NIRSmat{imat,1});

for imodule = numel(NIRS.Dt.fir.pp):-1:1
    if imodule == 2
        1;
    end
    for ifile = 1:numel(NIRS.Dt.fir.pp(imodule).p)
    %Find the new root dir 
    try
        fileloc_previous = NIRS.Dt.fir.pp(imodule).p{ifile};
    catch
        fileloc_previous = NIRS.Dt.fir.pp(imodule).p
    end
        [dir1,file1,ext1]=fileparts(fileloc_previous);
    %Descendre les r�pertoire pour trouver la partie qui doit �tre modifier
    %
    [toknew,remnew]=strtok(fliplr(NIRSmatnewdir),{'\','/'});
    [tok1,rem1]=strtok(fliplr(dir1),{'\','/'});
    ok = 1;
    while ok
     if  strcmp(toknew, tok1)
         dirrootnew = [fliplr(remnew),fliplr(toknew)];
        ok=0; %getout ! 
     else
         [toknew,remnew]=strtok(remnew,{'\','/'});
         dirrootnew = [fliplr(remnew),fliplr(toknew)];
     end 
     if isempty(toknew)
         ok=0; 
         disp([fliplr(tok1) ,' root folder do not exist'])
     end
         
    end 
     NIRS.Dt.fir.pp(imodule).p{ifile} = fullfile( dirrootnew , [file1,ext1]); 
    end                
end 
disp(['Folder adjustment location: ', dirrootnew]); 
if isfield(job.c_MultimodalPath,'b_MultimodalPath_yes')    
    dirrootnew = job.c_MultimodalPath.b_MultimodalPath_yes.e_MultimodalPath{1};

    if isfield(NIRS.Dt, 'AUX')     
        
                for iAUX=1:numel(NIRS.Dt.AUX)      
                    for imodule = 1:numel(NIRS.Dt.AUX(iAUX).pp)
                        for ifile = 1:numel(NIRS.Dt.AUX(iAUX).pp(imodule).p)
                                fileloc_previous = NIRS.Dt.AUX(iAUX).pp(imodule).p{ifile};
                                [dir1,file1,ext1]=fileparts(fileloc_previous);
                                NIRS.Dt.AUX(iAUX).pp(imodule).p{ifile} = fullfile( dirrootnew , [file1,ext1]);
                        end                    
                    end 
                end
    end
    if isfield(NIRS.Dt, 'EEG')  
        for iAUX=1:numel(NIRS.Dt.EEG)      
                    for imodule = 1:numel(NIRS.Dt.EEG(iAUX).pp)
                        for ifile = 1:numel(NIRS.Dt.EEG(iAUX).pp(imodule).p)
                                fileloc_previous = NIRS.Dt.EEG(iAUX).pp(imodule).p{ifile};
                                [dir1,file1,ext1]=fileparts(fileloc_previous);
                                NIRS.Dt.EEG(iAUX).pp(imodule).p{ifile} = fullfile( dirrootnew , [file1,ext1]);
                        end 
                    end 
        end
    end
    if isfield(NIRS.Dt, 'VIDEO')  
%         for iAUX=1:numel(NIRS.Dt.EEG)      
%                     for imodule = 1:numel(NIRS.Dt.EEG(iAUX).pp)
%                         for ifile = 1:numel(NIRS.Dt.EEG(iAUX).pp(imodule).p)
%                                 fileloc_previous = NIRS.Dt.EEG(iAUX).pp(imodule).p{ifile};
%                                 [dir1,file1,ext1]=fileparts(fileloc_previous);
%                                 NIRS.Dt.EEG(iAUX).pp(imodule).p{ifile} = fullfile( dirrootnew , [file1,ext1]);
%                         end 
%                     end 
%         end
    end
    
        if isfield(NIRS.Dt, 'VIDEO')  
%         for iAUX=1:numel(NIRS.Dt.EEG)      
%                     for imodule = 1:numel(NIRS.Dt.EEG(iAUX).pp)
%                         for ifile = 1:numel(NIRS.Dt.EEG(iAUX).pp(imodule).p)
%                                 fileloc_previous = NIRS.Dt.EEG(iAUX).pp(imodule).p{ifile};
%                                 [dir1,file1,ext1]=fileparts(fileloc_previous);
%                                 NIRS.Dt.EEG(iAUX).pp(imodule).p{ifile} = fullfile( dirrootnew , [file1,ext1]);
%                         end 
%                     end 
%         end
    end
end

save(job.NIRSmat{imat,1},'NIRS');
clear NIRS
end
out.NIRSmat = job.NIRSmat;