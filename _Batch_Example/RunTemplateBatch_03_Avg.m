%parameter to open by subject
fileXLS = 'C:\Data\ToolboxLIONirs\TutorialHandsOn2026_Visual\BatchExample\StudyINFO.xlsx';
[num,txt,raw] = xlsread(fileXLS);
labelcol = raw(1,:);
%BATCH TO RUN
jobfile = 'Batch_03_AVG.m';

fullpath = mfilename('fullpath');
[filepath,name,ext] = fileparts(fullpath);
%Enter Undefine information by participant
for isubject=3:size(raw,1) %ligne du xls à lancé ici

    % load template batch saved as as script TEMPLATE_BATCHdCONC.m');
    %try      
    jobid = cfg_util('initjob',fullfile(filepath,jobfile));
           
    %MODIFY THE SUBJECT information according to the good xls column
    Subject(1).NIRSmat = {fullfile([raw{isubject,Col_ID(labelcol,'AnalysisFolder')},'/dCONC'], 'NIRS.mat')}
      
    cfg_util('filljob', jobid, Subject(1).NIRSmat);


    %RUN THE JOB
    cfg_util('run', jobid)
    disp(['Subject : ', raw{isubject,1},' Done'])   
   
%    
    % catch
    %     disp(['ERROR isubject ' raw{isubject,1},' fail'] );
    % end
    
end