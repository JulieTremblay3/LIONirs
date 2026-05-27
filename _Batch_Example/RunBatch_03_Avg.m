%parameter to modify for each subject
fileXLS = 'C:\Data\ToolboxLIONirs\TutorialHandsOn2026_Visual\BatchExample\TemplateStudyINFO.xlsx';
T = readtable(fileXLS);
jobfile = 'TemplateBatch_03_AVG.m';                                         %Define the job to run
%Enter Undefine information by participant
for isubject=1:size(T,1) %ligne du xls à lancé ici      
    jobid = cfg_util('initjob',fullfile(filepath,jobfile));           
    %MODIFY THE SUBJECT information according to the good xls column
    NIRSmat = {fullfile([T.AnalysisFolder{isubject},'/dCONC'], 'NIRS.mat')} %Define nirs.mat entry point
    BRANCH= T.Branch{isubject} ;                                            %Define name for each average branch
    TRIGVAL = T.Trig(isubject);                                             %Define trig   
    cfg_util('filljob', jobid,  NIRSmat,  BRANCH, TRIGVAL);                 %Enter all undefined parameter in order
    cfg_util('run', jobid)                                                  %RUN THE JOB
    disp(['Subject : ', raw{isubject,1},' Done'])      
end