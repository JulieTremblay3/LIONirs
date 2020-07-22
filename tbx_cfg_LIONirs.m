%___________________________________________________________________
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine 
% www.lionlab.umontreal.ca
%___________________________________________________________________
function nirsHSJ = tbx_cfg_LIONirs


nirOutputDefaultPath = fullfile('C:\data\'); % Default output path for .nir files. Temporary solution.

addpath(fileparts(which(mfilename))); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 0 Montage project creation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_IOMTG         = cfg_exbranch;
E_IOMTG.name    = 'Montage configuration (.prj)';
E_IOMTG.tag     = 'E_IOMTG';
E_IOMTG.val     = {};
E_IOMTG.prog    = @nirs_run_E_IOMTG;
E_IOMTG.help    = {'Use to create project configuration file using source and detector position.'};


%General entry use by any modules
NIRSmat         = cfg_files; 
NIRSmat.name    = 'NIRS.mat';
NIRSmat.tag     = 'NIRSmat';      
NIRSmat.filter  = 'mat';
NIRSmat.ufilter = '^NIRS.mat$';    
NIRSmat.num     = [1 Inf];     % Number of inputs required 
NIRSmat.help    = {'Select NIRS.mat for the subject.'}; % help text displayed

prjfile         = cfg_files;
prjfile.name    = 'Helmet .prj'; 
prjfile.tag     = 'prjfile';       
prjfile.filter  = 'prj';
prjfile.ufilter = '.prj$';    
prjfile.num     = [1 Inf];     
prjfile.help    = {'Select project of the montage HSJ for the subject(s).'}; 


TOPOmat         = cfg_files; %Select NIRS.mat for this subject 
TOPOmat.name    = 'TOPOmat';       % The displayed name
TOPOmat.tag     = 'TOPOmat';       %file names
TOPOmat.filter  = 'mat';
TOPOmat.ufilter = 'TOPO.mat$';    
TOPOmat.num     = [1 Inf];     % Number of inputs required 
TOPOmat.help    = {'Select NIRS.mat for the subject(s).'}; % help text displayed


DelPreviousData         = cfg_menu; 
DelPreviousData.tag     = 'DelPreviousData';
DelPreviousData.name    = 'Delete previous .nir data file';
DelPreviousData.labels  = {'True','False'};
DelPreviousData.values  = {1,0};
DelPreviousData.val     = {0};
DelPreviousData.help    = {'Delete the previous data file.'}';
    
CreateNIRSCopy_false         = cfg_branch;  %to remove 
CreateNIRSCopy_false.tag     = 'CreateNIRSCopy_false';
CreateNIRSCopy_false.name    = 'Do not copy NIRS structure'; 
CreateNIRSCopy_false.help    = {'Do not copy NIRS structure.'
            'This will write over the previous NIRS.mat'}';

NewNIRSdir         = cfg_entry;             %to remove 
NewNIRSdir.name    = 'Directory for NIRS.mat';
NewNIRSdir.tag     = 'NewNIRSdir';       
NewNIRSdir.strtype = 's';
NewNIRSdir.val{1}  = 'Batch_00';
NewNIRSdir.num     = [1 Inf];     
NewNIRSdir.help    = {'Directory for NIRS.mat.'}'; 

CreateNIRSCopy         = cfg_branch;        %to remove 
CreateNIRSCopy.tag     = 'CreateNIRSCopy';
CreateNIRSCopy.name    = 'Create new directory and copy NIRS structure'; 
CreateNIRSCopy.val     = {NewNIRSdir};
CreateNIRSCopy.help    = {'Create new directory and copy NIRS structure there.'}';
        
%Common to most modules: for creating a new directory and copying NIRS.mat
% That the default option is false
NewDirCopyNIRS           = cfg_choice;      %to remove 
NewDirCopyNIRS.name      = 'Create new directory and copy NIRS.mat';
NewDirCopyNIRS.tag       = 'NewDirCopyNIRS';
NewDirCopyNIRS.values    = {CreateNIRSCopy_false CreateNIRSCopy}; 
NewDirCopyNIRS.val       = {CreateNIRSCopy_false}; 
NewDirCopyNIRS.help      = {'Choose whether to overwrite the NIRS.mat structure'
            'or to create a new directory'
            'and copy the NIRS.mat structure there'}'; 
        
%Common to most modules: for creating a new directory and copying NIRS.mat
%that the default option is true 
NewDirCopyNIRSTRUE           = cfg_choice; %to remove 
NewDirCopyNIRSTRUE.name      = 'Create new directory and copy NIRS.mat';
NewDirCopyNIRSTRUE.tag       = 'NewDirCopyNIRSTRUE';
NewDirCopyNIRSTRUE.values    = {CreateNIRSCopy_false CreateNIRSCopy}; 
NewDirCopyNIRSTRUE.val       = {CreateNIRSCopy}; 
NewDirCopyNIRSTRUE.help      = {'Choose whether to overwrite the NIRS.mat structure'
            'or to create a new directory'
            'and copy the NIRS.mat structure there'}'; 
        
        
%Common to ReadBoxy and Fast preprocessing   
cf_max         = cfg_entry;
cf_max.name    = 'Maximum coefficient of variation';
cf_max.tag     = 'cf_max';       
cf_max.strtype = 'r';
cf_max.num     = [1 Inf];
cf_max.val     = {0.1}; 
cf_max.help    = {'Reject channels with higher coefficient of variation (in normalized intensity).'};

amp_max         = cfg_entry;
amp_max.name    = 'Maximum amplitude';
amp_max.tag     = 'amp_max';       
amp_max.strtype = 'r';
amp_max.num     = [1 Inf];
amp_max.val     = {5}; 
amp_max.help    = {'Reject channels with higher amplitude (in normalized intensity) than entered threshold.'};


threshold_values         = cfg_branch;
threshold_values.tag     = 'threshold_values';
threshold_values.name    = 'Configuration options';
threshold_values.val     = {cf_max amp_max};
threshold_values.help    = {''}'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 1a configuration for ISS imagen system  (boxy mux 32)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputBOXY         = cfg_files;             % Select raw BOXY data files for this subject 
inputBOXY.name    = 'Select BOXY files';   % The displayed name
inputBOXY.tag     = 'inputBOXY';              %file names
inputBOXY.ufilter = '.0*';                 %files are labeled with extension .001, .002, .003, etc. 
inputBOXY.num     = [1 Inf];               % Number of inputs required 
inputBOXY.help    = {'Select raw BOXY data files for this subject record with and ISS imagent system (MUX32).'}; 

age1         = cfg_entry;
age1.name    = 'Subject age';
age1.tag     = 'age1';       
age1.strtype = 'r';
age1.num     = [1 Inf];
age1.val     = {25};
age1.help    = {'Age of the subject.',...
    'Used to caculate DPF to HbO/HbR conversion (Duncan et al 1995).',...
    'Only for 690,830,744,807 nm wavelenght.'};

raw_onset_files         = cfg_files;  
raw_onset_files.name    = 'Select onset files'; % The displayed name
raw_onset_files.tag     = 'raw_onset_files';          
raw_onset_files.num     = [0 Inf];     % Number of inputs required 
raw_onset_files.val{1}  = {''};
raw_onset_files.help    = {'Optional: Select raw onset files. '
    'Can be added at a later stage.'
    'Must specify one file for each data file, in same order.'}'; % help text displayed


subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {inputBOXY age1};
subj.help    = {'This module allows multi-subject processing, '
    'generating a NIRS.mat file for each subject. '
    'Note that a list of links to the NIRS.mat structures will be '
    'available as a virtual output for further processing'}';

generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Subject';
generic.help    = {'Note simple data organization that is recommended: '
    '1- directories with files for each subject must be in the same root directory. '
    '2- .prj file must have same name as expname and in folder \mtg\ '
    '3- anatomic files must be in folder \T1\'}';
generic.values  = {subj};
generic.num     = [1 Inf];


%path structure
output_path         = cfg_entry; %path
output_path.name    = 'Analysed files directory';
output_path.tag     = 'output_path';       
output_path.strtype = 's';
output_path.num     = [1 Inf];     
output_path.val     = {'c:/Data/Analysed/C01'}; 
output_path.help    = {'Path for .nir output files.',...
    'Should be something like ..\dataSPM\ (omit backslashes).'}; 

config_path         = cfg_branch;
config_path.tag     = 'config_path';
config_path.name    = 'Path Configuration options';
config_path.val     = {prjfile output_path}; 
config_path.help    = {''};

%Light wavelengths
Lambda         = cfg_entry; %Lambda
Lambda.name    = 'Laser wavelengths';
Lambda.tag     = 'Lambda';       
Lambda.strtype = 'r';
Lambda.num     = [1 Inf];     
Lambda.val     = {[830 690]};
Lambda.help    = {'Near Infrared laser wavelengths. Note order is critical and'
    'must correspond to recording order in raw data files.'}'; 
    
%Input frequency
i_freqISS         = cfg_entry; 
i_freqISS.name    = 'Input frequency';
i_freqISS.tag     = 'freq';       
i_freqISS.strtype = 'r';   
i_freqISS.num     = [1 Inf];     
i_freqISS.val     = {19.5312}; 
i_freqISS.help    = {'Input data frequency in Hertz.'}; 

%Minimum distance
distmin         = cfg_entry; 
distmin.name    = 'Minimum distance'; 
distmin.tag     = 'distmin';    
distmin.strtype = 'r';       
distmin.num     = [1 Inf];    
distmin.val     = {1}; 
distmin.help    = {'Cutoff: Minimum Cartesian channel distance in centimeters.'}; 

%Maximum distance
distmax         = cfg_entry; 
distmax.name    = 'Maximum distance'; 
distmax.tag     = 'distmax';      
distmax.strtype = 'r';       
distmax.num     = [1 1];     
distmax.val     = {6};
distmax.help    = {'Cutoff: Maximum Cartesian channel distance in centimeters.'}; 

%Number of MUX 
nb_Mux         = cfg_entry; %nb_Mux
nb_Mux.name    = 'MUX number';
nb_Mux.tag     = 'nb_Mux';       
nb_Mux.strtype = 'r';
nb_Mux.num     = [1 1];     
nb_Mux.val     = {32}; 
nb_Mux.help    = {'Number of MUX of ISS system it change the illumination order of the source, it is set and test in mux 32'}; 

%MaxSources
MaxSources         = cfg_entry; %MaxSources
MaxSources.name    = 'Maximum number of sources';
MaxSources.tag     = 'MaxSources';       
MaxSources.strtype = 'r';
MaxSources.num     = [1 1];     
MaxSources.val     = {64}; 
MaxSources.help    = {'Maximum number of sources'}; 

%Maximum number of detectors
nb_Det         = cfg_entry; %nb_Det
nb_Det.name    = 'Maximum number of detectors';
nb_Det.tag     = 'nb_Det';       
nb_Det.strtype = 'r';
nb_Det.num     = [1 1];     
nb_Det.val     = {16}; 
nb_Det.help    = {'Maximum number of detectors'}; 

%Maximum number of electrode
MaxElectrodes         = cfg_entry; %MaxElectrodes
MaxElectrodes.name    = 'Maximum number of electrodes';
MaxElectrodes.tag     = 'MaxElectrodes';       
MaxElectrodes.strtype = 'r';
MaxElectrodes.num     = [1 1];     
MaxElectrodes.val     = {16};
MaxElectrodes.help    = {'Maximum and habitual number of electrodes. This is only used'
    'in the interpolation, when one needs to fill in for missing electrode positions'}'; 

%use10_10system
use10_10system          = cfg_menu;
use10_10system.tag      = 'use10_10system';
use10_10system.name     = 'Use 10 10 system';
use10_10system.labels   = {'True','False: use 10 20 system'};
use10_10system.values   = {1,0};
use10_10system.val      = {1};
use10_10system.help     = {'10 10 system allows more precise location of optodes'
        'for viewing in .nir files.' }';

%Downsampling factor
resample         = cfg_entry; %resample
resample.name    = 'Downsampling factor';
resample.tag     = 'resample';       
resample.strtype = 'r';
resample.num     = [1 1];     
resample.val     = {1}; 
resample.help    = {'Downsampling factor: use 1 to keep all data.'
        'Note that no filtering will be done at this stage;'
        'therefore beware of aliasing artefacts.'}'; 

save_aux        = cfg_menu;
save_aux.tag    = 'save_aux';
save_aux.name   = 'Save aux channel';
save_aux.labels = {'True','False'};
save_aux.values = {1,0};
save_aux.val    =  {1};
save_aux.help   = {'Write auxiliary channel (AUX1,AUX2,AUX3,AUX4) in the NIRS.mat structure, if you use them a the recording'}';

cf1         = cfg_branch;
cf1.tag     = 'cf1';
cf1.name    = 'Configuration options';
cf1.val     = {Lambda i_freqISS distmin distmax save_aux...
    nb_Mux MaxSources nb_Det use10_10system MaxElectrodes resample};%  
cf1.help    = {'Configuration used for all subjects to be preprocessed.'};

boxyfastpre         = cfg_menu;
boxyfastpre.tag     = 'boxyfastpre';
boxyfastpre.name    = 'Allow fast preprocessing';
boxyfastpre.labels  = {'Yes','No'};
boxyfastpre.values  = {1,0};
boxyfastpre.val     = {1};
boxyfastpre.help    = {'Remove data series with too high coefficient of variation and amplitude.'}';

noisy           = cfg_menu;
noisy.tag       = 'noisy';
noisy.name      = 'Use the correction';
noisy.labels    = {'Yes','No'};
noisy.values    = {1,0};
noisy.val       = {1};
noisy.help      = {'Remove channels with high coefficient of variation.'}';


fastpreprocess         = cfg_branch;
fastpreprocess.tag     = 'fastpreprocess';
fastpreprocess.name    = 'STD deviation criteria';
fastpreprocess.val     = {noisy,cf_max};
fastpreprocess.help    = {'Criteria to reject noisy channels based on the high standard deviation of the normalised intensity.'};

STD_enable          = cfg_menu;
STD_enable.tag      = 'STD_enable';
STD_enable.name     = 'Use the correction';
STD_enable.labels   = {'No','Yes'};
STD_enable.values   = {0,1};
STD_enable.val      =  {0};
STD_enable.help     = {'Remove channels with high std variation.'}';

STD_amp         = cfg_entry;
STD_amp.tag     = 'STD_amp';
STD_amp.name    = 'STD amplitude';
STD_amp.strtype = 'r';       
STD_amp.num     = [1 1];     
STD_amp.val     = {0.1};
STD_amp.help    = {'Reject channels with higher std variation higher than STD amplitude (in normalized intensity).'};


STD_menu        = cfg_menu;
STD_menu.tag    = 'STD_menu';
STD_menu.name   = 'Apply option';
STD_menu.labels = {'Segment','Whole block'};
STD_menu.values = {1,2}; 
STD_menu.val    = {1};
STD_menu.help   = {'Look on different segment std deviation or whole block, segment will keep the channel if only a small part are lost.'}';

STD_amp_choice         = cfg_branch;
STD_amp_choice.tag     = 'STD_amp_choice';
STD_amp_choice.name    = 'STD deviation criteria';
STD_amp_choice.val     = {STD_enable STD_amp,STD_menu};
STD_amp_choice.help    = {'Criteria to reject noisy channels based on the high standard deviation of the normalised intensity.'};

DC_enable           = cfg_menu;
DC_enable.tag       = 'DC_enable';
DC_enable.name      = 'Use the correction';
DC_enable.labels    = {'No','Yes'};
DC_enable.values    = {0,1};
DC_enable.val       =  {0};
DC_enable.help      = {'Remove channels with low DC intensity.'};

DC_amp          = cfg_entry;
DC_amp.tag      = 'DC_amp';
DC_amp.name     = 'DC amplitude';
DC_amp.strtype  = 'r';       
DC_amp.num      = [1 1];     
DC_amp.val      = {100};
DC_amp.help     = {'Remove channel with DC intensity lower than a DC amplitude value.'}';

DC_amp_choice         = cfg_branch;
DC_amp_choice.tag     = 'DC_amp_choice';
DC_amp_choice.name    = 'DC intensity criteria';
DC_amp_choice.val     = {DC_enable DC_amp};
DC_amp_choice.help    = {'Criteria to reject noisy channels based on DC amplitude.'};



% Executable Branch
boxy1      = cfg_exbranch;      
boxy1.name = 'Read Boxy ISS';            
boxy1.tag  = 'boxy1'; 
boxy1.val  = {inputBOXY age1 prjfile output_path STD_amp_choice DC_amp_choice distmin  distmax };   
%boxy1.val  = {generic config_path cf1 STD_amp_choice DC_amp_choice};   
boxy1.prog = @nirs_run_readBoxyISS;  
boxy1.vout = @nirs_cfg_vout_readBoxyISS;
boxy1.help = {'Select raw BOXY data files for this subject.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_readBoxyISS(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 1b Configuration for .nirs from HOMER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputrawhomer         = cfg_files; %Select raw BOXY data files for this subject 
inputrawhomer.name    = 'Select .nirs'; % The displayed name
inputrawhomer.tag     = 'inputrawhomer';    %file names
inputrawhomer.ufilter = '.nirs';    %
inputrawhomer.num     = [1 Inf];     % Number of inputs required 
inputrawhomer.help    = {'Select raw .nirs file (HOMER FORMAT).'}; 



m_inputrawhomer      = cfg_menu;
m_inputrawhomer.tag  = 'm_inputrawhomer';
m_inputrawhomer.name = 'Options';
m_inputrawhomer.labels = {'Each file','Merge file', 'Merge and detrend'};
m_inputrawhomer.values = {0,1, 2};
m_inputrawhomer.val  = {1};
m_inputrawhomer.help = {'Read the data as each block separeted or merge and detrend each segment.'}';




% Executable Branch
E_rawhomer      = cfg_exbranch;      
E_rawhomer.name = 'Read .nirs (HoMER)';            
E_rawhomer.tag  = 'E_rawhomer'; 
E_rawhomer.val  = {inputrawhomer,age1, prjfile,output_path ,m_inputrawhomer};   
E_rawhomer.prog = @nirs_run_readhomerfile;  
E_rawhomer.vout = @nirs_cfg_vout_readhomerfile;
E_rawhomer.help = {'Select raw .nirs data files. Matlab structure containing field: ''d'' (raw data sample x channel), ''ml''(channel used sources, detector, weight wavelength), ''t'' (Time vector), data, ml, S vector will be used as aux trigger 1.'};

function vout = nirs_cfg_vout_readhomerfile(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end




inputNIRxscout         = cfg_files; %Select raw BOXY data files for this subject 
inputNIRxscout.name    = 'Select NIRx NIRScout data'; % The displayed name
inputNIRxscout.tag     = 'inputNIRxscout';    %file names
inputNIRxscout.ufilter = '.hdr';    %
inputNIRxscout.num     = [1 Inf];     % Number of inputs required 
inputNIRxscout.help    = {'Select raw NIRSCOUT recorded data.'}; 


b_foldername_NIRxscout         = cfg_branch;        %to remove 
b_foldername_NIRxscout.tag     = 'b_foldername_NIRxscout';
b_foldername_NIRxscout.name    = 'Folder name'; 
b_foldername_NIRxscout.val     = {};
b_foldername_NIRxscout.help    = {'Use folder name as identification.'}';

e_manualname_NIRxscout         = cfg_entry;             %to remove 
e_manualname_NIRxscout.name    = 'Name';
e_manualname_NIRxscout.tag     = 'e_manualname_NIRxscout';       
e_manualname_NIRxscout.strtype = 's';
e_manualname_NIRxscout.num     = [0 Inf];     
e_manualname_NIRxscout.help    = {'Specify data file name, if several file separate the name by a comma.',...
    'As example: file01, file02.'}'; 


b_manualname_NIRxscout         = cfg_branch;        %to remove 
b_manualname_NIRxscout.tag     = 'b_manualname_NIRxscout';
b_manualname_NIRxscout.name    = 'Manual entry'; 
b_manualname_NIRxscout.val     = {e_manualname_NIRxscout};
b_manualname_NIRxscout.help    = {'Specify name of the data file.'}';


b_defaultname_NIRxscout         = cfg_branch;        %to remove 
b_defaultname_NIRxscout.tag     = 'b_defaultname_NIRxscout';
b_defaultname_NIRxscout.name    = 'Default'; 
b_defaultname_NIRxscout.val     = {};
b_defaultname_NIRxscout.help    = {'Use default recording name of the raw file.'}';
        
%Common to most modules: for creating a new directory and copying NIRS.mat
% That the default option is false
c_nameconvention_NIRxscout           = cfg_choice;      %to remove 
c_nameconvention_NIRxscout.name      = 'Define name convention';
c_nameconvention_NIRxscout.tag       = 'c_nameconvention_NIRxscout';
c_nameconvention_NIRxscout.values    = {b_defaultname_NIRxscout b_manualname_NIRxscout b_foldername_NIRxscout }; 
c_nameconvention_NIRxscout.val       = {b_defaultname_NIRxscout}; 
c_nameconvention_NIRxscout.help      = {'Choose the convention you want to choose.'}'; 
        

b_shortdistunavailable         = cfg_branch;
b_shortdistunavailable.tag     = 'b_shortdistunavailable';
b_shortdistunavailable.name    = 'No';
b_shortdistunavailable.val     = {};
b_shortdistunavailable.help    = {'No short distance are recorded'}';

e_shortdistancedet         = cfg_entry;
e_shortdistancedet.name    = 'Detector associate to short distance probe';
e_shortdistancedet.tag     = 'e_shortdistancedet';       
e_shortdistancedet.strtype = 'r';
e_shortdistancedet.num     = [1 Inf];
e_shortdistancedet.val     = {[24]}; 
e_shortdistancedet.help    = {'Identify the short distance detector.'};

e_shortdistancesrs         = cfg_entry;
e_shortdistancesrs.name    = 'Source associate to short distance probe';
e_shortdistancesrs.tag     = 'e_shortdistancesrs';       
e_shortdistancesrs.strtype = 'r';
e_shortdistancesrs.num     = [1 Inf];
e_shortdistancesrs.val     = {[1,2,3,4,5,6,7,8]}; 
e_shortdistancesrs.help    = {'Identify the short distance localization on the helmet, as an example if the probes are placed in source 1 and source 2, write 1,2.'};

b_shortdistanceavailable         = cfg_branch;
b_shortdistanceavailable.tag     = 'b_shortdistanceavailable';
b_shortdistanceavailable.name    = 'Yes';
b_shortdistanceavailable.val     = {e_shortdistancesrs e_shortdistancedet};
b_shortdistanceavailable.help    = {'Short distance are available in the recording (det 16 will be associate).'}';

c_shortdistance         = cfg_choice;
c_shortdistance.tag     = 'c_shortdistance';
c_shortdistance.name    = 'Short distance probe';
c_shortdistance.values  = {b_shortdistanceavailable, b_shortdistunavailable};
c_shortdistance.val     = {b_shortdistanceavailable};
c_shortdistance.help    = {'NIRx provides a short distance probe with multiple short distances placed in the source location. If you use the probe to indicate the source number, you insert the probe. A zone file will be created to associate the short distance channels to the closest channels. They could be used as a regressor in the GLM. Use the zone as a regressor in the Component extract list.'};

% Executable Branch
E_readNIRxscout      = cfg_exbranch;      
E_readNIRxscout.name = 'Read NIRSCOUT';            
E_readNIRxscout.tag  = 'E_readNIRxscout'; 
E_readNIRxscout.val  = {inputNIRxscout,age1,prjfile,output_path,c_nameconvention_NIRxscout,STD_amp_choice, DC_amp_choice,distmin, distmax,c_shortdistance};   
E_readNIRxscout.prog = @nirs_run_readNIRxscout;  
E_readNIRxscout.vout = @nirs_cfg_vout_readNIRxscout;
E_readNIRxscout.help = {'Select NIRSCOUT raw data from NIRx.'};

function vout = nirs_cfg_vout_readNIRxscout(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

% Create a file selector component to select .snirf file.
inputSNIRF         = cfg_files;
inputSNIRF.name    = 'Select .snirf data'; % The displayed name
inputSNIRF.tag     = 'inputSNIRF';    %file names
inputSNIRF.ufilter = '.snirf';    %
inputSNIRF.num     = [1 Inf];     % Number of inputs required 
inputSNIRF.help    = {'Select .snirf file record data.'}; 

% Folder selector.
output_path_SNIRF         = cfg_entry; %path
output_path_SNIRF.name    = 'Path for output folder';
output_path_SNIRF.tag     = 'output_path_SNIRF';
output_path_SNIRF.strtype = 's';
output_path_SNIRF.num     = [1 Inf];
output_path_SNIRF.val     = {nirOutputDefaultPath};
output_path_SNIRF.help    = {'Path for the output folder where the .nir, .mat, .prj, etc. will be created'};

% Create the two branch required.
b_importProject         = cfg_branch; % Contains the project file selector.
b_importProject.tag     = 'b_importProject';
b_importProject.name    = 'Import .prj';
b_importProject.val     = {prjfile};
b_importProject.help    = {'Import an already existing project file'};

b_createProject         = cfg_branch; % Is an empty branch.
b_createProject.tag     = 'b_createProject';
b_createProject.name    = 'Create .prj';
b_createProject.val     = {};
b_createProject.help    = {'Create a new project file from the SNIRF file.'};

% Create an option on wether a project file should be created or imported.
c_createImportProjectSnirf          = cfg_choice;
c_createImportProjectSnirf.tag      = 'c_createImportProjectSnirf';
c_createImportProjectSnirf.name     = 'Create or import a project file';
c_createImportProjectSnirf.values   = {b_importProject, b_createProject};
c_createImportProjectSnirf.val      = {b_importProject};
c_createImportProjectSnirf.help     = {'Choose whether you want to import a project or create a new one'}';

% Executable Branch -- Create the drop down menu branch necessary to read .snirf file.
E_readSNIRF      = cfg_exbranch;      
E_readSNIRF.name = 'Read SNIRF';            
E_readSNIRF.tag  = 'E_readSNIRF'; 
E_readSNIRF.val  = {inputSNIRF,c_createImportProjectSnirf,output_path_SNIRF,STD_amp_choice, DC_amp_choice,distmin, distmax, c_shortdistance};
E_readSNIRF.prog = @nirs_run_readSNIRF; % Handle to the SNIRF import function.
E_readSNIRF.vout = @nirs_cfg_vout_readSNIRF;
E_readSNIRF.help = {'Import the SNIRF data into SPM.'};

function vout = nirs_cfg_vout_readSNIRF(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

% Create the two branch required.
b_importProject         = cfg_branch; % Contains the project file selector.
b_importProject.tag     = 'b_importProject';
b_importProject.name    = 'Import .prj';
b_importProject.val     = {prjfile};
b_importProject.help    = {'Import a project file that already exists'};

% Folder selector.
e_MultimodalPath          = cfg_files; %path
e_MultimodalPath.name     = 'Select the new folder';
e_MultimodalPath.tag      = 'e_MultimodalPath';
e_MultimodalPath.filter   = {'dir'};
e_MultimodalPath.ufilter  = '.*';    %
e_MultimodalPath.num      = [1 1];     % Number of inputs required 
e_MultimodalPath.help     = {'New folder used for multimodal files such as EEG, AUX, Video or Audio. If the location of the multimodal files has changed, enter the new directory location. This function expects all the multimodal files in the same folder. For a more specific adjustment use the Display GUI menu Setting/Multimodal files to define a new location.'};

b_MultimodalPath_yes        = cfg_branch; % Is an empty branch.
b_MultimodalPath_yes.tag    = 'b_MultimodalPath_yes';
b_MultimodalPath_yes.name   = 'New folder';
b_MultimodalPath_yes.val    = {e_MultimodalPath};
b_MultimodalPath_yes.help   = {'Folder adjustment for multimodal data.'};

b_MultimodalPath_no         = cfg_branch; % Is an empty branch.
b_MultimodalPath_no.tag     = 'b_MultimodalPath_no';
b_MultimodalPath_no.name    = 'No adjustement';
b_MultimodalPath_no.val     = {};
b_MultimodalPath_no.help    = {'No path adjustment for multimodal data.'};

% Choice to edit multimodal path
c_MultimodalPath          = cfg_choice;
c_MultimodalPath.tag      = 'c_MultimodalPath';
c_MultimodalPath.name     = 'Multimodal';
c_MultimodalPath.values   = {b_MultimodalPath_no, b_MultimodalPath_yes};
c_MultimodalPath.val      = {b_MultimodalPath_no};
c_MultimodalPath.help     = {''}';


E_NIRSmatdiradjust      = cfg_exbranch;      
E_NIRSmatdiradjust.name = 'Folder adjustment';            
E_NIRSmatdiradjust.tag  = 'E_NIRSmatdiradjust'; 
E_NIRSmatdiradjust.val  = {NIRSmat,c_MultimodalPath};   
E_NIRSmatdiradjust.prog = @nirs_run_NIRSmatdiradjust;  
E_NIRSmatdiradjust.vout = @nirs_cfg_vout_readNIRxscout;
E_NIRSmatdiradjust.help = {'NIRS.mat structure contains the information about all the localizations of the data files on the hard drive.',...
    'This function will modify all the subdirectories of the data files at the actuel localization  of the NIRS.mat.',...
    'Be careful, it is not yet tested for path that the dir contain the same name many times.'};

function vout = nirs_cfg_vout_NIRSmatdiradjust(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


e_NIRSmatdirnewbranch         = cfg_entry; %subfolder new branch
e_NIRSmatdirnewbranch.name    = 'Create a new branch of analysis in the subfolder';
e_NIRSmatdirnewbranch.tag     = 'e_NIRSmatdirnewbranch';       
e_NIRSmatdirnewbranch.strtype = 's';
e_NIRSmatdirnewbranch.num     = [1 Inf];     
e_NIRSmatdirnewbranch.val     = {'Branch'}; 
e_NIRSmatdirnewbranch.help    = {'A copy of the last step will be create to create a new branch of analysis.'}; 


m_newbranchcomponent         = cfg_menu; 
m_newbranchcomponent.name    = 'Components and corrections';
m_newbranchcomponent.tag     = 'm_newbranchcomponent';      
m_newbranchcomponent.labels  = {'Keep', 'Clear' };
m_newbranchcomponent.values  = {1,0};
m_newbranchcomponent.val     = {1}; 
m_newbranchcomponent.help    = {'When you create a new branch, you could keep the component from the previous operation tree or clear them to restart in the new branch.'};

E_NIRSmatcreatenewbranch      = cfg_exbranch;      
E_NIRSmatcreatenewbranch.name = 'New branch';            
E_NIRSmatcreatenewbranch.tag  = 'E_NIRSmatcreatenewbranch'; 
E_NIRSmatcreatenewbranch.val  = {NIRSmat,e_NIRSmatdirnewbranch,m_newbranchcomponent};   
E_NIRSmatcreatenewbranch.prog = @nirs_run_NIRSmatcreatenewbranch;  
E_NIRSmatcreatenewbranch.vout = @nirs_cfg_vout_NIRSmatcreatenewbranch;
E_NIRSmatcreatenewbranch.help = {'NIRS.mat structure contains the information about all the localizations of the data files on the hard drive.',...
    'This function will create a copy of the NIRS.mat structure to create a new independant branch of analysis.'};

function vout = nirs_cfg_vout_NIRSmatcreatenewbranch(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 2a DataOnset Create manual trig using ISS aux input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



E_aux2manualname         = cfg_entry;
E_aux2manualname.name    = 'File output';
E_aux2manualname.tag     = 'E_aux2manualname';        
E_aux2manualname.strtype = 's';
E_aux2manualname.num     = [1 inf];
E_aux2manualname.val     = {'Manualtrig'};
E_aux2manualname.help    = {'OutputName'};

E_aux2manualoption         = cfg_menu; 
E_aux2manualoption.name    = 'Option';
E_aux2manualoption.tag     = 'E_aux2manualoption';      
E_aux2manualoption.labels  = {'Separate file'};
E_aux2manualoption.values  = {1};
E_aux2manualoption.val     = {1}; 
E_aux2manualoption.help    = {'Separate file : identify each trig as a separe file.',...
    'Concatenate file : Special use if you plan to concatenate yout file in nirs file',...
    'this option will adjust the trig time assuming that each file from the project are place,'...
    'one after the another'};

%Create trigger for avg from GLMonset 
E_aux2manualtrig      = cfg_exbranch;
E_aux2manualtrig.name = 'AuxTrig to ManualTrig';
E_aux2manualtrig.tag  = 'E_aux2manualtrig';
E_aux2manualtrig.val  = {NIRSmat, E_aux2manualname};
E_aux2manualtrig.prog = @nirs_run_E_aux2manualtrig;
E_aux2manualtrig.vout = @nirs_cfg_vout_E_aux2manualtrig;
E_aux2manualtrig.help = {'Use ISS aux trig to create editable manual trig file to be use with manual trig',...
    'file will be create will be create in the current directory',...
    'under the name (fileid_trigid_trig.m)' };

function vout = nirs_cfg_vout_E_aux2manualtrig(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 2b DataOnset insert manual trig definition in the data E_manualtrig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_trigmode          = cfg_menu;
m_trigmode.tag      = 'm_trigmode';
m_trigmode.name     = 'Option';
m_trigmode.labels   = {'Replace trig', 'Add trig'};
m_trigmode.values   = {0,1};
m_trigmode.val      = {0};
m_trigmode.help     = {'Replace trig will erase the default trigger present during the recording to place only the trig you have define manually.',...
                   'Add trig will add new trig to the previous trig already present in the file definition.'};

m_trigfile         = cfg_files;
m_trigfile.name    = 'Enter new trig definition'; 
m_trigfile.tag     = 'm_trigfile';       %file names
m_trigfile.filter  = 'm';
m_trigfile.ufilter = 'trig.m$';    
m_trigfile.num     = [1 Inf];     % Number of inputs required 
m_trigfile.help    = {' Open ManualTrig file. The manual trig is a .m Matlab script that defines the new TrigValue=X',...
    'an integer used for the trigger identification and the timingfile{1}=[15,30]; which indicate the timing of the'...
    'trig event in seconds for the file 1. As the NIRS.mat structure could contain one or many files you must define a timingfile for each of the files present in your data. If as example the file 2 does not have any trig events defined the timingfile as empty value: timingfile{2} = [ ];.  Be aware that .m files names do not support space or accent in the name definition. ',...
    'Example: ',...
    'TrigValue = 3;',...
    'timingfile{1} = [15, 30, 60];',...
    'timingfile{2} = [ ];',...
    'timingfile{3} = [15, 30, 60];',...
    'timingfile{4} = [];'}; 


E_manualtrig      = cfg_exbranch;
E_manualtrig.name = 'ManualTrig to AuxTrig';
E_manualtrig.tag  = 'E_manualtrig';
E_manualtrig.val  = {NIRSmat m_trigfile m_trigmode};
E_manualtrig.prog = @nirs_run_E_manualtrig;
E_manualtrig.vout = @nirs_cfg_vout_E_manualtrig;
E_manualtrig.help = {'Use the manually edited file to enter trigger information. Triggers are used to identify the onset of the events to be segmented, normalized, model HRF or average the data. They are almost essential to data analysis.'};

function vout = nirs_cfg_vout_E_manualtrig(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

e_Concatenate_blocid            = cfg_entry;
e_Concatenate_blocid.tag        = 'e_Concatenate_blocid';
e_Concatenate_blocid.name       = 'File number (set to 0 take all)';
e_Concatenate_blocid.strtype    = 'r';
e_Concatenate_blocid.num        = [1 Inf];
e_Concatenate_blocid.val        = {[0]};
e_Concatenate_blocid.help       = {'Enter the file number to group as example 1,2,3; in the display file list will be identify such as Bloc001, set 0 to take all'};


m_Concatenate_option        = cfg_menu;
m_Concatenate_option.tag    = 'm_Concatenate_option';
m_Concatenate_option.name   = 'Options';
m_Concatenate_option.labels = {'Merge only','Merge and Detrend'};
m_Concatenate_option.values = {0,1};
m_Concatenate_option.val    = {0};
m_Concatenate_option.help   = {'Merge data off all nirs and block selected, merge only do not correct any offset while merge and detrend do a detrending operation on each block before to merge them'};


E_Concatenate_file      = cfg_exbranch;
E_Concatenate_file.name = 'Concatenate blocks in NIRS.mat';
E_Concatenate_file.tag  = 'E_Concatenate_file';
E_Concatenate_file.val  = {NIRSmat,m_Concatenate_option,e_Concatenate_blocid};
E_Concatenate_file.prog = @nirs_run_E_Concatenate_file;
E_Concatenate_file.vout = @nirs_cfg_vout_E_Concatenate_file;
E_Concatenate_file.help = {'Join several blocks in the NIRS.mat, support EEG or AUX multimodal information, additional file Allxxx.dat will be created. Do not support video'};

function vout = nirs_cfg_vout_E_Concatenate_file(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

f_nirsmatinfo         = cfg_files;
f_nirsmatinfo.name    = 'Enter NIRS.mat list to group (xls)'; 
f_nirsmatinfo.tag     = 'f_nirsmatinfo';       %file names
f_nirsmatinfo.filter  = {'xlsx','xls','txt'};
f_nirsmatinfo.ufilter = '.*';    
f_nirsmatinfo.num     = [1 Inf];     % Number of inputs required 
f_nirsmatinfo.help    = {'Create the list of subject to group in a new forder example: ./GrandAverage/List.xlsx. The excel file must have 2 columns, first, the NIRS.mat subject location and second the channel list order to respect or zone to use. The result will kept the same helmet of the first subject in the list.'}; 

m_Concatenate_Exclude           = cfg_menu;
m_Concatenate_Exclude.tag       = 'm_Concatenate_Exclude';
m_Concatenate_Exclude.name      = 'Set exclude channel to NaN';
m_Concatenate_Exclude.labels    = {'Exclude','Keep'};
m_Concatenate_Exclude.values    = {0,1};
m_Concatenate_Exclude.val       = {0};
m_Concatenate_Exclude.help      = {'Set exclude channel to NaN,',...
    'Exclude: rejected channel will be set to NaN. ',...
    'Keep: no special process concerning the rejected channels will be applied.'};

m_Concatenate_Normalized        = cfg_menu;
m_Concatenate_Normalized.tag    = 'm_Concatenate_Normalized';
m_Concatenate_Normalized.name   = 'Normalization';
m_Concatenate_Normalized.labels = {'No normalization','Min-Max Normalization','z-score Normalization'};
m_Concatenate_Normalized.values = {0,1,2};
m_Concatenate_Normalized.val    = {0};
m_Concatenate_Normalized.help   = {'Apply a normalization to reduce individual subject  variability. Only apply using the channel list',...
    'No normalization : do not apply any normalization. ',...
    'Min-Max Normalization: apply min max normalization fixing boundary between 0 and 1.',... 
    'Z-score Normalization: apply z-score normalization.',... 
    'Moeller, J., 2015. A word on standardization in longitudinal studies: don’t. Front Psychol 6. https://doi.org/10.3389/fpsyg.2015.01389'};

E_Concatenate_nirsmat      = cfg_exbranch;
E_Concatenate_nirsmat.name = 'Concatenate NIRS.mat (multi subject and zone)';
E_Concatenate_nirsmat.tag  = 'E_Concatenate_nirsmat';
E_Concatenate_nirsmat.val  = {f_nirsmatinfo,m_Concatenate_option, m_Concatenate_Exclude, m_Concatenate_Normalized};
E_Concatenate_nirsmat.prog = @nirs_run_E_Concatenate_nirsmat;
E_Concatenate_nirsmat.vout = @nirs_cfg_vout_E_Concatenate_nirsmat;
E_Concatenate_nirsmat.help = {'Join several NIRS.mat one after could be use to do group average, individual multimodal (EEG, AUX, Video, AUDIO) are not supported for this function '};

function vout = nirs_cfg_vout_E_Concatenate_nirsmat(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 3a Preprocessing Normalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_choiceNan         =  cfg_menu;
m_choiceNan.tag     =  'm_choiceNan';
m_choiceNan.name    = 'Exclude artifacts from the Io calculation';       
m_choiceNan.labels  = {'Yes','No'};
m_choiceNan.values  = {1,0};
m_choiceNan.val     = {1};
m_choiceNan.help    = {'Io will exclude artifacts (yellow marking).'};



b_choiceglobal         = cfg_branch;
b_choiceglobal.tag     = 'b_choiceglobal';
b_choiceglobal.name    = 'Normalization by files';
b_choiceglobal.val     = {m_choiceNan};
b_choiceglobal.help    = {'Define Io as each file average.',...
    'Apply on the whole file.'}';


timedurationinternan         = cfg_entry;
timedurationinternan.tag     = 'timedurationinternan';
timedurationinternan.name    = 'Selection of the periods between artifact identification (yellow marker)';       
timedurationinternan.strtype = 's';
timedurationinternan.val     = {'auto'}; 
timedurationinternan.num     = [1 Inf];
timedurationinternan.help    = {'Normalization will be done between each artifact periode'};

b_choiceinternan         = cfg_branch;
b_choiceinternan.tag     = 'b_choiceinternan';
b_choiceinternan.name    = 'Normalization inter-artifacts';
b_choiceinternan.val     = {timedurationinternan};
b_choiceinternan.help    = {'Normalize each segment by its mean, segment are separated by artifacted periods.'}';

trigger         = cfg_entry;
trigger.name    = 'Trigger';
trigger.tag     = 'trigger';
trigger.strtype = 'r';
trigger.num     = [1 Inf];
trigger.val     = {1}; 
trigger.help    = {'Enter trigger number.'};

pretime         = cfg_entry;
pretime.name    = 'PreTime';
pretime.tag     = 'pretime';       
pretime.strtype = 's';
pretime.num     = [1 Inf];
pretime.val     = {'5'}; 
pretime.help    = {'Time to include before the trigger.',...
    'Write keyword ''start'' to go to the beginning of the raw segment.',...
    'Define as a positive value, as example 5 seconds before the onset. Unintuitively, -30 will give you 30 seconds after the trigger.'};

posttime         = cfg_entry;
posttime.name    = 'PostTime';
posttime.tag     = 'posttime';       
posttime.strtype = 's';
posttime.num     = [1 Inf];
posttime.val     = {'30'};
posttime.help    = {'Time to include after the trigger.',...
    'Write keyword ''end'' to go to the end of the raw segment.'};


m_NormType          = cfg_menu;
m_NormType.tag      = 'm_NormType';
m_NormType.name     = 'How define Io';
m_NormType.labels   = {'Io = Pretime to 0','Io = Pretime to PostTime'};%,'Substract I-Io (Pretime','Do nothing, Raw segmented'};
m_NormType.values   = {0,1};%,2,3
m_NormType.val      = {0};
m_NormType.help     = {'Choose one of the two definition of Io'}';


b_choicenormstim         = cfg_branch;
b_choicenormstim.tag     = 'b_choicenormstim';
b_choicenormstim.name    = 'Normalization around trigger';
b_choicenormstim.val     = {trigger pretime posttime m_NormType m_choiceNan};
b_choicenormstim.help    = {'Normalization using pre post time trigger.'}';

c_normtype          = cfg_choice;
c_normtype.tag      = 'normtype';
c_normtype.name     = 'Normalization type';
c_normtype.values   = {b_choiceglobal,b_choicenormstim,b_choiceinternan};
c_normtype.val      = {b_choiceglobal};
c_normtype.help     = {''}';

% Executable Branch
E_normalization      = cfg_exbranch;
E_normalization.name = 'Normalization dOD=log(I/Io)';
E_normalization.tag  = 'normalization';
E_normalization.val  = {NIRSmat DelPreviousData c_normtype};
E_normalization.prog = @nirs_run_normalize;
E_normalization.vout = @nirs_cfg_vout_normalize;
E_normalization.help = {'Normalize raw data, transfer optical intensity in delta optical density',...
    'dOD  = log(I/Io)'};
%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_normalize(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



% Executable Branch
E_segment      = cfg_exbranch;
E_segment.name = 'Segment';
E_segment.tag  = 'segment';
E_segment.val  = {NIRSmat DelPreviousData trigger pretime posttime};
E_segment.prog = @nirs_run_segment;
E_segment.vout = @nirs_cfg_vout_segment;
E_segment.help = {'Segment raw data. Egal segment pretime posttime around trigger'};
%make NIRS.mat available as a dependency

function vout = nirs_cfg_vout_segment(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 4a Artifact detection Step Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PrintReport         = cfg_menu;
PrintReport.tag     = 'PrintReport';
PrintReport.name    = 'Save general report';
PrintReport.labels  = {'Yes','No'};
PrintReport.values  = {1,0};
PrintReport.val     = {1}; 
PrintReport.help    = {'Save an image of the percentage of rejected channels over time in the NIRS.mat folder.'}';

m_meandiff          = cfg_menu;
m_meandiff.tag      = 'm_meandiff';
m_meandiff.name     = 'Apply';
m_meandiff.labels   = {'Yes','No'};
m_meandiff.values   = {1,0};
m_meandiff.val      = {1}; 
m_meandiff.help     = {'Use the artifact detection using moving average.'}';

thresholdstep         = cfg_entry;
thresholdstep.name    = 'Z-score threshold';
thresholdstep.tag     = 'thresholdstep';
thresholdstep.strtype = 'r';
thresholdstep.num     = [1 1];
thresholdstep.val     = {3}; 
thresholdstep.help    = {'The difference (D) between M1 and M2 over time is converted into a standardized z-score distribution. A z-score threshold determines which variations in the moving average difference of two subsequent intervals are considered abnormal, thus are represent artifacts. Depending on the value defined by the user,  variations higher than this z-score identified.'};
   
too_small_step_dur         = cfg_entry;
too_small_step_dur.name    = 'Ignored interval shorter than x seconde:';
too_small_step_dur.tag     = 'too_small_step_dur';
too_small_step_dur.strtype = 'r';
too_small_step_dur.num     = [1 1];
too_small_step_dur.val     = {0.4}; 
too_small_step_dur.help    = {'This specifies the minimal duration (seconds) of an abnormal variation in the signal, as defined by the z-score threshold of the moving average, in order to be considered as an artifact. Intervals where the z-score threshold appears for a period shorter than this time are ignored.'};

movingaverage_nbpoint         = cfg_entry;
movingaverage_nbpoint.name    = 'Moving average step (seconde)';
movingaverage_nbpoint.tag     = 'movingaverage_nbpoint';
movingaverage_nbpoint.strtype = 'r';
movingaverage_nbpoint.num     = [1 1];
movingaverage_nbpoint.val     =  {1}; 
movingaverage_nbpoint.help    = {'Defines the duration (sec) of the time interval (step) to calculate moving average M1 and M2: Xn to Xn+step. Depending on the smoothing the user wants to apply, i.e. how sensitive detection of the signal’s variation should be, a higher or lower value should be indicated.'};


printreportthreshold        = cfg_menu;
printreportthreshold.tag    = 'printreportthreshold';
printreportthreshold.name   = 'Print threshold report for each channel';
printreportthreshold.labels = {'Yes','No'};
printreportthreshold.values = {1,0};
printreportthreshold.val    = {0}; 
printreportthreshold.help   = {'Select if you want to save a figure that displays the threshold detection made for each channel. The figure will be saved in …\ArtifactDetection_Report\EachCh folder at the NIRS.mat home location.'}';

b_meandiff         = cfg_branch;
b_meandiff.tag     = 'b_meandiff';
b_meandiff.name    = 'Artifact detection using moving average';
b_meandiff.val     = {m_meandiff thresholdstep  movingaverage_nbpoint too_small_step_dur printreportthreshold};
b_meandiff.help    = {'Moving average allows identifying discontinuity or strong perturbations in the signal. It is a criteria based on the signal’s mean variation for two subsequent time intervals: M1 from Xn to Xn+step and M2 from Xn+1 to Xn+step+1. The step is defined by the user and defines the period of time for the moving average time windows. Changes of light intensity between M1 and M2 are then specified as the difference (D) between both means: D = M2-M1. Difference D is transferred into z-scores normalized over the entire dataset to determine which periods have increased signal variations. The threshold, i.e. the z-score above which abnormal variations in the signal are identified as artifact intervals, is to be defined by the user. '}';

m_min_subinterval        = cfg_menu;
m_min_subinterval.tag    = 'm_min_subinterval';
m_min_subinterval.name   = 'Apply';
m_min_subinterval.labels = {'Yes','No'};
m_min_subinterval.values = {1,0};
m_min_subinterval.val    = {1};
m_min_subinterval.help   = {'Use minimum subinterval.'}';

min_subinterval         = cfg_entry;
min_subinterval.name    = 'Minimum subinterval duration (seconds)';
min_subinterval.tag     = 'min_subinterval';       
min_subinterval.strtype = 'r';
min_subinterval.num     = [1 1];
min_subinterval.val     = {2}; 
min_subinterval.help    = {'Select the minimum duration (in seconds) for a good subinterval.'};

b_min_subinterval         = cfg_branch;
b_min_subinterval.tag     = 'b_min_subinterval';
b_min_subinterval.name    = 'Minimal subinterval';
b_min_subinterval.val     = {m_min_subinterval min_subinterval};
b_min_subinterval.help    = {'This criterion helps to consider small artifacts as one event instead of two subparts. For each detected interval, ensure that more than 2 seconds separates it from the previous or next detected interval, otherwise consider them as one event.'};

m_corr          = cfg_menu;
m_corr.tag      = 'm_corr';
m_corr.name     = 'Apply';
m_corr.labels   = {'Yes','No'};
m_corr.values   = {1,0};
m_corr.val      = {1}; 
m_corr.help     = {'Use Correlation between channels for artifact interval.'}';

corr_thr         = cfg_entry;
corr_thr.name    = 'Correlation threshold';
corr_thr.tag     = 'corr_thr';       
corr_thr.strtype = 'r';
corr_thr.num     = [1 1];
corr_thr.val     = {0.80}; 
corr_thr.help    = {'Set the minimal correlation threshold between the channels. If the correlation is stronger then this threshold the time course will be marked as artifact.'};

b_corr         = cfg_branch;
b_corr.tag     = 'b_corr';
b_corr.name    = 'Correlation between channels for artifact interval';
b_corr.val     = {m_corr corr_thr};
b_corr.help    = {'Artifact detection is based on a correlation coefficient between channels. It determines the Pearson’s correlation coefficient threshold of channels to be considered as being affected by the same event. For each artifact interval that has previously been detected, channels that have a correlation equal or above the threshold with the time course of this artifact, are detected as well. This means that it asks for each artifact to find the channels showing the same signature, but have not been  detected based on the previous criteria.'}';

m_minpourcentagebad      = cfg_menu;
m_minpourcentagebad.tag  = 'm_minpourcentagebad';
m_minpourcentagebad.name = 'Apply';
m_minpourcentagebad.labels = {'Yes','No'};
m_minpourcentagebad.values = {1,0};
m_minpourcentagebad.val  = {1}; 
m_minpourcentagebad.help = {'Use Minimal percentage of bad channels to be marked as artifact.'}';

minpourcentagebad         = cfg_entry;
minpourcentagebad.name    = 'Minimal percentage of bad channels';
minpourcentagebad.tag     = 'minpourcentagebad';       
minpourcentagebad.strtype = 'r';
minpourcentagebad.num     = [1 1];
minpourcentagebad.val     = {5}; 
minpourcentagebad.help    = {'Define a minimal number of bad channels percentage to be marked at the same interval.',...
    'Max value 100, Min value 0'};

b_minpourcentagebad         = cfg_branch;
b_minpourcentagebad.tag     = 'b_minpourcentagebad';
b_minpourcentagebad.name    = 'Minimal percentage of bad channels to be mark as artifact';
b_minpourcentagebad.val     = {m_minpourcentagebad minpourcentagebad};
b_minpourcentagebad.help    = {'In some case, the signal is clean and very few channels are detected as artifact.',...
    'They are few chances that is due to a movement due to the impact on such a low number of channels.',...
    'You may choose to restore if less than 10% of the channels are detected as noisy for this period of time.'};


% Executable Branch
E_artefactdetection      = cfg_exbranch;
E_artefactdetection.name = 'Artifact Detection';
E_artefactdetection.tag  = 'E_artefactdetection';
E_artefactdetection.val  = {NIRSmat DelPreviousData  PrintReport b_meandiff b_minpourcentagebad b_min_subinterval b_corr};
E_artefactdetection.prog = @nirs_run_step_artefact_detection;
E_artefactdetection.vout = @nirs_run_vout_step_artefact_detection;
E_artefactdetection.help = {'Automatic artifact detection modules apply automatic artifact detection on either raw or normalized unfiltered data to detect abrupt variations, potentially related to an artifact. The user has the option to consider several criteria to improve detection.',...
    'These fourth criteria are:',...
    '‘Artifact detection using moving average’;',...
    '‘Minimal percentage of bad channels to be marked as artifact’;',...
    '‘Minimal subinterval’;',...
    '‘Correlation between channels for artifact’. ',...
    'A more detailed description of each criterion is presented above. '};




%make NIRS.mat available as a dependency
function vout = nirs_run_vout_step_artefact_detection(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end





%Input Frequency
i_Freq_cardiac         = cfg_entry; 
i_Freq_cardiac.name    = 'Frequency range to detect cardiac peak'; 
i_Freq_cardiac.tag     = 'i_Freq_cardiac';      
i_Freq_cardiac.strtype = 'r';       
i_Freq_cardiac.num     = [1 inf];     
i_Freq_cardiac.val     = {[0.8,2.3]};
i_Freq_cardiac.help    = {'Enter the range to detect the cardiac peak, the peak value between this interval. If you are outside of the range of the peak, please adjust these values; a normal cardiac pulse will be around 1 Hz for an adult population and 2 Hz for babies.'}; 

%Input Frequency
i_COHTRESHOLD_cardiac         = cfg_entry; 
i_COHTRESHOLD_cardiac.name    = 'Minimal coherence Cxy(f)'; 
i_COHTRESHOLD_cardiac.tag     = 'i_Freq_crossspectrum';      
i_COHTRESHOLD_cardiac.strtype = 'r';       
i_COHTRESHOLD_cardiac.num     = [1 inf];     
i_COHTRESHOLD_cardiac.val     = {[0.2]};
i_COHTRESHOLD_cardiac.help    = {'Threshold to determine if they are or not cardiac coherence.'}; 


%Input Frequency
i_minch_cardiac          = cfg_entry; 
i_minch_cardiac.name     = 'Minimal percentage of channel'; 
i_minch_cardiac.tag      = 'i_minch_cardiac';      
i_minch_cardiac.strtype  = 'r';       
i_minch_cardiac.num      = [1 inf];     
i_minch_cardiac.val      = {[10]};
i_minch_cardiac.help     = {'Reject the channel is less than 10% of the coherence among other channels obtain the minimal value.'}; 
i_cardiacwidth           = cfg_entry; 
i_cardiacwidth.name      = 'Range around the peak (Hz)'; 
i_cardiacwidth.tag       = 'i_cardiacwidth';      
i_cardiacwidth.strtype   = 'r';       
i_cardiacwidth.num       = [1 inf];     
i_cardiacwidth.val       = {[0]};
i_cardiacwidth.help      = {'Define the range to compute the coherence around the peak.'}; 



% Executable Branch
E_chcardiaccontrol      = cfg_exbranch;
E_chcardiaccontrol.name = 'Cardiac Detection';
E_chcardiaccontrol.tag  = 'E_chcardiaccontrol';
E_chcardiaccontrol.val  = {NIRSmat, i_Freq_cardiac i_COHTRESHOLD_cardiac i_minch_cardiac i_cardiacwidth };
E_chcardiaccontrol.prog = @nirs_run_chcardiaccontrol;
E_chcardiaccontrol.vout = @nirs_run_vout_chcardiaccontrol;
E_chcardiaccontrol.help = {'Detection of the cardiac beat using coherence measure among channels and the module rejects channels without any cardiac evidences '};




%make NIRS.mat available as a dependency
function vout = nirs_run_vout_chcardiaccontrol(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 3c Preprocessing Nullify bad intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paddingtime         = cfg_entry;
paddingtime.name    = 'Padding time';
paddingtime.tag     = 'paddingtime';
paddingtime.strtype = 'r';
paddingtime.num     = [1 1];
paddingtime.val     = {2};
paddingtime.help    = {'Time (in seconds) padded with NaN before and after the nullified intervals.'};

% Executable Branch
E_nullifybad      = cfg_exbranch;
E_nullifybad.name = 'Nullify bad intervals';
E_nullifybad.tag  = 'nullifybad';
E_nullifybad.val  = {NIRSmat DelPreviousData paddingtime};
E_nullifybad.prog = @nirs_run_E_nullifybad;
E_nullifybad.vout = @nirs_cfg_vout_E_nullifybad;
E_nullifybad.help = {'Nullify bad intervals (previously identified by Artifact detection modules).'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_E_nullifybad(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 3c Preprocessing Regress correlated bad intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 3d Preprocessing Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lowcutfreq         = cfg_entry;
lowcutfreq.name    = 'Low-Pass cutoff';
lowcutfreq.tag     = 'lowcutfreq';       
lowcutfreq.strtype = 's';
lowcutfreq.num     = [1 inf];
lowcutfreq.val     = {'0.2'};
lowcutfreq.help    = {'Cutoff of Low-Pass filter.',...
    'Enter a value if you want to apply a low pass filter, if you don''t, write ''No''.'};

highcutfreq         = cfg_entry;
highcutfreq.name    = 'High-Pass cutoff';
highcutfreq.tag     = 'highcutfreq';       
highcutfreq.strtype = 's';
highcutfreq.num     = [1 inf];
highcutfreq.val     = {'No'};
highcutfreq.help    = {'Cutoff for High-Pass filter.',...
    'Enter a value if you want to apply a high pass filterif you don''t, write ''No''.'};

paddingsymfilter         = cfg_menu;
paddingsymfilter.tag     = 'paddingsymfilter';
paddingsymfilter.name    = 'Symmetric Padding';
paddingsymfilter.help    = {'To avoid border effect on short segment filtering, apply symetrisation of the segment before filtering.'};
paddingsymfilter.labels  = {'True','False' }';                
paddingsymfilter.values  = {1 0};
paddingsymfilter.val     = {1};

interpolatebadfilter         = cfg_menu;
interpolatebadfilter.tag     = 'interpolatebadfilter';
interpolatebadfilter.name    = 'Interpolate bad interval';
interpolatebadfilter.help    = {'Bad interval or period in yellow are interpolate before filtering.'};
interpolatebadfilter.labels  = {'True','False' }';                
interpolatebadfilter.values  = {1 0};
interpolatebadfilter.val     = {0};

filterorder         = cfg_entry;
filterorder.name    = 'Filter order';
filterorder.tag     = 'filterorder';       
filterorder.strtype = 'r';
filterorder.num     = [1 Inf];
filterorder.val     = {4};
filterorder.help    = {'Butterworth filter order.'};

% Executable Branch
E_filter      = cfg_exbranch;
E_filter.name = 'Filter';
E_filter.tag  = 'bpfilt';
E_filter.val  = {NIRSmat DelPreviousData lowcutfreq highcutfreq filterorder paddingsymfilter interpolatebadfilter};
E_filter.prog = @nirs_run_filter;
E_filter.vout = @nirs_cfg_vout_filter;
E_filter.help = {'Bandpass filtering of the data. Bad intervals are interpolated prior to filtering.',...
    'Input & Output: Normalized intensity or concentrations.'};

function vout = nirs_cfg_vout_filter(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



% Executable Branch
E_detrend      = cfg_exbranch;
E_detrend.name = 'Detrend';
E_detrend.tag  = 'E_detrend';
E_detrend.val  = {NIRSmat DelPreviousData };
E_detrend.prog = @nirs_run_E_detrend;
E_detrend.vout = @nirs_cfg_vout_detrend;
E_detrend.help = {'Use detrend on each blocks.'};

function vout = nirs_cfg_vout_detrend(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 3e Preprocessing ODtoHbOHbR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PVF         = cfg_entry;
PVF.name    = 'Partial Volume Factors';
PVF.tag     = 'PVF';
PVF.strtype = 'r';
PVF.num     = [1 2];
PVF.val     = {[1 1]};
PVF.help    = {'Enter the partial volume factor values for each wavelength ',...
    'as a vector: [PVF(lambda_1) ... PVF(lambda_n)].'};

b_ODtoHbOHbR_DPF3         = cfg_entry;
b_ODtoHbOHbR_DPF3.name    = 'Manual';
b_ODtoHbOHbR_DPF3.tag     = 'b_ODtoHbOHbR_DPF3';
b_ODtoHbOHbR_DPF3.strtype = 'r';
b_ODtoHbOHbR_DPF3.num     = [1 2];
b_ODtoHbOHbR_DPF3.val     = {[6 6]};
b_ODtoHbOHbR_DPF3.help    = {'Adjust differential pathlength factor manually.'};

b_ODtoHbOHbR_DPF2         = cfg_branch;
b_ODtoHbOHbR_DPF2.name    = 'Duncan et al. 1996';
b_ODtoHbOHbR_DPF2.tag     = 'b_ODtoHbOHbR_DPF2';
b_ODtoHbOHbR_DPF2.val     = {};
b_ODtoHbOHbR_DPF2.help    = {'Adjust differential pathlenght factor in function of age.',...
    'DOI: 10.1203/00006450-199605000-00025'};



b_ODtoHbOHbR_DPF1         = cfg_branch;
b_ODtoHbOHbR_DPF1.name    = 'Scholkmann and Wolf 2013';
b_ODtoHbOHbR_DPF1.tag     = 'b_ODtoHbOHbR_DPF1';
b_ODtoHbOHbR_DPF1.val     = {};
b_ODtoHbOHbR_DPF1.help    = {'adjust differential pathlength factor depending on the wavelength and age of the subject.',...
    'DOI: 10.1117/1.JBO.18.10.105004'};


%Choice
C_ODtoHbOHbR_DPF         = cfg_choice;
C_ODtoHbOHbR_DPF.tag     = 'C_ODtoHbOHbR_DPF';
C_ODtoHbOHbR_DPF.name    = 'DPF method';
C_ODtoHbOHbR_DPF.values  = {b_ODtoHbOHbR_DPF1,b_ODtoHbOHbR_DPF2,b_ODtoHbOHbR_DPF3 };
C_ODtoHbOHbR_DPF.val     = {b_ODtoHbOHbR_DPF1};
C_ODtoHbOHbR_DPF.help    = {'adjust differential pathlenght factor (DPF) in Modify Beer Lambert Law calculation.'};


% Executable Branch
ODtoHbOHbR      = cfg_exbranch;       
ODtoHbOHbR.name = 'Modified Beer Lambert Law';             
ODtoHbOHbR.tag  = 'ODtoHbOHbR'; 
ODtoHbOHbR.val  = {NIRSmat DelPreviousData  PVF C_ODtoHbOHbR_DPF }; 
ODtoHbOHbR.prog = @nirs_run_ModifyBeerLambertLaw;  
ODtoHbOHbR.vout = @nirs_cfg_vout_ModifyBeerLambertLaw; 
ODtoHbOHbR.help = {'Convert Normalized intensity to HbO/HbR.',...
            'Input: dOD intensity',...
            'Output: Concentrations HbO and HbR from',...
            'Use the wavelength extenction coefficent',...
            'Citations:',...
            'W. B. Gratzer, Med. Res. Council Labs, Holly Hill,London',...
            'N. Kollias, Wellman Laboratories, Harvard Medical School, Boston',...
            '',...
            'Normalised intensity I/Io will be convert in dOD = log(I/Io)',...
            'then the Modified Beer Lamber Law will be apply directly.',...
            'Choose the DPF selection method'};          

function vout = nirs_cfg_vout_ModifyBeerLambertLaw(job)
    vout = cfg_dep;                     
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 4b Artifact detection Step Detection Normalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_norm_dur         = cfg_entry;
i_norm_dur.name    = 'Duration to calculate the offset (s)';
i_norm_dur.tag     = 'i_norm_dur';       
i_norm_dur.strtype = 'r';
i_norm_dur.num     = [1 1];
i_norm_dur.val     = {2.5}; 
i_norm_dur.help    = {'Set amplitude variation threshold between 2 good sub interval group the interval or separate it'};




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODULE 5a Processing Epoch Averaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savenirs         = cfg_menu;
savenirs.name    = 'Save in .nirs format';
savenirs.tag     = 'savenirs';      
savenirs.labels  = {'Yes','No'};
savenirs.values  = {1,0};
savenirs.val     = {1};
savenirs.help    = {'Save data in .nirs format (compatible with HomER).'};


avtype          = cfg_menu;
avtype.tag      = 'avtype';
avtype.name     = 'Average over';
avtype.labels   = {'Multiple files'};
avtype.values   = {1};
avtype.val      = {1};
avtype.help     = {'Use all the files in the NIRS.mat.'};

avg_datatype         = cfg_menu;
avg_datatype.name    = 'Data type';
avg_datatype.tag     = 'avg_datatype';      
avg_datatype.labels  = {'DC','AC','PH'};
avg_datatype.values  = {2,1,0};
avg_datatype.val     = {2};
avg_datatype.help    = {'Boxy data give DC, AC and PH component. Select the one you wish to average.'};


badintervalratio         = cfg_entry;
badintervalratio.tag     = 'badintervalratio';
badintervalratio.name    = 'Reject trial ratio';
badintervalratio.strtype = 'r';
badintervalratio.num     = [1 Inf];
badintervalratio.val     = {0.5};
badintervalratio.help    = {'Reject the trial if more than xx% of is duration is marked as a bad interval, set to 1 to keep trial.'};

badchannelratio         = cfg_entry;
badchannelratio.tag     = 'badchannelratio';
badchannelratio.name    = 'Reject channel ratio';
badchannelratio.strtype = 'r';
badchannelratio.num     = [1 inf];
badchannelratio.val     = {0.5};
badchannelratio.help    = {'Reject the channel if less than xx% of the trial are rejected, set to 0 to keep all channels.'};

helpmemoryprob          = cfg_menu;
helpmemoryprob.tag      = 'helpmemoryprob';
helpmemoryprob.name     = 'Do you need to solve help memory problem';
helpmemoryprob.labels   = {'No','No'};
helpmemoryprob.values   = {0,0};
helpmemoryprob.val      = {0};
helpmemoryprob.help     = {'Not available now, Help to solve help memory problem especially usefull for large .nir data file, each segment trial will be read independantly in spite of loading the whole file'};

e_base_pretime         = cfg_entry;
e_base_pretime.name    = 'PreTime baseline correction';
e_base_pretime.tag     = 'e_base_pretime';       
e_base_pretime.strtype = 'r';
e_base_pretime.num     = [1 1];
e_base_pretime.val     = {0}; 
e_base_pretime.help    = {'Pretime for the baseline correction (Trigger reference = time 0)'}; 

e_base_posttime         = cfg_entry;
e_base_posttime.name    = 'PostTime baseline correction';
e_base_posttime.tag     = 'e_base_posttime';       
e_base_posttime.strtype = 'r';
e_base_posttime.num     = [1 1];
e_base_posttime.val     = {0}; 
e_base_posttime.help    = {'Postime for the baseline correction (Trigger reference = time 0)'}; 

m_mean          = cfg_menu;
m_mean.tag      = 'm_mean';
m_mean.name     = 'Subtract value';
m_mean.labels   = {'Mean','Median'};
m_mean.values   = {0,1};
m_mean.val      = {0};
m_mean.help     = {'Subtract value using mean or median'};

b_manualbaseline_corr        = cfg_branch;
b_manualbaseline_corr.tag    = 'b_manualbaseline_corr';
b_manualbaseline_corr.name   = 'Manual';
b_manualbaseline_corr.val    = {m_mean,e_base_pretime,e_base_posttime};
b_manualbaseline_corr.help   = {'Manual time baseline correction for each stim'}; 

m_nobaseline_corr           = cfg_menu;
m_nobaseline_corr.tag       = 'm_nobaseline_corr';
m_nobaseline_corr.name      = 'No baseline correction';
m_nobaseline_corr.labels    = {'No'};
m_nobaseline_corr.values    = {0};
m_nobaseline_corr.val       = {0};
m_nobaseline_corr.help      = {'Keep the data original and just average without appling baseline correction'};

m_defaultbaseline_corr          = cfg_menu;
m_defaultbaseline_corr.tag      = 'm_defaultbaseline_corr';
m_defaultbaseline_corr.name     = 'Subtract PreTime';
m_defaultbaseline_corr.labels   = {'Yes'};
m_defaultbaseline_corr.values   = {1};
m_defaultbaseline_corr.val      = {1};
m_defaultbaseline_corr.help     = {'Subtract PreTime to trig to each event.'};

c_baseline_corr         = cfg_choice;
c_baseline_corr.tag     = 'c_baseline_corr';
c_baseline_corr.name    = 'Baseline Correction';
c_baseline_corr.values  = {m_defaultbaseline_corr,b_manualbaseline_corr,m_nobaseline_corr};
c_baseline_corr.val     = {m_defaultbaseline_corr}; %Default option
c_baseline_corr.help    = {''};

m_Tvalueoption        = cfg_menu;
m_Tvalueoption.tag    = 'm_Tvalueoption';
m_Tvalueoption.name   = 'Tvalue options';
m_Tvalueoption.labels = {'Against 0', 'Against mean baseline'};% 'Against worst baseline' 
m_Tvalueoption.values = {0,2};
m_Tvalueoption.val    = {2};
m_Tvalueoption.help   = {'Simple t-test against 0', 'Simple t-test against mean baseline value'}';

m_noreject_trial      = cfg_menu;
m_noreject_trial.tag  = 'm_noreject_trial';
m_noreject_trial.name = 'Keep trial';
m_noreject_trial.labels = {'Yes'};
m_noreject_trial.values = {1};
m_noreject_trial.val  = {1};
m_noreject_trial.help = {'Keep all trials, except if the trial is more noisy than the reject trial ratio'}; 


e_reject_outlier_threshold         = cfg_entry;
e_reject_outlier_threshold.name    = 'Threshold z-score';
e_reject_outlier_threshold.tag     = 'e_reject_outlier_threshold';       
e_reject_outlier_threshold.strtype = 'r';
e_reject_outlier_threshold.num     = [1 1];
e_reject_outlier_threshold.val     = {6};
e_reject_outlier_threshold.help    = {''};

m_reject_outlier_printreport        = cfg_menu;
m_reject_outlier_printreport.tag    = 'm_reject_outlier_printreport';
m_reject_outlier_printreport.name   = 'Print report';
m_reject_outlier_printreport.labels = {'No','Yes'};
m_reject_outlier_printreport.values = {0,1};
m_reject_outlier_printreport.val    = {1};
m_reject_outlier_printreport.help   = {'Print a report of rejected trial'}; 


b_reject_trial        = cfg_branch;
b_reject_trial.tag    = 'b_reject_trial';
b_reject_trial.name   = 'Reject outlier trial z-score';
b_reject_trial.val    = {e_reject_outlier_threshold,m_reject_outlier_printreport};
b_reject_trial.help   = {'Warning: Bad trials will be set as NaN, but no traces of which trials are rejected are writen in the vmrk file. Futher versions may include track of excluded channel in DisplayGui'};


c_rejecttrial         = cfg_choice;
c_rejecttrial.tag     = 'c_rejecttrial';
c_rejecttrial.name    = 'Reject outlier trial';
c_rejecttrial.values  = {m_noreject_trial,b_reject_trial};
c_rejecttrial.val     = {m_noreject_trial}; %Default option
c_rejecttrial.help    = {''};
%

choiceave         = cfg_branch;
choiceave.tag     = 'choiceave';
choiceave.name    = 'Averaging options: ';
choiceave.val     = {avtype trigger pretime posttime badintervalratio badchannelratio,avg_datatype,c_baseline_corr,c_rejecttrial,m_Tvalueoption};
choiceave.help    = {''}';

% Executable Branch
E_average      = cfg_exbranch;
E_average.name = 'Epoch averaging';
E_average.tag  = 'E_average';
E_average.val  = {NIRSmat DelPreviousData savenirs choiceave};
E_average.prog = @nirs_run_average;
E_average.vout = @nirs_cfg_vout_average;
E_average.help = {'Average over many epochs.',...
    'WARNING: For multiple subjects average, all subjects must have the exact same montage and triggers.',...
    'Input & Output: All types'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_average(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


timingfile         = cfg_files; %Select raw BOXY data files for this subject 
timingfile.name    = 'Select timing file'; % The displayed name
timingfile.tag     = 'timingfile';       %file names
timingfile.ufilter = '.m*';    %BOXY files are labeled .001, .002, .003, etc. 
%and it is very unlikely that there are more than 99 of them
timingfile.num     = [1 1];     % Number of inputs required 
timingfile.help    = {'Timing file.'}; 

% % Executable Branch
% PCA      = cfg_exbranch;
% PCA.name = 'PCA cleaning';
% PCA.tag  = 'PCAcleaning';
% PCA.val  = {NIRSmat DelPreviousData  timingfile};
% PCA.prog = @nirs_run_PCA;
% PCA.vout = @nirs_cfg_vout_PCA;
% PCA.help = {'À ENLEVER Apply PCA on selected time'};
% 
% function vout = nirs_cfg_vout_PCA(job)
%     vout = cfg_dep;                    
%     vout.sname      = 'NIRS.mat';       
%     vout.src_output = substruct('.','NIRSmat'); 
%     vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 6a WRITE NIRS
%%%%%%%%%%%%%%%%%%%%%%%%%

NIRSsession         = cfg_entry;
NIRSsession.name    = 'Session'; 
NIRSsession.tag     = 'NIRSsession';       
NIRSsession.strtype = 's';
NIRSsession.num     = [1 inf];
NIRSsession.val     = {'end'};
NIRSsession.help    = {'Enter the session number or the module in the processing NIRS file.'};

NIRSname         = cfg_entry;
NIRSname.name    = 'File output';
NIRSname.tag     = 'FileOutput';       
NIRSname.strtype = 's';
NIRSname.num     = [1 inf];
NIRSname.val     = {'name'};
NIRSname.help    = {'OutputName only use to name concatenate file, else initial name will be kept.'};


NIRS_export_Concatenatefile         = cfg_branch;
NIRS_export_Concatenatefile.tag     = 'NIRS_export_Concatenatefile';
NIRS_export_Concatenatefile.name    = 'Concatenate file';
NIRS_export_Concatenatefile.val     = {NIRSname};
NIRS_export_Concatenatefile.help    = {'Export .nirs file with all segments one after the other.'}';

NIRS_export_Separatefile         = cfg_branch;
NIRS_export_Separatefile.tag     = 'NIRS_export_Separatefile';
NIRS_export_Separatefile.name    = 'Separate file';
NIRS_export_Separatefile.val     = {};
NIRS_export_Separatefile.help    = {'Export each file in a separe .nirs file. The output file is located in the same folder of the NIRS.mat file.'}';

%to be delete
NIRS_exportoption         = cfg_menu;
NIRS_exportoption.name    = 'Option';
NIRS_exportoption.tag     = 'NIRS_exportoption';      
NIRS_exportoption.labels  = {'Separate file','Concatenate file' };
NIRS_exportoption.values  = {2,1};
NIRS_exportoption.val     = {2};
NIRS_exportoption.help    = {'Export each file in a separe .nirs file or in one nirs file with all segment one after the other. '};


c_NIRS_exportoption         = cfg_choice;
c_NIRS_exportoption.tag     = 'c_NIRS_exportoption';
c_NIRS_exportoption.name    = 'Options';
c_NIRS_exportoption.values  = {NIRS_export_Separatefile,NIRS_export_Concatenatefile };
c_NIRS_exportoption.val     = {NIRS_export_Separatefile}; %Default option
c_NIRS_exportoption.help    = {''};


% Executable Branch
E_writeNIRS      = cfg_exbranch;
E_writeNIRS.name = 'Write NIRS';
E_writeNIRS.tag  = 'WriteNIRS';
E_writeNIRS.val  = {NIRSsession,NIRSmat,c_NIRS_exportoption};
E_writeNIRS.prog = @nirs_run_writenirs;
E_writeNIRS.vout = @nirs_cfg_vout_writenirs;
E_writeNIRS.help = {'Write in .nirs format for Homer.'};

function vout = nirs_cfg_vout_writenirs(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end




%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 6b WRITE HMR
%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
E_writeHMR      = cfg_exbranch;
E_writeHMR.name = 'Write HMR';
E_writeHMR.tag  = 'E_writeHMR';
E_writeHMR.val  = {NIRSmat,NIRSname,prjfile};
E_writeHMR.prog = @nirs_run_writehmr;
E_writeHMR.vout = @nirs_cfg_vout_writehmr;
E_writeHMR.help = {'Write in NIRS.mat epoch average session in session for Homer .hmr.',...
    'This function will use the last epoch averaging file in the .nirs file.',...
    'It will cause and error if there are any epoch averaging in the NIRS.mat session.' };

function vout = nirs_cfg_vout_writehmr(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 6c WRITE NIR segment
%%%%%%%%%%%%%%%%%%%%%%%%%
NIR_FileIn          = cfg_files;
NIR_FileIn.name    = 'Enter .NIR to segment'; 
NIR_FileIn.tag     = 'NIR_FileIn';       %file names
NIR_FileIn.filter  = 'nir';
NIR_FileIn.ufilter = '.nir$';    
NIR_FileIn.num     = [1 Inf];     % Number of inputs required 
NIR_FileIn.help    = {'Select a .NIR file to extract a smaller time segment and create a new .NIR file.'}; 


NIR_START_TIME         = cfg_entry;
NIR_START_TIME.name    = 'Time start';
NIR_START_TIME.tag     = 'NIR_START_TIME';       
NIR_START_TIME.strtype = 'r';
NIR_START_TIME.num     = [1 1];
NIR_START_TIME.val     = {0};
NIR_START_TIME.help    = {'Time start selection.'};

NIR_STOP_TIME         = cfg_entry;
NIR_STOP_TIME.name    = 'Time stop';
NIR_STOP_TIME.tag     = 'NIR_STOP_TIME';       
NIR_STOP_TIME.strtype = 'r';
NIR_STOP_TIME.num     = [1 1];
NIR_STOP_TIME.val     = {0};
NIR_STOP_TIME.help    = {'Time end of the selection.'};


% Executable Branch
E_NIR_segment      = cfg_exbranch;
E_NIR_segment.name = 'Write NIR indivual file segment';
E_NIR_segment.tag  = 'Write_NIR_segment';
E_NIR_segment.val  = {NIRSmat, NIR_FileIn, NIR_START_TIME,NIR_STOP_TIME};
E_NIR_segment.prog = @nirs_run_E_NIR_segment;
E_NIR_segment.vout = @nirs_run_vout_E_NIR_segment;
E_NIR_segment.help = {'Use a .NIR file to write a smaller segment, can be use to specify only the baseline for the PCA artifact for example.'};

function vout = nirs_run_vout_E_NIR_segment(job)
    vout            = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 6d WRITE .snirf file
%%%%%%%%%%%%%%%%%%%%%%%%%

SNIRFname         = cfg_entry;
SNIRFname.name    = 'Output filename';
SNIRFname.tag     = 'FileOutput';       
SNIRFname.strtype = 's';
SNIRFname.num     = [1 inf];
SNIRFname.help    = {'Used to name the resulting output SNIRF file.'};

SNIRF_outpath         = cfg_entry; %path
SNIRF_outpath.name    = 'Folder path for output SNIRF file';
SNIRF_outpath.tag     = 'SNIRF_outpath';       
SNIRF_outpath.strtype = 's';
SNIRF_outpath.num     = [1 Inf];     
SNIRF_outpath.val     = {'C:\data\'}; 
SNIRF_outpath.help    = {'The path of the folder in which the SNIRF file is going to be created.',...
    'Ex: C:\data\user1\SNIRFfolder\'}; 

f_SNIRFfile         = cfg_files;
f_SNIRFfile.name    = 'Select existing SNIRF file'; 
f_SNIRFfile.tag     = 'f_SNIRFfile';       %file names
f_SNIRFfile.filter  = 'snirf';
f_SNIRFfile.ufilter = '.snirf$';    
f_SNIRFfile.num     = [1 1];     % Number of inputs required 
f_SNIRFfile.help    = {'Select the SNIRF file to which you wish to append your data'}; 

B_SNIRFAppend       = cfg_branch;
B_SNIRFAppend.name  = 'Append to an existing SNIRF file';
B_SNIRFAppend.tag   = 'B_SNIRFAppend';
B_SNIRFAppend.val   = {f_SNIRFfile};

B_SNIRFCreate       = cfg_branch;
B_SNIRFCreate.name  = 'Create a new SNIRF file';
B_SNIRFCreate.tag   = 'B_SNIRFCreate';
B_SNIRFCreate.val   = {SNIRFname SNIRF_outpath};
B_SNIRFCreate.help  = {'Export the data to a new SNIRF file.'};

B_SNIRFExportLocation         = cfg_choice;
B_SNIRFExportLocation.tag     = 'B_SNIRFExportLocation';
B_SNIRFExportLocation.name    = 'Choose SNIRF export type';
B_SNIRFExportLocation.val     = {B_SNIRFCreate};
B_SNIRFExportLocation.values  = {B_SNIRFCreate  B_SNIRFAppend};
B_SNIRFExportLocation.help    = {'Choose whether you want to export your data to a new SNIRF file or append it to an existing SNIRF file.'};

% Executable Branch
E_writeSNIRF      = cfg_exbranch;
E_writeSNIRF.name = 'Write SNIRF';
E_writeSNIRF.tag  = 'E_writeSNIRF';
E_writeSNIRF.val  = {NIRSmat,B_SNIRFExportLocation};
E_writeSNIRF.prog = @nirs_run_writeSNIRF;
E_writeSNIRF.vout = @nirs_cfg_vout_writeSNIRF;
E_writeSNIRF.help = {'Write session in Shared Near Infrared File Format Specification to export the data.',...
    'This will produce data format .snirf' };

function vout = nirs_cfg_vout_writeSNIRF(job)
    vout = cfg_dep;                    
    vout.sname      = 'SNIRF.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

% f_Pipeline         = cfg_files;
% f_Pipeline.name    = 'Select output folder to Batch pipeline.m'; 
% f_Pipeline.tag     = 'f_Pipeline';       %file names
% f_Pipeline.filter  = {'dir'};
% f_Pipeline.ufilter = '.*';    
% f_Pipeline.num     = [1 1];      % Number of inputs required 
% f_Pipeline.help    = {'Select the output location of your Batch .m script pipeline'}; 




% % Executable Branch
% E_createPipeline       = cfg_exbranch;
% E_createPipeline.name  = 'Create script pipeline';
% E_createPipeline.tag   = 'f_Pipeline';
% E_createPipeline.val   = {f_Pipeline};
% E_createPipeline.prog  = @nirs_run_createPipeline;
% %E_createPipeline.vout = @nirs_cfg_vout_createPipeline;
% E_writeSNIRF.help      = {'Write session in Shared Near Infrared File Format Specification to export the data',...
%     'This will produce data format .snirf' };

% function vout = nirs_cfg_vout_createPipeline(job)
%     vout = cfg_dep;                    
%     vout.sname      = 'SNIRF.mat';       
%     vout.src_output = substruct('.','NIRSmat'); 
%     vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 7 Display modules GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_GUI      = cfg_exbranch;
E_GUI.name = 'DisplayGUI';
E_GUI.tag  = 'E_GUI';
E_GUI.val  = {NIRSmat, DelPreviousData};
E_GUI.prog = @nirs_run_GUI;
E_GUI.vout = @nirs_cfg_vout_GUI;
E_GUI.help = {'This function opens the NIRS.mat structure for visualization from any step of the analysis, as well as optional multimodal data. It opens this interface to display the data with many visual features. Only the last open module could be edited for manual artifact adjustments.'};

function vout = nirs_cfg_vout_GUI(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 7b Display modules Video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HMRdata         = cfg_files; %Select NIRS.mat for this subject 
HMRdata.name    = 'FILE.hmr'; % The displayed name
HMRdata.tag     = 'HMRdata';       %file names
HMRdata.filter  = 'hmr';
HMRdata.ufilter = '.hmr$';    
HMRdata.num     = [1 Inf];     % Number of inputs required 
HMRdata.help    = {'Select NIRS.mat for the subject(s).'}; 

HMR_cffile          = cfg_entry;
HMR_cffile.tag      = 'HMR_cffile';
HMR_cffile.name     = 'File condition';
HMR_cffile.strtype  = 's';
HMR_cffile.num      = [1 Inf];     
HMR_cffile.val      = {'CAT'}; 
HMR_cffile.help     = {'File condition .hmr input'};

HMR_datatype           = cfg_menu;
HMR_datatype.tag       = 'HMR_datatype';
HMR_datatype.name      = 'Output type';     
HMR_datatype.labels    = {'dconc','dconc avgtime','dconc z-score','dconc z-score avgtime','dconc minmax','dconc minmax avgtime'}; 
HMR_datatype.values    = {1,2,3,4,5,6};
HMR_datatype.val       = {1};
HMR_datatype.help      = {'dconc, ... '};


HMR_datad1mat             = cfg_files; %Select NIRS.mat for this subject 
HMR_datad1mat.name        = 'D1mat'; % The displayed name
HMR_datad1mat.tag         = 'HMR_datad1mat';       %file names
HMR_datad1mat.filter      = 'mat';
HMR_datad1mat.ufilter     = '.mat$';    
HMR_datad1mat.num         = [0 Inf];     % Number of inputs required 
HMR_datad1mat.help        = {'Select a matrix of data to display in the same order than the subject in the hmr file. Let empty if you want to select concentration data type'}; 

HMR_savenormfig          = cfg_menu;
HMR_savenormfig.tag      = 'HMR_savenormfig';
HMR_savenormfig.name     = 'Save Normalised Figure';  
HMR_savenormfig.labels   = {'No','Yes'}; 
HMR_savenormfig.values   = {0,1};
HMR_savenormfig.val      = {1};
HMR_savenormfig.help     = {'Save concentration and normalised concentration figure'};

b_HMRdata          = cfg_branch;
b_HMRdata.tag      = 'b_HMRdata';
b_HMRdata.name     = 'HMR data';
%b_HMRdata.val     = {HMRdata HMR_cffile HMR_datatype HMR_datad1mat HMR_savenormfig};
b_HMRdata.val      = {HMRdata HMR_cffile HMR_datatype HMR_savenormfig};
b_HMRdata.help     = {''}';


b_HMRdatad1manual         = cfg_branch;
b_HMRdatad1manual.tag     = 'b_HMRdatad1manual';
b_HMRdatad1manual.name    = 'HMR ..mat data';
b_HMRdatad1manual.val     = {HMRdata HMR_cffile HMR_datatype HMR_datad1mat HMR_savenormfig};
b_HMRdatad1manual.help    = {''}';


NIRS_datatype           = cfg_menu;
NIRS_datatype.tag       = 'NIRS_datatype';
NIRS_datatype.name      = 'Output type';
NIRS_datatype.labels    = {'Epoch mean','Epoch mean avgtime','Epoch Tval','Epoch Tval avgtime'};
NIRS_datatype.values    = {1,2,3,4};
NIRS_datatype.val       = {1};
NIRS_datatype.help      = {'Epoch averaging mean, Epoch averaging mean over time the period start : start+steptime, Epoch averaging mean Tval '};

b_NIRSdata        = cfg_branch;
b_NIRSdata.tag    = 'b_NIRSdata';
b_NIRSdata.name   = 'NIRS epoch average data';
b_NIRSdata.val    = {NIRSmat NIRS_datatype HMR_savenormfig};
b_NIRSdata.help   = {''};

b_SPMdata        = cfg_branch;
b_SPMdata.tag    = 'b_SPMdata';
b_SPMdata.name   = 'NIRS GLM data';
b_SPMdata.val    = {NIRSmat};
b_SPMdata.help   = {''}';

vColordata         = cfg_files; %Select NIRS.mat for this subject 
vColordata.name    = 'Open topographic map'; % The displayed name
vColordata.tag     = 'vColordata';       
vColordata.filter  = '.*';
vColordata.ufilter = '.img|TOPO.mat';    
vColordata.num     = [1 Inf];     % Number of inputs required 
vColordata.help    = {'Select a serie of topographic file as .img or structure TOPO.mat. ',...
    'If you have selected .img file note that they should have be done on the same skin or cortex reconstruction.',...
    'Only one .prj can be use then. '}; 

b_vColordata        = cfg_branch;
b_vColordata.tag    = 'b_vColordata';
b_vColordata.name   = 'Load topographic map';
b_vColordata.val    = {vColordata};
b_vcolordata.help   = {''};

typedata         = cfg_choice;
typedata.tag     = 'typedata';
typedata.name    = 'Design';
typedata.val     = {b_NIRSdata};
typedata.help    = {''};
typedata.values  = {b_NIRSdata,b_HMRdata,b_HMRdatad1manual, b_SPMdata,b_vColordata};

%prjfile define above 
m_prjfile_mode          = cfg_menu;
m_prjfile_mode.tag      = 'm_prjfile_mode';
m_prjfile_mode.name     = 'Option prj';  
m_prjfile_mode.labels   = {'Paired', 'Same prj for all file'}; 
m_prjfile_mode.values   = {0,1};
m_prjfile_mode.val      = {0};
m_prjfile_mode.help     = {'The project file (.prj) use to created the topography can be paired which all file selected above, file1 use prj1, file2 use prj2...',...
    'You can also use the same prj for all the topography',...
    'In case of the selection Load topographic map is selected, the same prj skin and cortex have be use for all file'};

v_cmin          =  cfg_entry;
v_cmin.tag      = 'v_cmin';
v_cmin.name     = 'Min intensity';
v_cmin.val      = {-1};
v_cmin.strtype  = 'r';  
v_cmin.num      = [1 1]; 
v_cmin.help     = {'Min intensity'};

v_cmax          = cfg_entry;
v_cmax.tag      = 'v_cmax';
v_cmax.name     = 'Max intensity';
v_cmax.val      = {1};
v_cmax.strtype  = 'r';  
v_cmax.num      = [1 1]; 
v_cmax.help     = {'Max intensity'};

v_ctr          = cfg_entry;
v_ctr.tag      = 'v_ctr';
v_ctr.name     = 'Threshold intensity';
v_ctr.val      = {0.1};
v_ctr.strtype  = 'r';  
v_ctr.num      = [1 1]; 
v_ctr.help     = {'Threshold intensity'};

v_view          = cfg_entry;
v_view.tag      = 'v_view';
v_view.name     = 'View angle';
v_view.num      = [0 inf];
v_view.val      =  {'90,30;-90,30;180,30'};
v_view.strtype  = 's';  
v_view.help     = {['Write view angle list as : ''head rotation angle,inclinaison angle;...''  as example ''90,30;-90,30'' ',...
    'will produce a left view with 30 degree of inclinaison and a right view with 30 degree of inclinaison. ',...
    'Rotation angle : 0 front, 90 left, 180, back, -90 right view ',...
    'Inclinaison angle : 0 no inclinaison, 90 top view ',....
    'If you enter notting no figure or video will be produce. ']};

v_start          =  cfg_entry;
v_start.tag      = 'v_start';
v_start.name     = 'Start timing (s)';
v_start.val      = {0};
v_start.strtype  = 'r';  
v_start.num      = [1 1]; 
v_start.help     = {'Start timing (s)'};

v_stop          =  cfg_entry;
v_stop.tag      = 'v_stop';
v_stop.name     = 'Stop timing (s)';
v_stop.val      = {1};
v_stop.strtype  = 'r';  
v_stop.num      = [1 1]; 
v_stop.help     = {'Stop timing (s)'};

v_step          =  cfg_entry;
v_step.tag      = 'v_step';
v_step.name     = 'Step timing (s)';
v_step.val      = {1};
v_step.strtype  = 'r';  
v_step.num      = [1 1]; 
v_step.help     = {'Step timing (s)'};

v_outpath         = cfg_entry; %path
v_outpath.name    = 'path for .video output files';
v_outpath.tag     = 'v_outpath';       
v_outpath.strtype = 's';
v_outpath.num     = [1 Inf];     
v_outpath.val     = {'video'}; 
v_outpath.help    = {['Path for video output files: write video_cond1(omit backslashes) data will be save in subject folder like  ...subject\dataspm\video\ ',...
    'or write complete directory name as c:/data/video_cond1 to define manual the path']}; 

v_outtimeprecisionleft         = cfg_entry;
v_outtimeprecisionleft.name    = 'Number of left decimal output name time';
v_outtimeprecisionleft.tag     = 'v_outtimeprecisionleft';       
v_outtimeprecisionleft.strtype = 'r';
v_outtimeprecisionleft.num     = [1 Inf];     
v_outtimeprecisionleft.val     = {4}; 
v_outtimeprecisionleft.help    = {['Output file name will be write with this number of decimal exemple : 4=>0012.XX(s), 3 => 012.XXX(s)']}; 

v_outtimeprecision         = cfg_entry;
v_outtimeprecision.name    = 'Number of decimal output name time';
v_outtimeprecision.tag     = 'v_outtimeprecision';       
v_outtimeprecision.strtype = 'r';
v_outtimeprecision.num     = [1 Inf];     
v_outtimeprecision.val     = {2}; 
v_outtimeprecision.help    = {['Output file name will be write with this number of decimal exemple : 2=>0XX.23(s), 0 => 0XX(s)']}; 

v_skin        = cfg_menu;
v_skin.tag    = 'v_skin';
v_skin.name   = 'Select projection on skin or on cortex';
v_skin.labels = {'Skin','Cortex'};
v_skin.values = {0,1};
v_skin.val    = {0};
v_skin.help   = {'Display topography of the activation on the skin or the cortex'}';

v_showcover        = cfg_menu;
v_showcover.tag    = 'v_showcover';
v_showcover.name   = 'Show coverage';
v_showcover.labels = {'Yes','No'};
v_showcover.values = {1,0};
v_showcover.val    = {1};
v_showcover.help   = {'Show uncovered region by the montage in grey'}';

b_option         = cfg_branch;
b_option.tag     = 'b_option';
b_option.name    = 'Define display option';
b_option.help    = {'Display option for HSJ helmet.'};
b_option.val     = {prjfile m_prjfile_mode v_cmin,v_cmax,v_ctr,v_view,v_start,v_stop,v_step,v_outpath,v_outtimeprecisionleft,v_outtimeprecision,v_showcover,v_skin};

% Executable Branch
E_VIDEO      = cfg_exbranch;
E_VIDEO.name = 'Display Video';
E_VIDEO.tag  = 'E_VIDEO';
E_VIDEO.val  = {typedata,b_option};
E_VIDEO.prog = @nirs_run_E_VIDEO;
E_VIDEO.vout = @nirs_cfg_vout_E_VIDEO;
E_VIDEO.help = {'Display GUI.'};

function vout = nirs_cfg_vout_E_VIDEO(job)
    vout            = cfg_dep;                    
    vout.sname      = 'TOPO.mat';       
    vout.src_output = substruct('.','TOPOmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 7c  Display NIRS.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_viewNIRS      = cfg_exbranch;
E_viewNIRS.name = 'Display NIRS.mat';
E_viewNIRS.tag  = 'E_viewNIRS';
E_viewNIRS.val  = {NIRSmat};
E_viewNIRS.prog = @nirs_run_E_viewNIRS;
E_viewNIRS.vout = @nirs_cfg_vout_E_viewNIRS;
E_viewNIRS.help = {'Display NIRS.mat GUI.'};

function vout = nirs_cfg_vout_E_viewNIRS(~)
    vout            = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 9a Artifact extract component from noise event (bad interval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_extractcomponent         = cfg_entry;
m_extractcomponent.tag     = 'm_extractcomponent';
m_extractcomponent.name    = 'Extract';
m_extractcomponent.strtype = 's';
m_extractcomponent.num     = [1 Inf];
m_extractcomponent.val     = {'MVTPCA'}; 
m_extractcomponent.help    = {'MVTPCA serves as a label of identification to recognize extracted PCA components in the subtract or in the export function.'};





m_extractcomponentfigure        = cfg_menu;
m_extractcomponentfigure.tag    = 'm_extractcomponentfigure';
m_extractcomponentfigure.name   = 'Figure ';
m_extractcomponentfigure.labels = {'Yes', 'No'};
m_extractcomponentfigure.values = {1,0};
m_extractcomponentfigure.val    = {1}; 
m_extractcomponentfigure.help   = {'Choose whether or not you want to display the output figure.'};


i_extract_pourcentagech         = cfg_entry;
i_extract_pourcentagech.name    = 'Percentage minimal noisy channel to be an new component';
i_extract_pourcentagech.tag     = 'i_extract_pourcentagech';       
i_extract_pourcentagech.strtype = 'r';
i_extract_pourcentagech.num     = [1 1];
i_extract_pourcentagech.val     = {5}; 
i_extract_pourcentagech.help    = {'Consider this event to extract components if they have at least 5% of noisy channel for this time window.'};


i_extractnoiseupto_nbPCA         = cfg_entry;
i_extractnoiseupto_nbPCA.name    = 'Components to remove up to XX% of the variance in the data';
i_extractnoiseupto_nbPCA.tag     = 'i_extractnoiseupto_nbPCA';       
i_extractnoiseupto_nbPCA.strtype = 'r';
i_extractnoiseupto_nbPCA.num     = [1 Inf];
i_extractnoiseupto_nbPCA.val     = {97}; 
i_extractnoiseupto_nbPCA.help    = {'Find number of component to identify as artifact, one or more, in other to explain up to 97% of the variance in the data'};


i_extractnoise_nbPCA         = cfg_entry;
i_extractnoise_nbPCA.name    = 'Minimal percentage of variance to explain';
i_extractnoise_nbPCA.tag     = 'i_extractnoise_nbPCA';       
i_extractnoise_nbPCA.strtype = 'r';
i_extractnoise_nbPCA.num     = [1 Inf];
i_extractnoise_nbPCA.val     = {1}; 
i_extractnoise_nbPCA.help    = {' Find number of components to identify as an artifact in PCA setting the minimal xx % of explained variance.'};

b_extractcomponent_PCA        = cfg_branch;
b_extractcomponent_PCA.tag    = 'b_extractcomponent_PCA';
b_extractcomponent_PCA.name   = 'Identify PCA for artifact period';
b_extractcomponent_PCA.val    = {NIRSmat,  m_extractcomponent,m_extractcomponentfigure,i_extract_pourcentagech, i_extractnoiseupto_nbPCA, i_extractnoise_nbPCA};
b_extractcomponent_PCA.help   = {'First, identify noisy intervals using artifact detection or a manual revision. This function runs a PCA decomposition (targetPCA) on each bad interval (yellow segment in the DisplayGUI). The decomposition is performed on identified channels during a continuous bad interval. PCA decomposition sort components according to the explained variance. The component(s) explaining the highest variance during the artifactual interval is assumed to be mainly related to the artifact event. They will be stored in the component list with the label MVTPCA. We recommend using the module ‘Subtract Components’ that will subtract all the components identified with a specific label.'};

f_extractcomponent_physzone         = cfg_files;
f_extractcomponent_physzone.name    = 'Enter Regressor zone'; 
f_extractcomponent_physzone.tag     = 'f_extractcomponent_physzone';       %file names
f_extractcomponent_physzone.filter  = 'zone';
f_extractcomponent_physzone.ufilter = '.zone';    
f_extractcomponent_physzone.num     = [1 Inf];     % Number of inputs required regression 
f_extractcomponent_physzone.help    = {'Use to define channel use in the short distance regression: regressor zone1 will regress in zone1 '}; 

m_extractcomponent_physzone        = cfg_menu;
m_extractcomponent_physzone.tag    = 'm_extractcomponent_physzone';
m_extractcomponent_physzone.name   = 'Use ';
m_extractcomponent_physzone.labels = {'Mean' ,'PCA'};
m_extractcomponent_physzone.values = {0,1};
m_extractcomponent_physzone.val    = {0}; 
m_extractcomponent_physzone.help   = {'Use the mean or the first principal component of the regressor zone; usefull when the regressor zone contain many channels'};

b_extractcomponent_phys      = cfg_branch;
b_extractcomponent_phys.tag  = 'b_extractcomponent_phys';
b_extractcomponent_phys.name = 'Identify physiology regression';
b_extractcomponent_phys.val  = {NIRSmat,  f_extractcomponent_physzone ,m_extractcomponent_physzone};
b_extractcomponent_phys.help = {'Apply substraction of the short distance physiology as describe in Saager and Berger 2008, https://doi.org/10.1117/1.2940587',...
                                    'One signal was scaled to fit the other in a least-squares LS sense the scale estimation is save as (SHORTGLM) component and must be subtract'};

f_extractcomponent_glmlist         = cfg_files;
f_extractcomponent_glmlist.name    = 'List GLM to identify (xls)'; 
f_extractcomponent_glmlist.tag     = 'f_extractcomponent_glmlist';       %file names
f_extractcomponent_glmlist.filter  = {'xlsx','xls','txt'};
f_extractcomponent_glmlist.ufilter = '.*';
f_extractcomponent_glmlist.num     = [1 Inf];     % Number of inputs required 
f_extractcomponent_glmlist.help    = {'Enter the list xls to extract GLM, the list must include the following columns:',...
    '''NIRS.mat folder'': directory of the NIRS.mat to use;',...
    '''file'': number to identify the file to use;',...
    '''tStart'': time start in seconds;',...
    '''tStop'': time stop in seconds, the multiple regression will be applied on this time window;',...
    '''Label'': to include in the name of the component;',...
    'Regressor(s) (''X0'', ''X1'',''X2'',''X3'',''X4''...): the regressor could be a name in the aux list, a channel zone (regressor zone 1 apply to zone 1). '}; 

b_extractcomponent_glm          = cfg_branch;
b_extractcomponent_glm.tag      = 'b_extractcomponent_glm';
b_extractcomponent_glm.name     = 'Identify GLM';
b_extractcomponent_glm.val      = {f_extractcomponent_glmlist};
b_extractcomponent_glm.help     = {'Apply multiple linear regressions (regress.m)'};


f_extractcomponent_PARAFAClist         = cfg_files;
f_extractcomponent_PARAFAClist.name    = 'List Parafac to identify (xls)'; 
f_extractcomponent_PARAFAClist.tag     = 'f_component_PARAFAClist';       %file names
f_extractcomponent_PARAFAClist.filter  = {'xlsx','xls','txt'};
f_extractcomponent_PARAFAClist.ufilter = '.*';    
f_extractcomponent_PARAFAClist.num     = [1 Inf];     % Number of inputs required 
f_extractcomponent_PARAFAClist.help    = {'Enter the list xls to extract PARAFAC, the list must include the following columns:',...
    '''NIRS.mat folder'': directory of the NIRS.mat to use;',...
    '''file'': number to identify the file to use;',...
    '''tStart'' : time start in second;',...
    '''tStop'': time stop in second, the multiple regression will be applied on this time window;',...
    '''Label'': Label identification to write in the event (useful manage the export)'};

b_extractcomponent_PARAFAC          = cfg_branch;
b_extractcomponent_PARAFAC.tag      = 'b_extractcomponent_PARAFAC';
b_extractcomponent_PARAFAC.name     = 'Identify PARAFAC';
b_extractcomponent_PARAFAC.val      = {f_extractcomponent_PARAFAClist};
b_extractcomponent_PARAFAC.help     = {'Display option for HSJ helmet.'};


f_extractcomponent_AVGlist         = cfg_files;
f_extractcomponent_AVGlist.name    = 'List AVG to identify (xls)'; 
f_extractcomponent_AVGlist.tag     = 'f_component_AVGlist';       %file names
f_extractcomponent_AVGlist.filter  = {'xlsx','xls','txt'};
f_extractcomponent_AVGlist.ufilter = '.*';    
f_extractcomponent_AVGlist.num     = [1 Inf];     % Number of inputs required 
f_extractcomponent_AVGlist.help    = {'NIRS.mat folder: Directory to locate the data to extract.',... 
    '''File'':  1  : block in the NIRS.mat data file;',... 
    '''tStart'': to get the curve start point;',... 
    '''tStop: to get the curve stop point;',... 
    '''tStartavg'': to get the average starting point;',... 
    '''tStopavg'': to get the average stopping point;',... 
    '''Label'': Write label in the component name;',... 
    '''ZoneDisplay'': use the first zone channel to plot the average. Keep the zone file in the same folder as the excel ExtractAVG setting.'}; 

b_extractcomponent_AVG          = cfg_branch;
b_extractcomponent_AVG.tag      = 'b_extractcomponent_AVG';
b_extractcomponent_AVG.name     = 'Identify AVG';
b_extractcomponent_AVG.val      = {f_extractcomponent_AVGlist};
b_extractcomponent_AVG.help     = {'Find the average of a time period for each channel.'};

i_extractnoise_labelPARAFAC         = cfg_entry;
i_extractnoise_labelPARAFAC.name    = 'Extract';
i_extractnoise_labelPARAFAC.tag     = 'i_extractnoise_labelPARAFAC';       
i_extractnoise_labelPARAFAC.strtype = 's';
i_extractnoise_labelPARAFAC.num     = [1 Inf];
i_extractnoise_labelPARAFAC.val     = {'MVTPARAFAC'}; 
i_extractnoise_labelPARAFAC.help    = {'MVTPARAFAC serves as a label of identification to recognize extracted PARAFAC components in the subtract or in the export function.'};

i_extractnoise_nbPARAFAC         = cfg_entry;
i_extractnoise_nbPARAFAC.name    = 'Nb component to try';
i_extractnoise_nbPARAFAC.tag     = 'i_extractnoise_nbPARAFAC';       
i_extractnoise_nbPARAFAC.strtype = 'r';
i_extractnoise_nbPARAFAC.num     = [1 Inf];
i_extractnoise_nbPARAFAC.val     = {5}; 
i_extractnoise_nbPARAFAC.help    = {'Find the optimal number of component to decompose PARAFAC try up to x component to find best Concordia and minimise the error .'};

b_extractnoise_PARAFAC          = cfg_branch;
b_extractnoise_PARAFAC.tag      = 'b_extractnoise_PARAFAC';
b_extractnoise_PARAFAC.name     = 'Identify PARAFAC for artifact period';
b_extractnoise_PARAFAC.val      = {NIRSmat,  i_extractnoise_labelPARAFAC,m_extractcomponentfigure,i_extract_pourcentagech,i_extractnoise_nbPARAFAC};
b_extractnoise_PARAFAC.help     = {'Find the average of a time period for each channel.'};

c_extractcomponent         = cfg_choice;
c_extractcomponent.tag     = 'c_extractcomponent';
c_extractcomponent.name    = 'Option';
c_extractcomponent.val     = {b_extractcomponent_PCA};
c_extractcomponent.help    = {''};
c_extractcomponent.values  = {b_extractcomponent_PCA,b_extractnoise_PARAFAC, b_extractcomponent_phys, b_extractcomponent_glm,b_extractcomponent_PARAFAC,b_extractcomponent_AVG };


% Executable Branch
E_extractcomponent      = cfg_exbranch;
E_extractcomponent.name = 'Extract';
E_extractcomponent.tag  = 'E_extractcomponent';
E_extractcomponent.val  = {c_extractcomponent};
E_extractcomponent.prog = @nirs_run_E_extractcomponent;
E_extractcomponent.vout = @nirs_cfg_vout_E_extractcomponent;
E_extractcomponent.help = {'Identify data components using several data decomposition methods on target data. The component will be added to a list that could be visualized in the DisplayGUI (extract), subtract in case of an artifactual component using ‘Subtract component’ or export as a relevant activity for further statistics using ‘Export component’. '};
%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_E_extractcomponent(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 9 Artifact substract component decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_substract_and_offsetcorrection        = cfg_menu;
m_substract_and_offsetcorrection.tag    = 'm_substract_and_offsetcorrection';
m_substract_and_offsetcorrection.name   = 'Offset correction';
m_substract_and_offsetcorrection.labels = {'Yes', 'No'};
m_substract_and_offsetcorrection.values = {1,2};
m_substract_and_offsetcorrection.val    = {1}; 
m_substract_and_offsetcorrection.help   = {'Force offset adjustement after substracted artefact.'};

i_substractcomponent_label         = cfg_entry;
i_substractcomponent_label.name    = 'Component identification';
i_substractcomponent_label.tag     = 'i_substractcomponent_label';       
i_substractcomponent_label.strtype = 's';
i_substractcomponent_label.num     = [1 inf];
i_substractcomponent_label.val     = {'MVTPARAFAC'}; 
i_substractcomponent_label.help    = {'Identify the label of the component you want to substract,as example MVTPARAFAC MVTPCA for the automatic detectection'};

% Executable Branch
E_substractcomponent      = cfg_exbranch;
E_substractcomponent.name = 'Subtract';
E_substractcomponent.tag  = 'E_substractcomponent';
E_substractcomponent.val  = {NIRSmat,i_substractcomponent_label,m_substract_and_offsetcorrection};
E_substractcomponent.prog = @nirs_run_E_substractcomponent;
E_substractcomponent.vout = @nirs_cfg_vout_E_substractcomponent;
E_substractcomponent.help = {'Substract from the data the component identify by the label.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_E_substractcomponent(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODULE 9 Export component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_exportcomponent_label         = cfg_entry;
i_exportcomponent_label.name    = 'Identification';
i_exportcomponent_label.tag     = 'i_exportcomponent_label';       
i_exportcomponent_label.strtype = 's';
i_exportcomponent_label.num     = [1 inf];
i_exportcomponent_label.val     = {'M0'}; 
i_exportcomponent_label.help    = {'Identify the label to find in the component to export'};

m_exportcomponent_type        = cfg_menu;
m_exportcomponent_type.name   = 'Type';
m_exportcomponent_type.tag    = 'm_exportcomponent_type';       
m_exportcomponent_type.labels = {'GLM','PARAFAC','PCA'};
m_exportcomponent_type.values = {0,1,2};
m_exportcomponent_type.val    = {0};
m_exportcomponent_type.help   = {'Identify the label to find in the component to export'};

% Executable Branch
E_exportcomponent      = cfg_exbranch;
E_exportcomponent.name = 'Export component';
E_exportcomponent.tag  = 'E_exportcomponent';
E_exportcomponent.val  = {NIRSmat, NewDirCopyNIRSTRUE,i_exportcomponent_label,m_exportcomponent_type};
E_exportcomponent.prog = @nirs_run_E_exportcomponent;
E_exportcomponent.vout = @nirs_cfg_vout_E_exportcomponent;
E_exportcomponent.help = {'Substract from the data the component identify by the label'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_E_exportcomponent(job)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
v_component_list_remove         = cfg_entry;
v_component_list_remove.name    = 'Remove index';
v_component_list_remove.tag     = 'v_component_list_remove';       
v_component_list_remove.strtype = 'r';
v_component_list_remove.num     = [1 Inf];     
v_component_list_remove.val     = {0}; 
v_component_list_remove.help    = {'Allow removing few components that could be outliers for mean and avg topography, set zeros to keep all, else write down the index to remove. '};



f_component_list         = cfg_files;
f_component_list.name    = 'Enter list of component to export (xls)'; 
f_component_list.tag     = 'f_component_list';       %file names
f_component_list.filter  = {'xlsx','xls','txt'};
f_component_list.ufilter = '.*';    
f_component_list.num     = [1 Inf];     % Number of inputs required 
f_component_list.help    = {'The export will be organized using the channel list order.  Each subject must have a close localization on the head to be compared. Enter the list to export columns such as:  ‘NIRS.mat folder’, ‘Type’ (GLM, PARAFAC,PCA), ‘Label’, Component name to filter, ‘Channel List’ full file or .txt same folder as xls file and ‘Name’ use as the output name of the export.'};

E_exportcomponent_list      = cfg_exbranch;
E_exportcomponent_list.name = 'Export list';
E_exportcomponent_list.tag  = 'E_exportcomponent_list';
E_exportcomponent_list.val  = {f_component_list,v_component_list_remove};
E_exportcomponent_list.prog = @nirs_run_E_exportcomponent_list;
E_exportcomponent_list.vout = @nirs_cfg_vout_E_exportcomponent_list  ;
E_exportcomponent_list.help = {'Export component spatial contribution to perform group average or statistics.'};

function vout = nirs_cfg_vout_E_exportcomponent_list(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
f_component_zone         = cfg_files;
f_component_zone.name    = 'Enter list of component to export and there zone (xls)'; 
f_component_zone.tag     = 'f_component_zone';       %file names
f_component_zone.filter  = {'xlsx','xls','txt'};
f_component_zone.ufilter = '.*';    
f_component_zone.num     = [1 Inf];     % Number of inputs required 
f_component_zone.help    = {'Enter the list to export columns such as: ‘NIRS.mat folder’, ‘Type’ (GLM, PARAFAC,PCA),‘Label’ ,Component name to filter, ‘Zone List’ full file or .txt same folder as xls file and ‘Name’ the output name of the export .'}; 

e_zonetoextract_label         = cfg_entry;
e_zonetoextract_label.name    = 'Enter zone to extract';
e_zonetoextract_label.tag     = 'e_zonetoextract_label';       
e_zonetoextract_label.strtype = 's';
e_zonetoextract_label.num     = [1 inf];
e_zonetoextract_label.val     = {'zone1, zone2, zone3 '}; 
e_zonetoextract_label.help    = {'Identify the label to identify the zones to export, separate different labels using a comma.',...
    'As an example: zone1, zone2, zone3'};


f_zonetoextract_file         = cfg_files;
f_zonetoextract_file.name    = 'Enter .txt to list zone to extract'; 
f_zonetoextract_file.tag     = 'f_zonetoextract_file';       %file names
f_zonetoextract_file.filter  = 'txt';
f_zonetoextract_file.ufilter = '.txt';    
f_zonetoextract_file.num     = [1 Inf];     % Number of inputs required 
f_zonetoextract_file.help    = {'Use a text file to identify the zones to export.'}; 


c_zonetoextract         = cfg_choice;
c_zonetoextract.tag     = 'f_zonetoextract_file';
c_zonetoextract.name    = 'Identification zone to export';
c_zonetoextract.values  = {e_zonetoextract_label,f_zonetoextract_file};
c_zonetoextract.val     = {e_zonetoextract_label}; %Default option
c_zonetoextract.help    = {'List label of zones directly or through a .txt file.'};

E_exportcomponent_zone      = cfg_exbranch;
E_exportcomponent_zone.name = 'Export zone';
E_exportcomponent_zone.tag  = 'E_exportcomponent_zone';
E_exportcomponent_zone.val  = {f_component_zone,c_zonetoextract};
E_exportcomponent_zone.prog = @nirs_run_E_exportcomponent_zone;
E_exportcomponent_zone.vout = @nirs_cfg_vout_E_exportcomponent_zone;
E_exportcomponent_zone.help = {'Export component spatial contribution to perform group average or statistics.'};

function vout = nirs_cfg_vout_E_exportcomponent_zone(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


% Folder selector.
e_STATCOMPPath          = cfg_files; %path
e_STATCOMPPath.name     = 'Result folder';
e_STATCOMPPath.tag      = 'e_STATCOMPPath';
e_STATCOMPPath.filter   = {'dir'};
e_STATCOMPPath.ufilter  = '.*';    %
e_STATCOMPPath.num      = [1 1];     % Number of inputs required 
e_STATCOMPPath.help     = {'Result of the statistic will be saved in this folder.'};


m_TtestOneSample        = cfg_menu;
m_TtestOneSample.tag    = 'm_TtestOneSample';
m_TtestOneSample.name   = 'Use ';
m_TtestOneSample.labels = {'2 tails','1 tail left (neg) ' ,'1 tail right (pos) ' };
m_TtestOneSample.values = {1,2,3};
m_TtestOneSample.val    = {1}; 
m_TtestOneSample.help   = {'Use 2 tail or one tail assumption.'};


f_component         = cfg_files;
f_component.name    = 'Enter list of component'; 
f_component.tag     = 'f_component';       %file names
f_component.filter  = 'mat';
f_component.ufilter = '.mat';    
f_component.num     = [1 Inf];     % Number of inputs required 
f_component.help    = {'Enter the list of components to test statistically.'}; 


b_TtestOneSample        = cfg_branch;
b_TtestOneSample.tag    = 'b_TtestOneSample';
b_TtestOneSample.name   = 'One-sample t-test' ;
b_TtestOneSample.val    = {f_component,m_TtestOneSample};
b_TtestOneSample.help   = {'The one-sample t-test is a parametric test of the location parameter when the population standard deviation is unknown. T = mean(X)/(STD(x)*sqrt(n))'};


f_componentG1         = cfg_files;
f_componentG1.name    = 'Group 1 (components export)'; 
f_componentG1.tag     = 'f_componentG1';       %file names
f_componentG1.filter  = 'mat';
f_componentG1.ufilter = '.mat';    
f_componentG1.num     = [1 Inf];     % Number of inputs required 
f_componentG1.help    = {'Enter the list of components to test statistically.'}; 

f_componentG2         = cfg_files;
f_componentG2.name    = 'Group 2 (components export)'; 
f_componentG2.tag     = 'f_componentG2';       %file names
f_componentG2.filter  = 'mat';
f_componentG2.ufilter = '.mat';    
f_componentG2.num     = [1 Inf];     % Number of inputs required 
f_componentG2.help    = {'Enter the list of components to test statistically.',...
    'D1matrix or mean by subject from the export list.'}; 


b_TtestUnpaired        = cfg_branch;
b_TtestUnpaired.tag    = 'b_TtestUnpaired';
b_TtestUnpaired.name   = 'Unpaired t-test' ;
b_TtestUnpaired.val    = {f_componentG1,f_componentG2,m_TtestOneSample};
b_TtestUnpaired.help   = {'The implementation an unpaired t-test.'};

f_anovan         = cfg_files;
f_anovan.name    = 'Group identification'; 
f_anovan.tag     = 'f_anovan';       %file names
f_anovan.filter  = {'xlsx','xls','txt'};
f_anovan.ufilter =  '.*';    
f_anovan.num     = [1 Inf];     % Number of inputs required 
f_anovan.help    = {''}; 

 

b_ANOVAN        = cfg_branch;
b_ANOVAN.tag    = 'b_ANOVAN';
b_ANOVAN.name   = 'By channel' ;
b_ANOVAN.val    = {f_anovan};
b_ANOVAN.help   = {'Apply the anonan on each channel of the xls list. Enter the observation to compute the anovan.',...
    'First row: dir',...
    'Second row: compenent file name',...
    'Third row and following for factor identification.'};


b_ANOVANzone        = cfg_branch;
b_ANOVANzone.tag    = 'b_ANOVANzone';
b_ANOVANzone.name   = 'By zone' ;
b_ANOVANzone.val    = {f_anovan};
b_ANOVANzone.help   = {'Apply the anonan on each specific region define by the xls file. Enter the observation to compute the anovan.',...
    'First row: dir',...
    'Second row: compenent file name',...
    'Third row for the zone file',...
    'Fourth row for the region label to used',...
    'Fifth and following for factor identification '};

c_ANOVAN         = cfg_choice;
c_ANOVAN.tag     = 'c_ANOVAN';
c_ANOVAN.name    = 'Anovan' ;
c_ANOVAN.values  = {b_ANOVAN,b_ANOVANzone};
c_ANOVAN.val     = {b_ANOVAN};
c_ANOVAN.help    = {'N-way analysis of variance, the anova could be performed by channel (channelist) or by zone.'};


c_statcomponent         = cfg_choice; 
c_statcomponent.tag     = 'c_statcomponent';
c_statcomponent.name    = 'Choose the statistical test';
c_statcomponent.values  = {b_TtestOneSample, b_TtestUnpaired,c_ANOVAN};
c_statcomponent.val     = {b_TtestOneSample}; %Default option
c_statcomponent.help    = {'Three statistical tests are available: One sample t-test, Unpaired t-test, and Anovan. The first step is to choose which test you want to perform.'};

E_statcomponent      = cfg_exbranch;
E_statcomponent.name = 'Stats components';
E_statcomponent.tag  = 'E_statcomponent';
E_statcomponent.val  = {c_statcomponent, e_STATCOMPPath};
E_statcomponent.prog = @nirs_run_E_statcomponent;
E_statcomponent.vout = @nirs_cfg_vout_E_statcomponent;
E_statcomponent.help = {'Apply basic statistic on component exported'};

function vout = nirs_cfg_vout_E_statcomponent(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generic data EEG file
EEG_files         = cfg_files;  
EEG_files.name    = 'EEG_files'; % The displayed name
EEG_files.tag     = 'EEG_files';          
EEG_files.num     = [0 Inf];     % Number of inputs required 
EEG_files.val{1}  = {''};
EEG_files.help    = {'Open EEG files.'}; % help text displayed

%%%%%%MODULE 6
% Executable Branch
E_readEEG      = cfg_exbranch;      
E_readEEG.name = 'Read EEG' ;            
E_readEEG.tag  = 'E_readEEG'; 
E_readEEG.val  = {NIRSmat,EEG_files};   
E_readEEG.prog = @nirs_run_readEEG;  
E_readEEG.vout = @nirs_cfg_vout_readEEG;
E_readEEG.help = {'Read simultaneous EEG data; Note trig must be synchronized with the NIRS.',... 
    'Use before the normalization step the block will be cut and associate to the trig in nirs data.',...
    'Use after the normalization step you will need to match the blocks to the nirs data manually.'};


%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_readEEG(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%Generic data EEG file
e_readEEGMarker_Marker         = cfg_entry;
e_readEEGMarker_Marker.name    = 'EEGMarker=NIRSTrigger';
e_readEEGMarker_Marker.tag     = 'e_readEEGMarker_Marker';       
e_readEEGMarker_Marker.strtype = 's';
e_readEEGMarker_Marker.num     = [1 Inf];     
e_readEEGMarker_Marker.val     = {'MARKER1:1,MARKER2:2,MARKER3:3'}; 
e_readEEGMarker_Marker.help    = {'Write the marker to import use the exact label',...
    'As an example: MARKER1:1,MARKER2=2,MARKER3=3'}; 


E_readEEGMarker      = cfg_exbranch;      
E_readEEGMarker.name = 'Import EEG marker' ;            
E_readEEGMarker.tag  = 'E_readEEG'; 
E_readEEGMarker.val  = {NIRSmat, e_readEEGMarker_Marker};   
E_readEEGMarker.prog = @nirs_run_readEEGMarker;  
E_readEEGMarker.vout = @nirs_cfg_vout_readEEGMarker;
E_readEEGMarker.help = {'Import EEG markers as fNIRS triggers. fNIRS uses trigger information (i.e. integer value) to identify the event. The toolbox uses trigger information  segmentation, to model a hemodynamic response or to average fNIRS data.'};


%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_readEEGMarker(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generic data EEG file
Video_files         = cfg_files;  
Video_files.name    = 'Video files'; % The displayed name
Video_files.tag     = 'Video_files';          
Video_files.num     = [0 Inf];     % Number of inputs required 
Video_files.val{1}  = {''};
Video_files.help    = {'Open Video or Audio files.',...
'The issue with 32 bit codec is that it needs 32 bits Matlab version.'}; % help text displayed


video_files_sync         = cfg_menu;
video_files_sync.name    = 'Ref Synchronisation trigger';
video_files_sync.tag     = 'video_files_sync';       
video_files_sync.labels  = {'EEG ','NIRS ','AUX '};
video_files_sync.values  = {0,1,2};
video_files_sync.val     = {0};
video_files_sync.help    = {'Video is recording simultaneously to EEG or NIRS data recording, if they have a delay please edit the video to removing exceeding part at the beginning.'};


b_videooffset_no         = cfg_branch;
b_videooffset_no.tag     = 'b_videooffset_no';
b_videooffset_no.name    = 'No';
b_videooffset_no.val     = {};
b_videooffset_no.help    = {'The video start are synchronise with the trigger'}';

i_videooffset_yes         = cfg_entry;
i_videooffset_yes.name    = 'Lag(seconde) = Video - Ref';
i_videooffset_yes.tag     = 'i_videooffset';       
i_videooffset_yes.strtype = 'r';
i_videooffset_yes.num     = [1 inf];
i_videooffset_yes.val     = {0}; 
i_videooffset_yes.help    = {'As example if the in video the event occur at 5 seconde and the trig in EEG happen at 10 you should adjust the lag to -5 as AUDIO(5)- EEG(10)=LAG(-5)'};;

b_videooffset_yes         = cfg_branch;
b_videooffset_yes.tag     = 'b_videooffset_yes';
b_videooffset_yes.name    = 'Yes';
b_videooffset_yes.val     = {i_videooffset_yes};
b_videooffset_yes.help    = {'They have a delay between the video and the trigger if the video start after use a positive value, if the video start before use a negative value'}';

c_videooffset        = cfg_choice;
c_videooffset.tag    = 'c_videooffset';
c_videooffset.name   = 'Offset';
c_videooffset.values = {b_videooffset_yes, b_videooffset_no};
c_videooffset.val    = {b_videooffset_no};
c_videooffset.help   = {'Is there is lag between Video and Ref Synchronisation trigger, select ''yes''. Otherwise, seelct ''No'''}';

%%%%%%MODULE 6
% Executable Branch
E_readVideo      = cfg_exbranch;      
E_readVideo.name = 'Read Video' ;            
E_readVideo.tag  = 'E_readVideo'; 
E_readVideo.val  = {NIRSmat,Video_files,video_files_sync,c_videooffset};   
E_readVideo.prog = @nirs_run_readVideo;  
E_readVideo.vout = @nirs_cfg_vout_readVideo;
E_readVideo.help = {'Read simultaneous video data, the synchronization could be done using EEG trig, off'};


%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_readVideo(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%Generic data EEG file
Audio_files         = cfg_files;  
Audio_files.name    = 'Audio files'; % The displayed name
Audio_files.tag     = 'Audio_files';          
Audio_files.num     = [0 Inf];     % Number of inputs required 
Audio_files.val{1}  = {''};
Audio_files.help    = {'Open Audio file'}; % help text displayed


Audio_files_sync         = cfg_menu;
Audio_files_sync.name    = 'Ref Synchronisation trigger';
Audio_files_sync.tag     = 'Audio_files_sync';       
Audio_files_sync.labels  = {'EEG ','NIRS ','AUX '};
Audio_files_sync.values  = {0,1,2};
Audio_files_sync.val     = {0};
Audio_files_sync.help    = {'Audio is record silmultaneously to EEG or NIRS data recording, if they have a delay please edit offset'};


b_Audiooffset_no         = cfg_branch;
b_Audiooffset_no.tag     = 'b_Audiooffset_no';
b_Audiooffset_no.name    = 'No';
b_Audiooffset_no.val     = {};
b_Audiooffset_no.help    = {'The Audio start are synchronise with the trigger'}';

i_Audiooffset         = cfg_entry;
i_Audiooffset.name    = 'Lag(seconde) = Audio - Ref';
i_Audiooffset.tag     = 'i_Audiooffset';       
i_Audiooffset.strtype = 'r';
i_Audiooffset.num     = [1 inf];
i_Audiooffset.val     = {0}; 
i_Audiooffset.help    = {'As example if the in audio the event occur at 5 seconds and the trig in EEG happen at 10 you should adjust the lag to -5 as AUDIO(5)- EEG(10)=LAG(-5)'};

b_Audiooffset_yes         = cfg_branch;
b_Audiooffset_yes.tag     = 'b_Audiooffset_yes';
b_Audiooffset_yes.name    = 'Yes';
b_Audiooffset_yes.val     = {i_Audiooffset};
b_Audiooffset_yes.help    = {'They have a delay between the audio and the trigger if the video start after use a positive value, if the audio start before use a negative value'}';

c_Audiooffset        = cfg_choice;
c_Audiooffset.tag    = 'c_Audiooffset';
c_Audiooffset.name   = 'Offset';
c_Audiooffset.values = {b_Audiooffset_yes, b_Audiooffset_no};
c_Audiooffset.val    = {b_Audiooffset_no};
c_Audiooffset.help   = {'Do you have lag between Video and Ref Synchronisation trigger ?'}';


% Executable Branch
E_readAudio      = cfg_exbranch;      
E_readAudio.name = 'Read Audio' ;            
E_readAudio.tag  = 'E_readAudio'; 
E_readAudio.val  = {NIRSmat,Audio_files,Audio_files_sync,c_Audiooffset};   
E_readAudio.prog = @nirs_run_readAudio;  
E_readAudio.vout = @nirs_cfg_vout_readAudio;
E_readAudio.help = {'Read simultaneous audio data, The synchronisation could be done using EEG AUX or NIRS trig, off'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_readAudio(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


%Generic data AUX file
AUX_files         = cfg_files;  
AUX_files.name    = 'AUX_files'; % The displayed name
AUX_files.tag     = 'AUX_files';          
AUX_files.num     = [0 Inf];     % Number of inputs required 
AUX_files.val{1}  = {''};
AUX_files.help    = {'Open AUX files. '}; % help text displayed

%%%%%%MODULE 6
% Executable Branch
E_readAUX      = cfg_exbranch;      
E_readAUX.name = 'Read AUX' ;            
E_readAUX.tag  = 'E_readAUX'; 
E_readAUX.val  = {NIRSmat,AUX_files};   
E_readAUX.prog = @nirs_run_readAUX;  
E_readAUX.vout = @nirs_cfg_vout_readAUX;
E_readAUX.help = {'Read simultaneous AUX data; Note trig must be synchronized with the NIRS',... 
    'Use before the normalization step the block will be cut and associate to the trig in nirs data.',...
    'Use after the normalization step you will need to match the blocks to the nirs data manually.'};



%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_readAUX(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

m_replaceaux          = cfg_menu;
m_replaceaux.tag      = 'm_replaceaux';
m_replaceaux.name     = 'Option';  
m_replaceaux.labels   = {'Replace aux','Add aux'}; 
m_replaceaux.values   = {0,1};
m_replaceaux.val      = {0};
m_replaceaux.help     = {'Replace aux will erase any reference to previous aux file, will add aux will add the new hrf as an additional auxilairy file please do not use the same name if you want to add a new one.'};

e_TimetoPeak1         = cfg_entry; 
e_TimetoPeak1.name    = 'Time to peak of the first gamma density'; 
e_TimetoPeak1.tag     = 'e_TimetoPeak1';      
e_TimetoPeak1.strtype = 'r';       
e_TimetoPeak1.num     = [1 inf];     
e_TimetoPeak1.val     = {[5.4]};
e_TimetoPeak1.help    = {'Time to peak of the first gamma density.'}; 


%General entry use by any modules
e_AUXdir         = cfg_files; 
e_AUXdir.name    = 'Multimodal folder';
e_AUXdir.tag     = 'e_AUXdir';      
e_AUXdir.filter  = 'dir';
e_AUXdir.ufilter = '.*';    
e_AUXdir.dir = '';
e_AUXdir.num = [1 1];  % Number of inputs required 
e_AUXdir.check = [];
e_AUXdir.help    = {'Auxilairy output folder'}; % help text displayed


e_FWHM1         = cfg_entry; 
e_FWHM1.name    = 'FWHM1: approximate FWHM of the first gamma density'; 
e_FWHM1.tag     = 'e_FWHM1';      
e_FWHM1.strtype = 'r';       
e_FWHM1.num     = [1 inf];     
e_FWHM1.val     = {[5.2]};
e_FWHM1.help    = {'FWHM1: approximate FWHM of the first gamma density'}; 

 
e_TimetoPeak2         = cfg_entry; 
e_TimetoPeak2.name    = 'Time to peak of the second gamma density'; 
e_TimetoPeak2.tag     = 'e_TimetoPeak2';      
e_TimetoPeak2.strtype = 'r';       
e_TimetoPeak2.num     = [1 inf];     
e_TimetoPeak2.val     = {[10.8]};
e_TimetoPeak2.help    = {'Time to peak of the second gamma density'}; 


e_FWHM2         = cfg_entry; 
e_FWHM2.name    = 'FWHM2: approximate FWHM of the second gamma density'; 
e_FWHM2.tag     = 'e_FWHM2';      
e_FWHM2.strtype = 'r';       
e_FWHM2.num     = [1 inf];     
e_FWHM2.val     = {[7.35]};
e_FWHM2.help    = {'FWHM2: approximate FWHM of the second gamma density'};

e_DIP         = cfg_entry; 
e_DIP.name    = 'DIP: coefficient of the second gamma density'; 
e_DIP.tag     = 'e_DIP';      
e_DIP.strtype = 'r';       
e_DIP.num     = [1 inf];     
e_DIP.val     = {[0.35]};
e_DIP.help    = {'DIP: coefficient of the second gamma density'}; 


e_HRFlabel         = cfg_entry; 
e_HRFlabel.name    = 'HRF label'; 
e_HRFlabel.tag     = 'e_HRFlabel';      
e_HRFlabel.strtype = 's';       
e_HRFlabel.num     = [1 inf];     
e_HRFlabel.val     = {'HRFtask'};
e_HRFlabel.help    = {'Enter the label identification.'}; 

e_HRFtrigger         = cfg_entry; 
e_HRFtrigger.name    = 'Trigger onset'; 
e_HRFtrigger.tag     = 'e_HRFtrigger';      
e_HRFtrigger.strtype = 'r';       
e_HRFtrigger.num     = [1 inf];     
e_HRFtrigger.val     = {[2 4]};
e_HRFtrigger.help    = {'Trigger onset'}; 

 
e_HRFduration         = cfg_entry; 
e_HRFduration.name    = 'Duration (s)'; 
e_HRFduration.tag     = 'e_HRFduration';      
e_HRFduration.strtype = 'r';       
e_HRFduration.num     = [1 inf];     
e_HRFduration.val     = {[30]};
e_HRFduration.help    = {'Trigger onset'};


b_HRFtriggeronset        = cfg_branch;
b_HRFtriggeronset.tag    = 'b_HRFtriggeronset';
b_HRFtriggeronset.name   = 'HRF trigger onset';
b_HRFtriggeronset.val    = {e_HRFtrigger, e_HRFduration e_HRFlabel e_TimetoPeak1, e_FWHM1,e_TimetoPeak2,e_FWHM2,e_DIP,e_AUXdir };
b_HRFtriggeronset.help   = {'Model HRF response using trigger onset and user define fix duration.'};

c_createAUXauto         = cfg_choice;
c_createAUXauto.tag     = 'c_createAUXauto';
c_createAUXauto.name    = 'Choose AUX to generate';
c_createAUXauto.values  = {b_HRFtriggeronset};
c_createAUXauto.val     = {b_HRFtriggeronset}; %Default option
c_createAUXauto.help    = {''};

% Executable Branch
E_createAUXauto      = cfg_exbranch;      
E_createAUXauto.name = 'Create AUX' ;            
E_createAUXauto.tag  = 'E_createAUXauto'; 
E_createAUXauto.val  = {NIRSmat,m_replaceaux,c_createAUXauto};   
E_createAUXauto.prog = @nirs_run_createAUXauto;  
E_createAUXauto.vout = @nirs_cfg_vout_createAUXauto;
E_createAUXauto.help = {'Read simultaneous AUX data; Note trig must be synchronized with the NIRS',... 
    'Use before the normalization step the block will be cut and associate to the trig in nirs data.',...
    'Use after the normalization step you will need to match the blocks to the nirs data manually.'};



%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_createAUXauto(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



% Executable Branch
E_createAUX      = cfg_exbranch;      
E_createAUX.name = 'Use GUI_AUXedit' ;            
E_createAUX.tag  = 'E_createAUX'; 
E_createAUX.val  = {};   
E_createAUX.prog = @nirs_run_GUIAUXEDIT;  
E_createAUX.help = {'This is a utilitary GUI to create auxiliary information for GLM regression or visualization.'};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executable Branch
I_zonecorrlist         = cfg_files; %Select NIRS.mat for this subject 
I_zonecorrlist.name    = 'List of zone node'; % The displayed name
I_zonecorrlist.tag     = 'I_zonecorrlist';       %file names
I_zonecorrlist.val     = {};
I_zonecorrlist.filter  = 'zone';
I_zonecorrlist.ufilter = '.zone';    
I_zonecorrlist.num     = [1 Inf];     % Number of inputs required 
I_zonecorrlist.help    = {'Use a zone structure that identifies channels belonging to a region of interest. This structure could be created and saved in the DisplayGUI.'};


I_chcorrlist         = cfg_files; %Select NIRS.mat for this subject 
I_chcorrlist.name    = 'List of channels as nodes'; % The displayed name
I_chcorrlist.tag     = 'I_chcorrlist';       %file names
I_chcorrlist.val     = {};
I_chcorrlist.filter  = 'txt';
I_chcorrlist.ufilter = '.txt';    
I_chcorrlist.num     = [1 Inf];     % Number of inputs required 
I_chcorrlist.help    = {'List of channels to be used as a node using a text file, one node per line, example: F a1b2;',...
    'To create automatically this list use the module Utility /Create channel list.'};

b_nodelist         =  cfg_choice; %Select NIRS.mat for this subject 
b_nodelist.name    = 'Node list'; % The displayed name
b_nodelist.tag     = 'b_nodelist';       %file names
b_nodelist.val     = {I_chcorrlist };
b_nodelist.values  = {I_chcorrlist I_zonecorrlist};  
b_nodelist.help    = {'Nodes could be defined by a region of interest (zones) or each channel. Technical note: to allow subject comparison to ensure channels or zones spatial localization analogous.'};



%Input Frequency
i_Freq_crossspectrum         = cfg_entry; 
i_Freq_crossspectrum.name    = 'Frequency range to obtain Cxy(f)'; 
i_Freq_crossspectrum.tag     = 'i_Freq_crossspectrum';      
i_Freq_crossspectrum.strtype = 'r';       
i_Freq_crossspectrum.num     = [1 inf];     
i_Freq_crossspectrum.val     = {[0.03,0.08]};
i_Freq_crossspectrum.help    = {'Frequency to extract the matrix output.'}; 

%Zone channel plot
i_ch_crossspectrum         = cfg_files; 
i_ch_crossspectrum.name    = 'Zone to display spectrum'; 
i_ch_crossspectrum.tag     = 'i_ch_crossspectrum';      
i_ch_crossspectrum.filter  = 'zone';
i_ch_crossspectrum.ufilter = '.zone';          
i_ch_crossspectrum.num     = [1 inf];     
i_ch_crossspectrum.val     = {};
i_ch_crossspectrum.help    = {'Plot and save the spectrum in for zone, it will be saved in the same folder as the matrix.'}; 


%Input Frequency
i_TrialLenght_crossspectrum         = cfg_entry; 
i_TrialLenght_crossspectrum.name    = 'Trial length (s)'; 
i_TrialLenght_crossspectrum.tag     = 'i_TrialLenght_crossspectrum';      
i_TrialLenght_crossspectrum.strtype = 'r';       
i_TrialLenght_crossspectrum.num     = [1 inf];     
i_TrialLenght_crossspectrum.val     = {60};
i_TrialLenght_crossspectrum.help    = {'Length of the segments adjust to change FFT frequency measured.'}; 

%Input Frequency
i_RandomSample_crossspectrum         = cfg_entry; 
i_RandomSample_crossspectrum.name    = 'Number of random sample'; 
i_RandomSample_crossspectrum.tag     = 'i_RandomSample_crossspectrum';      
i_RandomSample_crossspectrum.strtype = 'r';       
i_RandomSample_crossspectrum.num     = [1 inf];     
i_RandomSample_crossspectrum.val     = {200};
i_RandomSample_crossspectrum.help    = {'Number of random samples to compute the cross spectrum, they are take using pseudorandom segment in the whole block.'}; 

%Input Frequency
i_OutlierControl_crossspectrum         = cfg_entry; 
i_OutlierControl_crossspectrum.name    = 'Zscore Outlier Control'; 
i_OutlierControl_crossspectrum.tag     = 'i_OutlierControl_crossspectrum';      
i_OutlierControl_crossspectrum.strtype = 'r';       
i_OutlierControl_crossspectrum.num     = [1 inf];     
i_OutlierControl_crossspectrum.val     = {5};
i_OutlierControl_crossspectrum.help    = {'Zscore in the frequency domain specific to each channel and frequency help to remove outlier.'}; 
 
m_savefft_crossspectrum          = cfg_menu;
m_savefft_crossspectrum.tag      = 'm_savefft_crossspectrum';
m_savefft_crossspectrum.name     = 'Save FFT spectrum';  
m_savefft_crossspectrum.labels   = {'No','Yes'}; 
m_savefft_crossspectrum.values   = {0,1};
m_savefft_crossspectrum.val      = {0};
m_savefft_crossspectrum.help     = {'Choose whether the spectrum is saved or not.'};

b_crossspectrum        = cfg_branch;
b_crossspectrum.tag    = 'b_crossspectrum';
b_crossspectrum.name   = 'Coherence';
b_crossspectrum.val    = {i_Freq_crossspectrum,i_ch_crossspectrum,i_TrialLenght_crossspectrum,i_RandomSample_crossspectrum,i_OutlierControl_crossspectrum,m_savefft_crossspectrum};
b_crossspectrum.help   = {'Coherence is a statistic representing the relationship between two signals and is also an extension of correlation to the frequency domain (Kida, 2016). Coherence is known as magnitude squared coherence is defined as the complex conjugate product of the Fourier transforms data X(f)* Y*T(f). x(t) and y(t) are two time series, Gxy(f) is the cross-spectral density between x and y, and Gxx(f) and Gyy(f) are the auto spectral densities of x and y, respectively. The coherence is implemented to use one long continuous segment of the recording. In case you record multiple sessions you may join them using the concatenate module. A large number of segments (Number of random samples) of duration (Length of the segment) will be picked randomly (circular bootstrap). Any segments that belong to a specific artifact period will be excluded from the coherence calculation. The segment will be randomly segmented to calculate coherence based on many segments. An FFT is computed on each random segment and the coherence is measured based on the specified frequency range (The frequency range to obtain Cxy(f)) to obtain a connectivity matrix representative of the whole recording. '};


estd_Phase         = cfg_entry; %path
estd_Phase.name    = 'Threshold phase standard deviation';
estd_Phase.tag     = 'estd_Phase';       
estd_Phase.strtype = 's';
estd_Phase.num     = [1 Inf];     
estd_Phase.val     = {'50'}; 
estd_Phase.help    = {['Reject phase with standard deviation higher than ']}; 

b_Phase        = cfg_branch;
b_Phase.tag    = 'b_Phase';
b_Phase.name   = 'Phase joint probability ISS raw';
b_Phase.val    = {estd_Phase };
b_Phase.help   = {'Look at the joint phase probability of the ISS phase measure (unwrap and ppod as previous step) circ_kurtosis,',...
    'References:Pewsey, Metrika, 2004 Fisher, Circular Statistics, p. 34 Circular Statistics Toolbox for Matlab By Philipp Berens, 2009',...
    'berens@tuebingen.mpg.de'};

edownsample_Granger         = cfg_entry; %path
edownsample_Granger.name    = 'Downsample';
edownsample_Granger.tag     = 'edownsample_Granger';       
edownsample_Granger.strtype = 's';
edownsample_Granger.num     = [1 Inf];     
edownsample_Granger.val     = {'10'}; 
edownsample_Granger.help    = {['Downsample factor usefull if you want to do granger causality efficiently on slow frequency signal ']}; 


enb_Granger         = cfg_entry; %path
enb_Granger.name    = 'Lag order';
enb_Granger.tag     = 'enb_Granger';       
enb_Granger.strtype = 's';
enb_Granger.num     = [1 Inf];     
enb_Granger.val     = {'3'}; 
enb_Granger.help    = {['Nb of compendant for the Multivariate analysis']}; 

b_Granger        = cfg_branch;
b_Granger.tag    = 'b_Granger';
b_Granger.name   = 'Multivariate Granger Causality';
b_Granger.val    = {enb_Granger, edownsample_Granger};
b_Granger.help   = {'Perform time domain Multivariate granger causality, Use MVGC Toolbox',...
                    'this toolbox are Copyright (C) Lionel Barnett and Anil K. Seth, 2012.',...
                    'Please include in Matlab set path the mvgc_v1.0 toolbox,it cause some conflic with other built in matlab function,',...
                    'then we recomand to remove the folder from the set path after used to avoid potential conflics'};
            
                     
b_PearsonBootstrap         = cfg_branch;
b_PearsonBootstrap.tag     = 'b_PearsonBootstrap';
b_PearsonBootstrap.name    = 'Circular bootstrap';
b_PearsonBootstrap.val     = {i_TrialLenght_crossspectrum,i_RandomSample_crossspectrum,i_OutlierControl_crossspectrum};
b_PearsonBootstrap.help    = {'Use circular bootstrap to compute cross-corelation analysis.'};
                
m_Pearson           = cfg_menu;
m_Pearson.tag       = 'm_Pearson';
m_Pearson.name      = 'Segments';
m_Pearson.labels    = {'Zero lag Cross-Correlation (Pearson)'};
m_Pearson.values    = {1};
m_Pearson.val       = {1};
m_Pearson.help      = {'Apply the zero lag Cross-Correlation on each segmented file of the data set'};


c_Pearson          = cfg_choice;
c_Pearson.tag     = 'c_Pearson';
c_Pearson.name    = 'Option';
c_Pearson.values  = {m_Pearson,b_PearsonBootstrap};
c_Pearson.val     = {m_Pearson}; %Default option
c_Pearson.help    = {''};

b_Pearson       = cfg_branch;
b_Pearson.tag    = 'b_Pearson';
b_Pearson.name   = 'Cross-correlation analysis';
b_Pearson.val    = {c_Pearson};
b_Pearson.help   = {'Compute zeros lag cross-correlation between channel or region of interest.'};                  



m_Hilbert           = cfg_menu;
m_Hilbert.tag       = 'm_Hilbert';
m_Hilbert.name      = 'Use file segment';
m_Hilbert.labels    = {'Hilbert transform'};
m_Hilbert.values    = {1};
m_Hilbert.val       = {1};
m_Hilbert.help      = {'Use each data file segmented to compute the Hilbert transform.'};



b_HilbertBootstrap         = cfg_branch;
b_HilbertBootstrap.tag     = 'b_HilbertBootstrap';
b_HilbertBootstrap.name    = 'Circular bootstrap';
b_HilbertBootstrap.val     = {i_TrialLenght_crossspectrum,i_RandomSample_crossspectrum,i_OutlierControl_crossspectrum};
b_HilbertBootstrap.help    = {'Use circular bootstrap to segment randomly multiple segment in the fNIRS files to compute the Hilbert transform.'};



c_Hilbert         = cfg_choice;
c_Hilbert.tag     = 'c_Hilbert';
c_Hilbert.name    = 'Option';
c_Hilbert.values  = {m_Hilbert,b_HilbertBootstrap};
c_Hilbert.val     = {m_Hilbert}; %Default option
c_Hilbert.help    = {''};

b_Hilbert        = cfg_branch;
b_Hilbert.tag    = 'b_Hilbert';
b_Hilbert.name   = 'Hilbert join phase probability';
b_Hilbert.val    = {c_Hilbert};
b_Hilbert.help   = {'See B. Molavi et al., “Analyzing the resting state functional connectivity in the human language system using near infrared spectroscopy.”,',...
                    'Front. Hum. Neurosci. 7(January), 921 (2013) [doi:10.3389/fnhum.2013.00921].',...
                    'Use fonction from the CircStat2012a toolbox' };



e_startwaveletcluster         = cfg_entry; %path
e_startwaveletcluster.name    = 'Start Time (s)';
e_startwaveletcluster.tag     = 'e_startwaveletcluster';       
e_startwaveletcluster.strtype = 's';
e_startwaveletcluster.num     = [1 Inf];     
e_startwaveletcluster.val     = {'100'}; 
e_startwaveletcluster.help    = {['Start time for wavelet decomposition']}; 

e_stopwaveletcluster         = cfg_entry; %path
e_stopwaveletcluster.name    = 'Stop Time (s)';
e_stopwaveletcluster.tag     = 'e_stopwaveletcluster';       
e_stopwaveletcluster.strtype = 's';
e_stopwaveletcluster.num     = [1 Inf];     
e_stopwaveletcluster.val     = {'300'}; 
e_stopwaveletcluster.help    = {['Stop time for wavelet decomposition']}; 

e_widthwaveletcluster         = cfg_entry; %path
e_widthwaveletcluster.name    = 'Width morlet coefficient';
e_widthwaveletcluster.tag     = 'e_widthwaveletcluster';       
e_widthwaveletcluster.strtype = 's';
e_widthwaveletcluster.num     = [1 Inf];     
e_widthwaveletcluster.val     = {'3'}; 
e_widthwaveletcluster.help    = {['Width parameter morlet coefficient']}; 

b_waveletcluster        = cfg_branch;
b_waveletcluster.tag    = 'b_waveletcluster';
b_waveletcluster.name   = 'Wavelet Coherence';
b_waveletcluster.val    = {e_startwaveletcluster,e_stopwaveletcluster,i_Freq_crossspectrum, e_widthwaveletcluster };
b_waveletcluster.help   = {''};


I_chcorrlist_type         = cfg_choice;
I_chcorrlist_type.tag     = 'I_chcorrlist_type';
I_chcorrlist_type.name    = 'Connectivity to use';
I_chcorrlist_type.val     = {b_crossspectrum};
I_chcorrlist_type.help    = {''};
I_chcorrlist_type.values  = {b_Pearson, b_Hilbert,b_crossspectrum};%,b_Granger,b_Phase,b_waveletcluster};%,b_Hilbert,b_Granger, b_Phase, b_crossspectrum
% 
% I_chcorrlist_type.labels = {'Pearson', 'Pearson with zscore','Hilbert phase joint probability', 'Granger','Phase ISS','Analyzer Correlation/Autocorrelation'};
% I_chcorrlist_type.values = {1 2 3 4 5 6};
% I_chcorrlist_type.def  = {1};
% I_chcorrlist_type.help = {'Use pearson connectivity matrix, an average on all blocks'};
%     
I_chcorrlistoutpath         = cfg_entry; %path
I_chcorrlistoutpath.name    = 'Path connectivity matrix';
I_chcorrlistoutpath.tag     = 'I_chcorrlistoutpath';       
I_chcorrlistoutpath.strtype = 's';
I_chcorrlistoutpath.num     = [1 Inf];     
I_chcorrlistoutpath.val     = {'C:\data\Data_NIRS\BebeResting\connectivityMAT\'}; 
I_chcorrlistoutpath.help    = {['Path for output files.']}; 


I_ConnectivityMATName         = cfg_entry;
I_ConnectivityMATName.name    = 'File name';
I_ConnectivityMATName.tag     = 'I_ConnectivityMATName';       
I_ConnectivityMATName.strtype = 's';
I_ConnectivityMATName.num     = [0 inf];
I_ConnectivityMATName.val     = {};
I_ConnectivityMATName.help    = {'OutputName for the matrix. Use identifier as subject number.'};


E_chcorrMat      = cfg_exbranch;      
E_chcorrMat.name = 'Connectivity matrix using channel list (NIRS.mat)' ;            
E_chcorrMat.tag  = 'E_chcorrMat'; 
E_chcorrMat.val  = {NIRSmat, b_nodelist I_chcorrlist_type I_chcorrlistoutpath I_ConnectivityMATName};   
E_chcorrMat.prog = @nirs_run_chcorrMat;  
E_chcorrMat.vout = @nirs_cfg_vout_chcorrMat;
E_chcorrMat.help = {'Create a connectivity matrix using channels as node information or zones as node information.'};
%make NIRS.mat available as a dependency


function vout = nirs_cfg_vout_chcorrMat(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

% Executable Branch
E_createseedlist      = cfg_exbranch;      
E_createseedlist.name = 'Create channel list' ;            
E_createseedlist.tag  = 'E_createseedlist'; 
E_createseedlist.val  = {NIRSmat};   
E_createseedlist.prog = @nirs_run_createseedlist;  
E_createseedlist.vout = @nirs_cfg_vout_createseedlist;
E_createseedlist.help = {'Create a list of channels present in the NIRS file; this list could be used as seed for this helmet.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_createseedlist(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

f_qualityreport          = cfg_files; %path
f_qualityreport.name     = 'Result folder';
f_qualityreport.tag      = 'f_qualityreport';
f_qualityreport.filter   = {'dir'};
f_qualityreport.ufilter  = '.*';    %
f_qualityreport.num      = [1 1];     % Number of inputs required 
f_qualityreport.help     = {'QualityReport.xlsx result will be saved in this folder.'};


% Executable Branch
E_qualityreport      = cfg_exbranch;      
E_qualityreport.name = 'Quality report' ;            
E_qualityreport.tag  = 'E_qualityreport'; 
E_qualityreport.val  = {NIRSmat, f_qualityreport};   
E_qualityreport.prog = @nirs_run_qualityreport;  
E_qualityreport.vout = @nirs_cfg_vout_qualityreport;
E_qualityreport.help = {'Create a small report of the time of data corrected in the signal.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_qualityreport(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

m_zonematformat           = cfg_menu;
m_zonematformat.tag       = 'm_zonematformat';
m_zonematformat.name      = 'Label format';
m_zonematformat.labels    = {'NIRx', 'ISS'};
m_zonematformat.values    = {1,2};
m_zonematformat.val       = {1};
m_zonematformat.help      = {'Write the list of detector source according to NIRx label D01 S1 for detector D01 and source S1 or according to ISS label A a1b2 for detector A and source a1b2.'};

%General entry use by any modules
zonemat         = cfg_files; 
zonemat.name    = 'zone.mat';
zonemat.tag     = 'zonemat';      
zonemat.filter  = 'zone';
zonemat.ufilter = '.zone$';    
zonemat.num     = [1 Inf];     % Number of inputs required 
zonemat.help    = {'Select zone to create channel list for verifaction, zone label as well as detector source combinaison will be list'}; % help text displayed

% Executable Branch
E_zone2channellist      = cfg_exbranch;      
E_zone2channellist.name = 'Transfert .zone into a .txt list' ;            
E_zone2channellist.tag  = 'E_zone2channellist'; 
E_zone2channellist.val  = {zonemat, m_zonematformat};   
E_zone2channellist.prog = @nirs_run_E_zone2channellist;  
E_zone2channellist.help = {'Create a list with the labels, RBG color and the channel of the zone.'};

%General entry use by any modules
zonetxt         = cfg_files; 
zonetxt.name    = 'zone.txt';
zonetxt.tag     = 'zonetxt';      
zonetxt.filter  = 'txt';
zonetxt.ufilter = '.txt$';    
zonetxt.num     = [1 Inf];     % Number of inputs required 
zonetxt.help    = {'Select .txt list that contain zonelabel  to create channel list for verifaction, zone label as well as detector source combinaison will be list'}; % help text displayed

% Executable Branch
E_channellist2zone      = cfg_exbranch;      
E_channellist2zone.name = 'Transfert a .txt list into .zone' ;            
E_channellist2zone.tag  = 'E_channellist2zone'; 
E_channellist2zone.val  = {NIRSmat,zonetxt, m_zonematformat};   
E_channellist2zone.prog = @nirs_run_E_channellist2zone;  
E_channellist2zone.vout = @nirs_cfg_vout_channellist2zone;
E_channellist2zone.help = {'Create a zone using list of detectors and sources belonging to this zone, use Label: zonename ;  RGBcolor: 255 0 0 ;  then a line for each Detector Source combinaison D01 E1.',...
    'The zone will be create at the NIRS.mat location to be associate to this subject.'};

%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_channellist2zone(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



E_GUI_lookmatrices = cfg_exbranch;      
E_GUI_lookmatrices.name = 'Matrix visualisation GUI' ;            
E_GUI_lookmatrices.tag  = 'E_GUI_lookmatrices'; 
E_GUI_lookmatrices.val  = {};   
E_GUI_lookmatrices.prog = @nirs_run_E_GUI_lookmatrices;  
E_GUI_lookmatrices.help = {'Open a GUI to visualized fNIRS connectivity matrix.'};


m_fishertransform           = cfg_menu;
m_fishertransform.tag       = 'm_fishertransform';
m_fishertransform.name      = 'Fisher transform';
m_fishertransform.labels    = {'Yes', 'No'};
m_fishertransform.values    = {1,2};
m_fishertransform.val       = {1};
m_fishertransform.help      = {'Use the fisher transform 1/2 ln((1+p)/(1-p)), when the transformation is applied to the sample correlation coefficient, the sampling distribution of the resulting variable is approximately normal, with a variance that is stable over different values of the underlying true correlation.'};


m_nodeunit           = cfg_menu;
m_nodeunit.tag       = 'm_nodeunit';
m_nodeunit.name      = 'Nodes';
m_nodeunit.labels    = {'Channels', 'Zones'};
m_nodeunit.values    = {1,2};
m_nodeunit.val       = {1};
m_nodeunit.help      = {'Apply the statistic on each node. Define nodes as each channel or each average of zone channels.'};

m_TtestOneSample_matrix        = cfg_menu;
m_TtestOneSample_matrix.tag    = 'm_TtestOneSample_matrix';
m_TtestOneSample_matrix.name   = 'Use ';
m_TtestOneSample_matrix.labels = {'2 tails','1 tail left (neg)' ,'1 tail right (pos)' };
m_TtestOneSample_matrix.values = {1,2,3};
m_TtestOneSample_matrix.val    = {3}; 
m_TtestOneSample_matrix.help   = {'Use one tail or two tail t-test.'};

f_matrix         = cfg_files;
f_matrix.name    = 'Enter list connectivity matrix'; 
f_matrix.tag     = 'f_matrix';       %file names
f_matrix.filter  = {'xlsx','xls','txt'};
f_matrix.ufilter = '.*';    
f_matrix.num     = [1 Inf];     % Number of inputs required 
f_matrix.help    = {'Enter the list of connectivity matrix to test statistically.',...
    'column 1 dir, column 2 name,  column 3 zone id, column 4 group 1 include 0 exclude'}; 

e_TtestOneSampleGR         = cfg_entry; %path
e_TtestOneSampleGR.name    = 'Group identification';
e_TtestOneSampleGR.tag     = 'e_TtestOneSampleGR';       
e_TtestOneSampleGR.strtype = 'r';
e_TtestOneSampleGR.num     = [1 Inf];     
e_TtestOneSampleGR.val     = {1}; 
e_TtestOneSampleGR.help    = {['Enter the group to apply the one sample t-test (refer to groupe value in the xls file).']}; 

b_TtestOneSamplematrix        = cfg_branch;
b_TtestOneSamplematrix.tag    = 'b_TtestOneSamplematrix';
b_TtestOneSamplematrix.name   = 'One-sample t-test' ;
b_TtestOneSamplematrix.val    = {m_TtestOneSample_matrix,e_TtestOneSampleGR};
b_TtestOneSamplematrix.help   = {'The one-sample t-test is a parametric test of the location parameter when the population standard deviation is unknown.',...
    'T = mean(X)/(STD(x)*sqrt(n))'};

e_TtestOneSampleGR2         = cfg_entry; %path
e_TtestOneSampleGR2.name    = 'Group 2 identification';
e_TtestOneSampleGR2.tag     = 'e_TtestOneSampleGR2';       
e_TtestOneSampleGR2.strtype = 'r';
e_TtestOneSampleGR2.num     = [1 Inf];     
e_TtestOneSampleGR2.val     = {2}; 
e_TtestOneSampleGR2.help    = {['Enter the second group compared (refer to groupe value in the xls file)']}; 

e_npermutation         =  cfg_entry;
e_npermutation.name    = 'Nb permutations';
e_npermutation.tag     = 'e_npermutation';       
e_npermutation.strtype = 's';
e_npermutation.num     = [0 inf];
e_npermutation.val     = {'500'};
e_npermutation.help    = {'Define the number of non parametric permutation to used'};


b_PermutationTest = cfg_branch;
b_PermutationTest.tag    = 'b_PermutationTest';
b_PermutationTest.name   = 'Unpaired permutation t-test' ;
b_PermutationTest.val    = {e_npermutation, e_TtestOneSampleGR, e_TtestOneSampleGR2};
b_PermutationTest.help   = {'Compared 2 groups using permutation '};

b_exportNBSformat        = cfg_branch;
b_exportNBSformat.tag    = 'b_exportNBSformat';
b_exportNBSformat.name   = 'Export NBS format' ;
b_exportNBSformat.val    = {};
b_exportNBSformat.help   = {'Export to NBS network based statistic format. A Zalesky(2010) [doi: 10.1016/j.neuroimage.2010.06.041] '};


b_Covariable_Mat         =  cfg_entry;
b_Covariable_Mat.name    = 'Covariable';
b_Covariable_Mat.tag     = 'b_Covariable_Mat';       
b_Covariable_Mat.strtype = 's';
b_Covariable_Mat.num     = [0 inf];
b_Covariable_Mat.val     = {'Name column'};
b_Covariable_Mat.help    = {'Use the exact column title to reconized which covariable to correlate with the connectivity scores. Separate by a comma when they are many covariables to explore'};


b_PearsonCorr_Mat        = cfg_branch;
b_PearsonCorr_Mat.tag    = 'b_PearsonCorr_Mat';
b_PearsonCorr_Mat.name   = 'Pearson correlation' ;
b_PearsonCorr_Mat.val    = {b_Covariable_Mat};
b_PearsonCorr_Mat.help   = {'Compute the pearson correlation coefficient between connectivity score and one covariable in the excel file'};

b_GLM_Mat        = cfg_branch;
b_GLM_Mat.tag    = 'b_GLM_Mat';
b_GLM_Mat.name   = 'GLM' ;
b_GLM_Mat.val    = {b_Covariable_Mat};
b_GLM_Mat.help   = {'Apply the general linear model, specify covariable as regressor y=b1*x1+b2*x2+...+c, do not forget to include a covariable for the constante'};


c_statmatrix         = cfg_choice;
c_statmatrix.tag     = 'c_statmatrix';
c_statmatrix.name    = 'Choose the statistical test';
c_statmatrix.values  = {b_TtestOneSamplematrix,b_PermutationTest,b_PearsonCorr_Mat, b_GLM_Mat, b_exportNBSformat};
c_statmatrix.val     = {b_TtestOneSamplematrix}; %Default option
c_statmatrix.help    = {''};

% Folder selector.
e_statmatrixPath          = cfg_files; %path
e_statmatrixPath.name     = 'Result folder';
e_statmatrixPath.tag      = 'e_statmatrixPath';
e_statmatrixPath.filter   = {'dir'};
e_statmatrixPath.ufilter  = '.*';    %
e_statmatrixPath.num      = [1 1];     % Number of inputs required 
e_statmatrixPath.help     = {'Result of the statistics will be saved in this folder.'};

E_statmatrix    = cfg_exbranch;
E_statmatrix.name = 'Stats Matrices';
E_statmatrix.tag  = 'E_statmatrix';
E_statmatrix.val  = {f_matrix,m_fishertransform,m_nodeunit,c_statmatrix, e_statmatrixPath};
E_statmatrix.prog = @nirs_run_E_statmatrix;
E_statmatrix.vout = @nirs_cfg_vout_E_statmatrix;
E_statmatrix.help = {'Apply basic statistic on connectivity matrices exported'};

function vout = nirs_cfg_vout_E_statmatrix(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module IOMTG (creation .prj)
M_createMTG        = cfg_choice; 
M_createMTG.name   = '3dMTG';
M_createMTG.tag    = 'M_createMTG';
M_createMTG.values = {E_IOMTG}; 
M_createMTG.help   = {'Use to create the montage configuration file.'};

M_readMultimodal        = cfg_choice; 
M_readMultimodal.name   = 'Read Multimodal';
M_readMultimodal.tag    = 'M_readMultimodal';
M_readMultimodal.values = {E_readEEG E_readEEGMarker E_readAUX E_createAUXauto E_readVideo E_readAudio E_createAUX}; %
M_readMultimodal.help   = {'Read EEG or Auxiliary format to do visualised multimodal NIRS in the GUI'};

%Module Read data 
M_readNIRS        = cfg_choice; 
M_readNIRS.name   = 'Read data';
M_readNIRS.tag    = 'M_readNIRS';
M_readNIRS.values = {E_readNIRxscout,E_rawhomer, E_readSNIRF,boxy1,  M_readMultimodal}; 
M_readNIRS.help   = {'These modules read NIRS data in different formats.'};

%Module segment or concatenate data 
M_Segment        =   cfg_choice; 
M_Segment.name   =  'Segment/Onset';
M_Segment.tag    = 'M_Segment';
M_Segment.values = {E_segment ,E_Concatenate_file,E_Concatenate_nirsmat,E_aux2manualtrig,E_manualtrig}; 
M_Segment.help   = {'These modules segment or combine data.'};



%Module 3 Preprocessing NIRS data (3a,3b,3c,3d,3e,3f)
M_preprocessing        = cfg_choice; 
M_preprocessing.name   = 'Preprocessing NIRS data';
M_preprocessing.tag    = 'M_preprocessing'; 
M_preprocessing.values = {E_artefactdetection E_chcardiaccontrol E_normalization E_filter ODtoHbOHbR  E_detrend E_average E_nullifybad}; %preprocess enlever fast preprocessing JT
M_preprocessing.help   = {'These apply basic operation on fNIRS data'};



%Module 6 Write external file
M_datawritenirs        =  cfg_choice; 
M_datawritenirs.name   = 'Write external file';
M_datawritenirs.tag    = 'M_datawritenirs';
M_datawritenirs.values = {E_writeNIRS, E_writeHMR, E_writeSNIRF, E_NIR_segment }; 
M_datawritenirs.help   = {'These module convers nir file in .nirs'};

%Module 7 Data Display GUI'
M_GUI        = cfg_choice;
M_GUI.name   = 'Data Display';
M_GUI.tag    = 'M_GUI';
M_GUI.values = {E_GUI E_VIDEO }; 
M_GUI.help   = {'Graphical user interface for data display.'};

%Module 8 Connectivity

%Module 9 Component
M_dataComponent = cfg_choice;
M_dataComponent.name   = 'Decomposition';
M_dataComponent.tag    = 'M_dataComponent';
M_dataComponent.values = { E_extractcomponent E_substractcomponent E_exportcomponent_list E_exportcomponent_zone E_statcomponent}; 
M_dataComponent.help   = {'Apply operation on selected component identify in display GUI data decomposition'};


%Module Connectivity
M_Connectivity        = cfg_choice; 
M_Connectivity.name   = 'Connectivity';
M_Connectivity.tag    = 'M_Connectivity';
M_Connectivity.values = { E_chcorrMat E_GUI_lookmatrices E_statmatrix}; %
M_Connectivity.help   = {'Connectivity fonction.'};

%Module Utility
M_Utility        = cfg_choice; 
M_Utility.name   = 'Utility NIRSmat';
M_Utility.tag    = 'M_Utility';
M_Utility.values = {E_NIRSmatdiradjust, E_NIRSmatcreatenewbranch  E_createseedlist E_qualityreport E_zone2channellist E_channellist2zone E_viewNIRS M_datawritenirs}; %
M_Utility.help   = {'Utility on NIRSmat fonction.'};

nirsHSJ        = cfg_choice;
nirsHSJ.name   = 'LIONirs';
nirsHSJ.tag    = 'nirsHSJ'; %Careful, this tag nirsHSJ must be the same as
%the name of the toolbox and when called by spm_jobman in nirsHSJ.mM_dataprocessing 
nirsHSJ.values = {M_createMTG M_readNIRS M_Segment  M_preprocessing  M_dataComponent  M_Connectivity M_Utility E_GUI}; 

end
