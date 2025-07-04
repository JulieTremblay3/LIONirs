% Toolbox lionirs main function menu
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine 
% www.lionlab.umontreal.ca
%___________________________________________________________________
function nirsHSJ = tbx_cfg_LIONirs;


addpath(fileparts(which(mfilename))); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Module 0 Montage project creation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_IOMTG         = cfg_exbranch;
E_IOMTG.name    = 'Montage configuration (.prj)';
E_IOMTG.tag     = 'E_IOMTG';
E_IOMTG.val     = {};
E_IOMTG.prog    = @nirs_run_E_IOMTG;
E_IOMTG.help    = {'Use this module to create project configuration file using source and detector position.'};


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
prjfile.help    = {'Associate recording helmet configuration .prj created in 3DMTG.'}; 


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
age1.help    = {'Write the subject�s age in years.',...
    'Used to calculate DPF to HbO/HbR conversion (Duncan et al 1995).',...
    'Only for wavelengths of 690,830,744,807 nm.'};

raw_onset_files         = cfg_files;  
raw_onset_files.name    = 'Select onset files'; % The displayed name
raw_onset_files.tag     = 'raw_onset_files';          
raw_onset_files.num     = [0 Inf];     % Number of inputs required 
raw_onset_files.val{1}  = {''};
raw_onset_files.help    = {'Optional: Select raw onset files. '
    'Can be added at a later stage.'
    'Must specify one file for each data file, in the same order.'}'; % help text displayed


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
output_path.name    = 'Analyzed files directory';
output_path.tag     = 'output_path';       
output_path.strtype = 's';
output_path.num     = [1 Inf];     
output_path.val     = {'c:/Data/Analyzed/C01'}; 
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
distmin.help    = {'Defines a minimal geometric distance between a source and a detector for data to be included (in centimeter).'}; 

%Maximum distance
distmax         = cfg_entry; 
distmax.name    = 'Maximum distance'; 
distmax.tag     = 'distmax';      
distmax.strtype = 'r';       
distmax.num     = [1 1];     
distmax.val     = {6};
distmax.help    = {'Defines the maximal geometric distance between a source and a detector for data to be included (in centimeter).'}; 

%Number of MUX 
nb_Mux         = cfg_entry; %nb_Mux
nb_Mux.name    = 'MUX number';
nb_Mux.tag     = 'nb_Mux';       
nb_Mux.strtype = 'r';
nb_Mux.num     = [1 1];     
nb_Mux.val     = {32}; 
nb_Mux.help    = {'Number of MUX of ISS system it changes the illumination order of the source, it is set and test in mux 32'}; 

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
    'in the interpolation when one needs to fill in for missing electrode positions'}'; 

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
        'therefore beware of aliasing artifacts.'}'; 

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


STD_enable          = cfg_menu;
STD_enable.tag      = 'STD_enable';
STD_enable.name     = 'Use the correction';
STD_enable.labels   = {'No','Yes'};
STD_enable.values   = {0,1};
STD_enable.val      =  {0};
STD_enable.help     = {'Remove channels with high STD variation.'}';

STD_amp         = cfg_entry;
STD_amp.tag     = 'STD_amp';
STD_amp.name    = 'STD';
STD_amp.strtype = 'r';       
STD_amp.num     = [1 1];     
STD_amp.val     = {0.1};
STD_amp.help    = {'Standard deviation threshold value.'};


STD_menu        = cfg_menu;
STD_menu.tag    = 'STD_menu';
STD_menu.name   = 'Apply option';
STD_menu.labels = {'Segment','Whole block'};
STD_menu.values = {1,2}; 
STD_menu.val    = {1};
STD_menu.help   = {'Look on different segment std deviation or whole block, segment will keep the channel if only a small parts are lost.'}';

STD_amp_choice         = cfg_branch;
STD_amp_choice.tag     = 'STD_amp_choice';
STD_amp_choice.name    = 'Standard deviation (STD) criteria';
STD_amp_choice.val     = {STD_enable STD_amp,STD_menu};
STD_amp_choice.help    = {'Defines STD minimal to avoid channels measuring an abnormally noisy and high variation signal in the data. High STD channels will be marked as rejected. '};

DC_enable           = cfg_menu;
DC_enable.tag       = 'DC_enable';
DC_enable.name      = 'Use the correction';
DC_enable.labels    = {'No','Yes'};
DC_enable.values    = {0,1};
DC_enable.val       =  {0};
DC_enable.help      = {'Remove channels with low DC intensity.'};

DC_amp          = cfg_entry;
DC_amp.tag      = 'DC_amp';
DC_amp.name     = 'Amplitude (DC';
DC_amp.strtype  = 'r';       
DC_amp.num      = [1 1];     
DC_amp.val      = {100};
DC_amp.help     = {'Remove channel with light intensity lower than a value. DC refer to light intensity level in contrast to AC or Phase measure also available in ISS sytem'}';

DC_amp_choice         = cfg_branch;
DC_amp_choice.tag     = 'DC_amp_choice';
DC_amp_choice.name    = 'Light intensity (DC) criteria';
DC_amp_choice.val     = {DC_enable DC_amp};
DC_amp_choice.help    = {'Defines a minimal light intensity to ensure good signal detection. Low light intensity channels will be marked as rejected.'};



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

 






inputNIRxscout         = cfg_files; %Select raw BOXY data files for this subject 
inputNIRxscout.name    = 'Select NIRx NIRScout data'; % The displayed name
inputNIRxscout.tag     = 'inputNIRxscout';    %file names
inputNIRxscout.ufilter = '.hdr';    %
inputNIRxscout.num     = [1 Inf];     % Number of inputs required 
inputNIRxscout.help    = {'Recording files (.evt, .hdr, .wl1, .wl2, .inf, .tpl, .set) are placed in the same folder during recording sessions. Indicate the .hdr file containing the raw data. Trigger, wavelengths, and parameters will then be read and associated with the helmet project .prj. Ideally you do not rename the files, or if you do, rename them all identically and keep them in the same folder.'}; 


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
e_manualname_NIRxscout.help    = {'Specify data file name, if several files, separate the names by a comma.',...
    'As an example: file01, file02.'}'; 


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

%use10_10system
b_shortdistancenewproject          = cfg_menu;
b_shortdistancenewproject.tag      = 'b_shortdistancenewproject';
b_shortdistancenewproject.name     = 'Create project SD';
b_shortdistancenewproject.labels   = {'True','False'};
b_shortdistancenewproject.values   = {1,0};
b_shortdistancenewproject.val      = {0};
b_shortdistancenewproject.help     = {'NIRx short distance probes record additional detectors next to a source.',...
    'During data recording, these multiple detector probes are encoded in the data of last detector D24.',...  
    'If you set this option to yes, a new project is created to display additional detector ',...
    'close to the associated source where they were place.If you set this option to not, the actual project with detector D24 is kept' }';


b_shortdistanceavailable         = cfg_branch;
b_shortdistanceavailable.tag     = 'b_shortdistanceavailable';
b_shortdistanceavailable.name    = 'Yes';
b_shortdistanceavailable.val     = {e_shortdistancesrs e_shortdistancedet b_shortdistancenewproject };
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
E_readNIRxscout.help = {'Read data acquired with NIRScout on a  NIRx system.'};

function vout = nirs_cfg_vout_readNIRxscout(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


inputNIRxSport        = cfg_files; %Select raw BOXY data files for this subject 
inputNIRxSport.name    = 'Select NIRSport data'; % The displayed name
inputNIRxSport.tag     = 'inputNIRxSport';    %file names
inputNIRxSport.ufilter = '.hdr';    %
inputNIRxSport.num     = [1 Inf];     % Number of inputs required 
inputNIRxSport.help    = {'Recording files (.hdr, .wl1, .wl2, .lsl ) are placed in the same folder during recording sessions. Indicate the .hdr file containing the raw data. Trigger, wavelengths, and parameters will then be read and associated with the helmet project .prj. Ideally you do not rename the files, or if you do, rename them all identically and keep them in the same folder.'}; 



%Maximum distance
Pruning_nirsport         = cfg_entry; 
Pruning_nirsport.name    = 'Default ChannelMask'; 
Pruning_nirsport.tag     = 'Pruning_nirsport';      
Pruning_nirsport.strtype = 's';       
Pruning_nirsport.num     = [0 inf];     
Pruning_nirsport.val     = {' '};
Pruning_nirsport.help    = {['Optional, NIRsport already prune channel according to channel mask and NIRsite montage configuration.']}; 

 
%Maximum distance
distmax_nirsport         = cfg_entry; 
distmax_nirsport.name    = 'Maximum distance'; 
distmax_nirsport.tag     = 'distmax_nirsport';      
distmax_nirsport.strtype = 's';       
distmax_nirsport.num     = [0 inf];     
distmax_nirsport.val     = {'8'};
distmax_nirsport.help    = {['Optional, NIRsport already prune channel according to channel mask and NIRsite montage configuration. If all channel have been kept example you could use the maximal distance to help you prune the channel, else keep this information empty. ' ...
    'Defines the maximal geometric distance between a source and a detector for data to be included (in centimeter).' ...
    'It could also be a fix channel mask. ']}; 

channelmask_nirsport         = cfg_entry; 
channelmask_nirsport.name    = 'Custom ChannelMask'; 
channelmask_nirsport.tag     = 'channelmask_nirsport';      
channelmask_nirsport.strtype = 's';       
channelmask_nirsport.num     = [0 inf];     
channelmask_nirsport.val     = {''};
channelmask_nirsport.help    = {['Custom channel mask channel mask apply to .wl1 and .wl2 use nirs converter and ensure that all raw data are present 16x16 channel',...
    ' 1  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  1  1  1  1  1  0  0  0  1  0  0  0  0  0  0  0  0  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1']}; 




% Create an option on wether a project file should be created or imported.
c_PruningNIRSport          = cfg_choice;
c_PruningNIRSport.tag      = 'c_PruningNIRSport';
c_PruningNIRSport.name     = 'Pruning';
c_PruningNIRSport.values   = {Pruning_nirsport distmax_nirsport channelmask_nirsport};
c_PruningNIRSport.val      = {Pruning_nirsport};
c_PruningNIRSport.help     = {'Choose whether you want to import a project or create a new one.'}';


% Executable Branch
E_readNIRSport     = cfg_exbranch;      
E_readNIRSport.name = 'Read NIRSport';            
E_readNIRSport.tag  = 'E_readNIRSport'; 
E_readNIRSport.val  = {inputNIRxSport,age1,prjfile,c_PruningNIRSport , output_path,c_nameconvention_NIRxscout};   
E_readNIRSport.prog = @nirs_run_readNIRSport;  
E_readNIRSport.vout = @nirs_cfg_vout_readNIRSport;
E_readNIRSport.help = {'Read data acquired with NIRSport on a  NIRx system. Tested for 16x16 montage'};

function vout = nirs_cfg_vout_readNIRSport(job)
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

% Create the two branch required.
b_importProject         = cfg_branch; % Contains the project file selector.
b_importProject.tag     = 'b_importProject';
b_importProject.name    = 'Import .prj';
b_importProject.val     = {prjfile};
b_importProject.help    = {'Import an already existing project file.'};

b_createProject         = cfg_branch; % Is an empty branch.
b_createProject.tag     = 'b_createProject';
b_createProject.name    = 'Create .prj';
b_createProject.val     = {};
b_createProject.help    = {'Create a new project file from the .snirf file.'};

% Create an option on wether a project file should be created or imported.
c_createImportProjectSnirf          = cfg_choice;
c_createImportProjectSnirf.tag      = 'c_createImportProjectSnirf';
c_createImportProjectSnirf.name     = 'Create or import a project file';
c_createImportProjectSnirf.values   = {b_importProject, b_createProject};
c_createImportProjectSnirf.val      = {b_importProject};
c_createImportProjectSnirf.help     = {'Choose whether you want to import a project or create a new one.'}';

% Executable Branch -- Create the drop down menu branch necessary to read .snirf file.
E_readSNIRF      = cfg_exbranch;      
E_readSNIRF.name = 'Read .snirf';            
E_readSNIRF.tag  = 'E_readSNIRF'; 
E_readSNIRF.val  = {inputSNIRF,age1,c_createImportProjectSnirf,output_path};
E_readSNIRF.prog = @nirs_run_readSNIRF; % Handle to the SNIRF import function.
E_readSNIRF.vout = @nirs_cfg_vout_readSNIRF;
E_readSNIRF.help = {'Import the .snirf please install https://github.com/fNIRS/snirf_homer3 to support snirf class'};

function vout = nirs_cfg_vout_readSNIRF(job)
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


% Executable Branch
E_rawhomer      = cfg_exbranch;      
E_rawhomer.name = 'Read .nirs (HoMER)';            
E_rawhomer.tag  = 'E_rawhomer'; 
E_rawhomer.val  = {inputrawhomer,age1, c_createImportProjectSnirf,output_path};   
E_rawhomer.prog = @nirs_run_readhomerfile;  
E_rawhomer.vout = @nirs_cfg_vout_readhomerfile;
E_rawhomer.help = {'Select raw .nirs data files. Matlab structure containing fields: ''d'' (raw data sample x channel), ''ml''(channel used sources, detector, weight wavelength), ''t'' (Time vector), data, ml, S vector will be used as AUX trigger 1.'};

function vout = nirs_cfg_vout_readhomerfile(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


input_GenericDataExportBV         = cfg_files; %Select raw BOXY data files for this subject 
input_GenericDataExportBV.name    = 'Select .dat'; % The displayed name
input_GenericDataExportBV.tag     = 'input_GenericDataExportBV';    %file names
input_GenericDataExportBV.ufilter = '.dat';    %
input_GenericDataExportBV.num     = [1 Inf];     % Number of inputs required 
input_GenericDataExportBV.help    = {'Select Generic Data Export from BrainVision Analyser'}; 


% Executable Branch
E_GenericDataExportBV      = cfg_exbranch;      
E_GenericDataExportBV.name = 'Read .dat (Generic Data export BV)';            
E_GenericDataExportBV.tag  = 'E_GenericDataExportBV'; 
E_GenericDataExportBV.val  = {input_GenericDataExportBV, age1,output_path};   
E_GenericDataExportBV.prog = @nirs_run_readGenericDataExportBV;  
E_GenericDataExportBV.vout = @nirs_cfg_vout_GenericDataExportBV;
E_GenericDataExportBV.help = {'Select raw generic data export from brain vision analyser. Patch pour load EEG file and use target PCA for cardiac correction without ECG. '};

function vout = nirs_cfg_vout_GenericDataExportBV(job)
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
e_MultimodalPath.num      = [0 1];     % Number of inputs required 
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
E_NIRSmatdiradjust.help = {'NIRS.mat structure contains the information about all the localization of the data files on the hard drive.',...
    'This function will modify all the subdirectories of the data files at the actual localization  of the NIRS.mat.',...
    'Be careful, it is not yet tested for paths that the dir contains the same name many times.'};

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
e_NIRSmatdirnewbranch.help    = {'A copy of the last step will be created to create a new branch of analysis.'}; 


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
E_NIRSmatcreatenewbranch.help = {'NIRS.mat structure contains the information about all the localization of the data files on the hard drive.',...
    'This function will create a copy of the NIRS.mat structure to create a new independent branch of analysis.'};

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
E_aux2manualname.help    = {'Indicate a name for the output file (OutputName).'};

E_aux2manualoption         = cfg_menu; 
E_aux2manualoption.name    = 'Option';
E_aux2manualoption.tag     = 'E_aux2manualoption';      
E_aux2manualoption.labels  = {'Separate file'};
E_aux2manualoption.values  = {1};
E_aux2manualoption.val     = {1}; 
E_aux2manualoption.help    = {'Separate file : identify each trig as a separated file.',...
    'Concatenate file : Special use if you plan to concatenate your files in nirs file.',...
    'This option will adjust the trig time assuming that each file from the project is placed,'...
    'one after the other'};

%Create trigger for avg from GLMonset 
E_aux2manualtrig      = cfg_exbranch;
E_aux2manualtrig.name = 'AuxTrig to ManualTrig';
E_aux2manualtrig.tag  = 'E_aux2manualtrig';
E_aux2manualtrig.val  = {NIRSmat, E_aux2manualname};
E_aux2manualtrig.prog = @nirs_run_E_aux2manualtrig;
E_aux2manualtrig.vout = @nirs_cfg_vout_E_aux2manualtrig;
E_aux2manualtrig.help = {'Use ISS AUX trig to create editable manual trig file to be used with manual trig.',...
    'A file will be created in the current directory under the name (fileid_trigid_trig.m).' };

function vout = nirs_cfg_vout_E_aux2manualtrig(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


f_correlationsignal         = cfg_files; %Select raw BOXY data files for this subject 
f_correlationsignal.name    = 'Select .fig signal model to find'; % The displayed name
f_correlationsignal.tag     = 'f_correlationsignal';    %file names
f_correlationsignal.ufilter = '.fig';    %
f_correlationsignal.num     = [1 Inf];     % Number of inputs required 
f_correlationsignal.help    = {'Select .fig with the model of curve to find in the data.'}; 

e_onset        = cfg_entry; %subfolder new branch
e_onset.name    = 'Trigger';
e_onset.tag     = 'e_onset';       
e_onset.strtype = 's';
e_onset.num     = [1 Inf];     
e_onset.val     = {'8'}; 
e_onset.help    = {'Use integer to be mark as a trigger.'}; 

e_correlationsignal_min       = cfg_entry; %subfolder new branch
e_correlationsignal_min.name    = 'Correlation minimal';
e_correlationsignal_min.tag     = 'e_correlationsignal_min';       
e_correlationsignal_min.strtype = 's';
e_correlationsignal_min.num     = [1 Inf];     
e_correlationsignal_min.val     = {'0.6'}; 
e_correlationsignal_min.help    = {'Use minimal correlation value to be mark as a local peak of correlation to mark the event'}; 


e_correlationsignal_LPF       = cfg_entry; %subfolder new branch
e_correlationsignal_LPF.name    = 'Low Pass Filter';
e_correlationsignal_LPF.tag     = 'e_correlationsignal_LPF';       
e_correlationsignal_LPF.strtype = 's';
e_correlationsignal_LPF.num     = [1 Inf];     
e_correlationsignal_LPF.val     = {'No'}; 
e_correlationsignal_LPF.help    = {'Use a numerical value to apply Low pass filter before correlation. The filter will only be apply to mark the event with the correlation'}; 


e_correlationsignal_HPF       = cfg_entry; %subfolder new branch
e_correlationsignal_HPF.name    = 'High Pass Filter';
e_correlationsignal_HPF.tag     = 'e_correlationsignal_HPF';       
e_correlationsignal_HPF.strtype = 's';
e_correlationsignal_HPF.num     = [1 Inf];     
e_correlationsignal_HPF.val     = {'No'}; 
e_correlationsignal_HPF.help    = {'Use a numerical value to apply High pass filter before correlation.'}; 

%General entry use by any modules
zonecorrelation         = cfg_files; 
zonecorrelation.name    = 'zone';
zonecorrelation.tag     = 'zonecorrelation';      
zonecorrelation.filter  = {'ele','zone','txt'};
zonecorrelation.ufilter = '.*';    
zonecorrelation.num     = [1 Inf];     % Number of inputs required 
zonecorrelation.help    = {'Select a .zone create in the GUI to design channel where to apply correlation.',...
    'Or .ele file listing few electrode to find there label. Fp1 0 0 0'};

%Create marker based on correlation on signal based on recording
E_createonset_correlationsignal     = cfg_exbranch;
E_createonset_correlationsignal.name = 'Mark based on correlation';
E_createonset_correlationsignal.tag  = 'E_createonset_correlationsignal';
E_createonset_correlationsignal.val  = {NIRSmat,DelPreviousData, f_correlationsignal, e_onset,e_correlationsignal_min ,e_correlationsignal_LPF  , e_correlationsignal_HPF, zonecorrelation};
E_createonset_correlationsignal.prog = @nirs_run_E_createonset_correlationsignal;
E_createonset_correlationsignal.vout = @nirs_cfg_vout_E_createonset_correlationsignal;
E_createonset_correlationsignal.help = {'Add trig based correlation on feature on the signal' };

function vout = nirs_cfg_vout_E_createonset_correlationsignal(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%Create marker based on correlation on signal based on recording
E_eegnirs_MarkMuscular     = cfg_exbranch;
E_eegnirs_MarkMuscular.name = 'Identify EEG/NIRS muscular artifact ';
E_eegnirs_MarkMuscular.tag  = 'E_eegnirs_MarkMuscular';
E_eegnirs_MarkMuscular.val  = {NIRSmat,DelPreviousData};
E_eegnirs_MarkMuscular.prog = @nirs_run_E_eegnirs_MarkMuscular;
E_eegnirs_MarkMuscular.vout = @nirs_cfg_vout_E_eegnirs_MarkMuscular;
E_eegnirs_MarkMuscular.help = {'Mark as artefact EEG muscular event, they are identify as beta gamma burst higher than zscore distribution at p<0.05 and are also see as fast and local hemodynamic variation in NIRS, espacialy jaw and chewing artefact' };

function vout = nirs_cfg_vout_E_eegnirs_MarkMuscular(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end




%Create marker based on correlation on signal based on recording
E_GoNoGotrig     = cfg_exbranch;
E_GoNoGotrig.name = 'GoNoGotrig';
E_GoNoGotrig.tag  = 'E_GoNoGotrig';
E_GoNoGotrig.val  = {NIRSmat,output_path};
E_GoNoGotrig.prog = @nirs_run_E_GoNoGotrig;
E_GoNoGotrig.vout = @nirs_cfg_vout_E_GoNoGotrig;
E_GoNoGotrig.help = {'Function special for GonoGo malnutrition project , check the good answer in the EEG trig ' };

function vout = nirs_cfg_vout_E_GoNoGotrig(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



e_onset_trig        = cfg_entry; %subfolder new branch
e_onset_trig.name    = 'Trigger';
e_onset_trig.tag     = 'e_onset_trig';       
e_onset_trig.strtype = 's';
e_onset_trig.num     = [1 Inf];     
e_onset_trig.val     = {'8'}; 
e_onset_trig.help    = {'Use integer to be mark as a trigger in the display GUI to indicate the onset event.'}; 

e_CardiacExpectedFrequency         = cfg_entry; %subfolder new branch
e_CardiacExpectedFrequency.name    = 'Expected cardiac beat (Hz)';
e_CardiacExpectedFrequency.tag     = 'e_CardiacExpectedFrequency';       
e_CardiacExpectedFrequency.strtype = 's';
e_CardiacExpectedFrequency.num     = [1 Inf];     
e_CardiacExpectedFrequency.val     = {'2'}; 
e_CardiacExpectedFrequency.help    = {'Use to help to determine the time window to dectect the local maximal peak of the cardiac artifact'}; 

e_Cardiac_HighPass         = cfg_entry; %subfolder new branch
e_Cardiac_HighPass.name    = 'High Pass detection (Hz)';
e_Cardiac_HighPass.tag     = 'e_Cardiac_HighPass';       
e_Cardiac_HighPass.strtype = 's';
e_Cardiac_HighPass.num     = [1 Inf];     
e_Cardiac_HighPass.val     = {'10'}; 
e_Cardiac_HighPass.help    = {'Apply a High pass filter to ease the detection of the high frequency artefact of the cardiac pulse and use the local maximal to dectection on this channel '}; 


%Create marker based on correlation on signal based on recording
E_markCardiac_TargetPCA     = cfg_exbranch;
E_markCardiac_TargetPCA.name = 'Mark cardiac artifact (EEG)';
E_markCardiac_TargetPCA.tag  = 'E_markCardiac_TargetPCA';
E_markCardiac_TargetPCA.val  = {NIRSmat,DelPreviousData,e_CardiacExpectedFrequency, e_Cardiac_HighPass, e_onset_trig};
E_markCardiac_TargetPCA.prog = @nirs_run_E_markCardiac_TargetPCA;
E_markCardiac_TargetPCA.vout = @nirs_cfg_vout_E_markCardiac_TargetPCA;
E_markCardiac_TargetPCA.help = {'Mark cardiac artifac in EEG (without ECG) using target PCA.',...
    'Apply a high pass filter to facilitate the identification of the cardiac artifact.',...
    'use the local maximal peak and select a small time window to apply the target PCA and substract the first component' };


function vout = nirs_cfg_vout_E_markCardiac_TargetPCA(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

e_Cardiac_Window            = cfg_entry; %subfolder new branch
e_Cardiac_Window.name    = 'Target PCA time window (ms)';
e_Cardiac_Window.tag     = 'e_Cardiac_Window';       
e_Cardiac_Window.strtype = 's';
e_Cardiac_Window.num     = [1 Inf];     
e_Cardiac_Window.val     = {'100'}; 
e_Cardiac_Window.help    = {'Time window around the local maximum to apply the target PCA'}; 


%Create marker based on correlation on signal based on recording
E_correctCardiac_TargetPCA     = cfg_exbranch;
E_correctCardiac_TargetPCA.name = 'Correct using Target PCA on marker';
E_correctCardiac_TargetPCA.tag  = 'E_correctCardiac_TargetPCA';
E_correctCardiac_TargetPCA.val  = {NIRSmat,DelPreviousData,e_Cardiac_Window,e_onset_trig};
E_correctCardiac_TargetPCA.prog = @nirs_run_E_correctCardiac_TargetPCA;
E_correctCardiac_TargetPCA.vout = @nirs_cfg_vout_E_correctCardiac_TargetPCA;
E_correctCardiac_TargetPCA.help = {'Correct cardiac artifac in EEG (without ECG) using target PCA.',...
    'Apply a high pass filter to facilitate the identification of the cardiac artifact.',...
    'use the local maximal peak and select a small time window to apply the target PCA and substract the first component' };


function vout = nirs_cfg_vout_E_correctCardiac_TargetPCA(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



E_correctCardiac_exportBVFolder        = cfg_files;
E_correctCardiac_exportBVFolder.name    = 'Select output folder'; 
E_correctCardiac_exportBVFolder.tag     = 'E_correctCardiac_exportBVFolder';       %file names
E_correctCardiac_exportBVFolder.filter   = {'dir'};
E_correctCardiac_exportBVFolder.ufilter = '.*';    
E_correctCardiac_exportBVFolder.num     = [1 Inf];       % Number of inputs required 
E_correctCardiac_exportBVFolder.help    = {'Select the output folder where to export data.'}; 



%Create marker based on correlation on signal based on recording
E_correctCardiac_exportBV     = cfg_exbranch;
E_correctCardiac_exportBV.name = 'Save generic data export BV';
E_correctCardiac_exportBV.tag  = 'E_correctCardiac_exportBV';
E_correctCardiac_exportBV.val  = {NIRSmat,DelPreviousData,E_correctCardiac_exportBVFolder};
E_correctCardiac_exportBV.prog = @nirs_run_E_correctCardiac_exportBV;
E_correctCardiac_exportBV.vout = @nirs_cfg_vout_E_correctCardiac_exportBV;
E_correctCardiac_exportBV.help = {'Save Data corrected to be reimported in BrainVision ' };


function vout = nirs_cfg_vout_E_correctCardiac_exportBV(job)
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
m_trigmode.help     = {'Replace trig will erase the default trigger present during the recording to place only the trig you have defined manually.',...
                   'Add trig will add new trig to the previous trig already present in the file definition.'};

m_trigfile         = cfg_files;
m_trigfile.name    = 'Enter new trig definition'; 
m_trigfile.tag     = 'm_trigfile';       %file names
m_trigfile.filter  = 'm';
m_trigfile.ufilter = 'trig.m$';    
m_trigfile.num     = [1 Inf];     % Number of inputs required 
m_trigfile.help    = {'Open ManualTrig file.',...
    'The manual trig is a .m Matlab script that defines the new TrigValue=X an integer used for the trigger identification and the timingfile{1}=[15,30]; which indicates the timing of the trig event in seconds for the file 1.',...
    'As the NIRS.mat structure could contain one or many files you must define a timingfile for each of the files present in your data.',...
    'If, as an example, the file 2 does not have any trig events defined, set the timingfile as empty value: timingfile{2} = [ ];.  Be aware that .m file name does not support space or accent in the name definition. ',...
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




m_trigfilelsl         = cfg_files;
m_trigfilelsl.name    = 'Enter new trig definition'; 
m_trigfilelsl.tag     = 'm_trigfilelsl';       %file names
m_trigfilelsl.filter  = 'lsl.tri';
m_trigfilelsl.ufilter = 'lsl.tri';    
m_trigfilelsl.num     = [1 Inf];     % Number of inputs required 
m_trigfilelsl.help    = {'Open file'}; 


E_manualtriglsl      = cfg_exbranch;
E_manualtriglsl.name = 'Trigger lsl.tri';
E_manualtriglsl.tag  = 'E_manualtriglsl';
E_manualtriglsl.val  = {NIRSmat m_trigfilelsl m_trigmode};
E_manualtriglsl.prog = @nirs_run_E_manualtriglsl;
E_manualtriglsl.vout = @nirs_cfg_vout_E_manualtriglsl;
E_manualtriglsl.help = {['Read lsl.tri file to add trigger as example the file contain one line by trigger as follow ',...
    'Timestamp;sample in nirs ; triggervalue',...
    '2024-07-26T10:35:02.901947;1481;0',...
    '2024-07-26T10:35:40.891001;1867;3',...
    '2024-07-26T10:36:43.853145;2508;2',...
    '2024-07-26T10:36:56.863759;2640;1']};

function vout = nirs_cfg_vout_E_manualtriglsl(job)
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
e_Concatenate_blocid.help       = {'Enter the file number to group, as an example 1,2,3. Set 0 to take them all.',...
    'In the display, a file list will be identified such as Bloc001.'};


m_Concatenate_option        = cfg_menu;
m_Concatenate_option.tag    = 'm_Concatenate_option';
m_Concatenate_option.name   = 'Options';
m_Concatenate_option.labels = {'Merge only','Merge and Detrend'};
m_Concatenate_option.values = {0,1};
m_Concatenate_option.val    = {0};
m_Concatenate_option.help   = {'Merge data off all nirs and block selected.',...
    'Merge only does not correct any offset while merge and detrend does a detrending operation on each block before to merge them.'};


f_Concatenate_outdir        = cfg_files;
f_Concatenate_outdir.name    = 'Select output folder'; 
f_Concatenate_outdir.tag     = 'f_writeNIRSdir';       %file names
f_Concatenate_outdir.filter  = {'dir'};
f_Concatenate_outdir.ufilter = '.*';    
f_Concatenate_outdir.num     = [0 1];      % Number of inputs required 
f_Concatenate_outdir.val    = {''};
f_Concatenate_outdir.help    = {'Select the output folder where to save data. By default if you let empty it will save in the current NIRS.mat file'}; 



E_Concatenate_file      = cfg_exbranch;
E_Concatenate_file.name = 'Concatenate several data files in the same NIRS.mat';
E_Concatenate_file.tag  = 'E_Concatenate_file';
E_Concatenate_file.val  = {NIRSmat,m_Concatenate_option,e_Concatenate_blocid, f_Concatenate_outdir  };
E_Concatenate_file.prog = @nirs_run_E_Concatenate_file;
E_Concatenate_file.vout = @nirs_cfg_vout_E_Concatenate_file;
E_Concatenate_file.help = {'Join several data files in the NIRS.mat, not support EEG, AUX or video multimodal information. Use concatenate NIRS.mat to combine multiple NIRS.mat file '};

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
f_nirsmatinfo.help    = {'Create the list of subjects to group in a new folder example: ./GrandAverage/List.xlsx.',...
    'The excel file must have 2 columns: first, the NIRS.mat subject location and second, the channel list order to respect or zone to use. The result will keep the same helmet of the first subject in the list.'}; 

m_Concatenate_Exclude           = cfg_menu;
m_Concatenate_Exclude.tag       = 'm_Concatenate_Exclude';
m_Concatenate_Exclude.name      = 'Set exclude channel to NaN';
m_Concatenate_Exclude.labels    = {'Exclude','Keep'};
m_Concatenate_Exclude.values    = {0,1};
m_Concatenate_Exclude.val       = {0};
m_Concatenate_Exclude.help      = {'Set excluded channels to NaN.',...
    'Exclude: rejected channel will be set to NaN. ',...
    'Keep: no special process concerning the rejected channels will be applied.'};

m_Concatenate_Normalized        = cfg_menu; 
m_Concatenate_Normalized.tag    = 'm_Concatenate_Normalized';
m_Concatenate_Normalized.name   = 'Normalization';
m_Concatenate_Normalized.labels = {'No normalization','Min-Max Normalization','z-score Normalization'};
m_Concatenate_Normalized.values = {0,1,2};
m_Concatenate_Normalized.val    = {0};
m_Concatenate_Normalized.help   = {'Apply a normalization to reduce individual subject  variability. Only apply using the channel list.',...
    'No normalization: does not apply any normalization. ',...
    'Min-Max Normalization: applies min max normalization fixing boundary between 0 and 1.',... 
    'Z-score Normalization: applies z-score normalization.',... 
    'Moeller, J., 2015. A word on standardization in longitudinal studies: don�t. Front Psychol 6. https://doi.org/10.3389/fpsyg.2015.01389'};

E_Concatenate_nirsmat      = cfg_exbranch;
E_Concatenate_nirsmat.name = 'Concatenate NIRS.mat (multi subject and zone)';
E_Concatenate_nirsmat.tag  = 'E_Concatenate_nirsmat';
E_Concatenate_nirsmat.val  = {f_nirsmatinfo,m_Concatenate_option, m_Concatenate_Exclude, m_Concatenate_Normalized};
E_Concatenate_nirsmat.prog = @nirs_run_E_Concatenate_nirsmat;
E_Concatenate_nirsmat.vout = @nirs_cfg_vout_E_Concatenate_nirsmat;
E_Concatenate_nirsmat.help = {'Join several NIRS.mat one after the other. Could be used to do group average.',...
    'Individual multimodal (EEG, AUX, Video, AUDIO) are not supported for this function.'};

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
timedurationinternan.help    = {'Normalization will be done between each artifact period'};

b_choiceinternan         = cfg_branch;
b_choiceinternan.tag     = 'b_choiceinternan';
b_choiceinternan.name    = 'Normalization inter-artifacts';
b_choiceinternan.val     = {timedurationinternan};
b_choiceinternan.help    = {'Normalize each segment by its mean, segments are separated by artifacted periods.'}';

trigger         = cfg_entry;
trigger.name    = 'Trigger';
trigger.tag     = 'trigger';
trigger.strtype = 'r';
trigger.num     = [0 Inf];
trigger.val     = {1}; 
trigger.help    = {'Enter trigger number. If you used Comments to select an interval to segment let the this field empty'};

pretime         = cfg_entry;
pretime.name    = 'PreTime';
pretime.tag     = 'pretime';       
pretime.strtype = 's';
pretime.num     = [0 Inf];
pretime.val     = {'5'}; 
pretime.help    = {'Time to include before the trigger.',...
    'Write keyword ''start'' to go to the beginning of the raw segment.',...
    'Define as a positive value, as example 5 seconds before the onset. Unintuitively, -30 will give you 30 seconds after the trigger.',...
    'To write 0 use ''0''',...
    'Write the keyword ''Comments evt01Start to stop the segment at the comment named evt01Start',...
    'Write the same keyword for pretime and postime if you like to event start and be segment for the whole duration of the event.'};

posttime         = cfg_entry;
posttime.name    = 'PostTime';
posttime.tag     = 'posttime';       
posttime.strtype = 's';
posttime.num     = [0 Inf];
posttime.val     = {'30'};
posttime.help    = {'Time to include after the trigger.',...
    'Write keyword ''end'' to go to the end of the raw segment.',...
    'Write keyword ''trig9'' to go until this next specific trigger. Test only for one trig case yet',...
    'Write the keyword ''Comments evt01End to stop the segment at the comment named evt01End',...
    'Write the same keyword as pretime if you like to use the duration of the event marked in the comments'};


m_NormType          = cfg_menu;
m_NormType.tag      = 'm_NormType';
m_NormType.name     = 'How define Io';
m_NormType.labels   = {'Io = Pretime to 0','Io = Pretime to PostTime'};%,'Substract I-Io (Pretime','Do nothing, Raw segmented'};
m_NormType.values   = {0,1};%,2,3
m_NormType.val      = {0};
m_NormType.help     = {'Choose one of the two definitions of Io'}';



b_choicenormstim         = cfg_branch;
b_choicenormstim.tag     = 'b_choicenormstim';
b_choicenormstim.name    = 'Normalization around trigger';
b_choicenormstim.val     = {trigger pretime posttime m_NormType m_choiceNan};
b_choicenormstim.help    = {'Normalization using pre and post time triggers.'}';

c_normtype          = cfg_choice;
c_normtype.tag      = 'normtype';
c_normtype.name     = 'Normalization type';
c_normtype.values   = {b_choiceglobal,b_choicenormstim,b_choiceinternan};
c_normtype.val      = {b_choiceglobal};
c_normtype.help     = {''}';

% Executable Branch
E_normalization      = cfg_exbranch;
E_normalization.name = 'Normalization dOD = log(I/Io)';
E_normalization.tag  = 'normalization';
E_normalization.val  = {NIRSmat DelPreviousData c_normtype};
E_normalization.prog = @nirs_run_normalize;
E_normalization.vout = @nirs_cfg_vout_normalize;
E_normalization.help = {'Normalize raw data, transfer optical intensity in delta optical density;',...
    'dOD  = log(I/Io).'};
%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_normalize(~)
vout = cfg_dep;                    
vout.sname      = 'NIRS.mat';       
vout.src_output = substruct('.','NIRSmat'); 
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end


m_SegmentTrig          = cfg_menu;
m_SegmentTrig.tag      = 'm_SegmentTrig';
m_SegmentTrig.name     = 'Option trig';
m_SegmentTrig.labels   = {'Use all','Only first'};
m_SegmentTrig.values   = {0,1};
m_SegmentTrig.val      = {0};
m_SegmentTrig.help     = {'By defaults, use all trig to segment. However, in some cases use only first trig to segment the whole bloc of the data'};


% Executable Branch
segment      = cfg_exbranch;
segment.name = 'Segment';
segment.tag  = 'segment';
segment.val  = {NIRSmat DelPreviousData trigger pretime posttime m_SegmentTrig };
segment.prog = @nirs_run_segment;
segment.vout = @nirs_cfg_vout_segment;
segment.help = {'This module segments data around triggers (defining pretime and post time). This step is essential to synchronize multimodal data such as EEG, auxiliary (AUX)  or video. '};
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
PrintReport.help    = {'Save an image of the percentage of rejected channels over time in the NIRS.mat folder. Each criterion will save a summary report if they are applied.'}';

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
thresholdstep.help    = {'The difference (D) between M1 and M2 over time is converted into a standardized z-score distribution. A z-score threshold determines which variations in the moving average difference of two subsequent intervals are considered abnormal, thus represent artifacts. Depending on the value defined by the user,  variations higher than this z-score identified.'};
   
m_thresholdstepzscore        = cfg_menu;
m_thresholdstepzscore.tag      = 'm_thresholdstepzscore';
m_thresholdstepzscore.name     = 'Z-score option';
m_thresholdstepzscore.labels   = {'All-time points','Valid time points','Old version - using median value threshold'};
m_thresholdstepzscore.values   = {1,0,2};
m_thresholdstepzscore.val      = {1}; 
m_thresholdstepzscore.help     = {'The z-score is used to define the threshold for artifact detection. You may calculate the z-score for each channel using all-time points or using only valid time points (excluding ones identified as an artifact).'};


too_small_step_dur         = cfg_entry;
too_small_step_dur.name    = 'Ignored interval shorter than x second';
too_small_step_dur.tag     = 'too_small_step_dur';
too_small_step_dur.strtype = 'r';
too_small_step_dur.num     = [1 1];
too_small_step_dur.val     = {0.4}; 
too_small_step_dur.help    = {'This specifies the minimal duration (seconds) of an abnormal variation in the signal, as defined by the z-score threshold of the moving average, in order to be considered as an artifact. Intervals where the z-score threshold appears for a period shorter than this time are ignored.'};

movingaverage_nbpoint         = cfg_entry;
movingaverage_nbpoint.name    = 'Moving average step (seconds)';
movingaverage_nbpoint.tag     = 'movingaverage_nbpoint';
movingaverage_nbpoint.strtype = 'r';
movingaverage_nbpoint.num     = [1 1];
movingaverage_nbpoint.val     =  {1}; 
movingaverage_nbpoint.help    = {'Defines the duration (sec) of the time interval (step) to calculate moving average M1 and M2: Xn to Xn+step. Depending on the smoothing the user wants to apply, i.e. how sensitive detection of the signal�s variation should be, a higher or lower value should be indicated.'};


 
printreportthreshold        = cfg_menu;
printreportthreshold.tag    = 'printreportthreshold';
printreportthreshold.name   = 'Print threshold report for each channel';
printreportthreshold.labels = {'Yes','No'};
printreportthreshold.values = {1,0};
printreportthreshold.val    = {0}; 
printreportthreshold.help   = {'Select if you want to save a figure that displays the threshold detection made for each channel. The figure will be saved in �\ArtifactDetection_Report\EachCh folder at the NIRS.mat home location.'}';

b_meandiff         = cfg_branch;
b_meandiff.tag     = 'b_meandiff';
b_meandiff.name    = 'Artifact detection using moving average';
b_meandiff.val     = {m_meandiff thresholdstep m_thresholdstepzscore movingaverage_nbpoint too_small_step_dur printreportthreshold};
b_meandiff.help    = {'Moving average allows identifying discontinuity or strong perturbations in the signal. It is a criterion based on the signal mean variation for two subsequent time intervals: M1 from Xn to Xn+step and M2 from Xn+1 to Xn+step+1. The step is defined by the user and defines the period of time for the moving average time windows. Changes of light intensity between M1 and M2 are then specified as the difference (D) between both means: D = M2-M1. Difference D is transferred into z-scores normalized over the entire dataset to determine which periods have increased signal variations. The threshold, i.e. the z-score above which abnormal variations in the signal are identified as artifact intervals, is to be defined by the user. '}';

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
corr_thr.help    = {'Set the minimal correlation threshold between the channels. If the correlation is stronger then this threshold, the time course will be marked as artifacts.'};

b_corr         = cfg_branch;
b_corr.tag     = 'b_corr';
b_corr.name    = 'Correlation between channels for artifact interval';
b_corr.val     = {m_corr corr_thr};
b_corr.help    = {'Artifact detection is based on a correlation coefficient between channels. It determines the Pearson�s correlation coefficient threshold of channels to be considered as being affected by the same event. For each artifact interval that has previously been detected, channels that have a correlation equal or above the threshold with the time course of this artifact, are detected as well. This means that it asks for each artifact to find the channels showing the same signature, but have not been  detected based on the previous criteria.'}';

m_minpourcentagebad      = cfg_menu;
m_minpourcentagebad.tag  = 'm_minpourcentagebad';
m_minpourcentagebad.name = 'Apply';
m_minpourcentagebad.labels = {'Yes','No'};
m_minpourcentagebad.values = {1,0};
m_minpourcentagebad.val  = {1}; 
m_minpourcentagebad.help = {'Use Minimal percentage of bad channels to be marked as artifacts.'}';

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
b_minpourcentagebad.help    = {'In some case, the signal is clean and very few channels are detected as artifacts.',...
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
    'These four criteria are: �Artifact detection using moving average�; �Minimal percentage of bad channels to be marked as artifact�; �Minimal subinterval�; �Correlation between channels for artifact�.',...
    'A more detailed description of each criterion is presented above.'};




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


m_cardiacwavelenght           = cfg_menu; 
m_cardiacwavelenght.tag       = 'm_cardiacwavelenght';    
m_cardiacwavelenght.name      = 'Wavength to apply rejection'; 
m_cardiacwavelenght.labels    = {'First','Second', 'Both','Any'};
m_cardiacwavelenght.values    = {1,2,3,4};  
m_cardiacwavelenght.val       = {[4]};
m_cardiacwavelenght.help      = {'Define the wavelenght were rejected channel will be selected. The cardiac coherence is perform on both optical wavelenght',...
    'You could reject channel base on the first wavelenght only, the second only, both of them need to be good or any of them need to be good to keep the channel as valid'}; 

% Executable Branch
E_chcardiaccontrol      = cfg_exbranch;
E_chcardiaccontrol.name = 'Cardiac Detection';
E_chcardiaccontrol.tag  = 'E_chcardiaccontrol';
E_chcardiaccontrol.val  = {NIRSmat, i_Freq_cardiac i_COHTRESHOLD_cardiac i_minch_cardiac i_cardiacwidth m_cardiacwavelenght };
E_chcardiaccontrol.prog = @nirs_run_chcardiaccontrol;
E_chcardiaccontrol.vout = @nirs_run_vout_chcardiaccontrol;
E_chcardiaccontrol.help = {'Detection of the cardiac beat using coherence measure among channels and the module rejects channels without any cardiac evidence.'};




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
E_nullifybad.help = {'This module replaces all artifacts identified in the data (yellow parts previously identified by Artifact detection modules) by a missing value (NaN).'};

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
    'Enter a value if you want to apply a low-pass filter, if you don''t, write ''No''.'};

highcutfreq         = cfg_entry;
highcutfreq.name    = 'High-Pass cutoff';
highcutfreq.tag     = 'highcutfreq';       
highcutfreq.strtype = 's';
highcutfreq.num     = [1 inf];
highcutfreq.val     = {'No'};
highcutfreq.help    = {'Cutoff for High-Pass filter.',...
    'Enter a value if you want to apply a high-pass filter, if you don''t, write ''No''.'};

paddingsymfilter         = cfg_menu;
paddingsymfilter.tag     = 'paddingsymfilter';
paddingsymfilter.name    = 'Symmetric Padding';
paddingsymfilter.labels  = {'True','False' }';                
paddingsymfilter.values  = {1 0};
paddingsymfilter.val     = {1};
paddingsymfilter.help    = {'To avoid border effect on short segment filtering, apply symmetrization of the segment before filtering.'};


interpolatebadfilter         = cfg_menu;
interpolatebadfilter.tag     = 'interpolatebadfilter';
interpolatebadfilter.name    = 'Interpolate bad interval';
interpolatebadfilter.labels  = {'True','False' }';                
interpolatebadfilter.values  = {1 0};
interpolatebadfilter.val     = {0};
interpolatebadfilter.help    = {'Bad interval or period in yellow are interpolated before filtering if ''yes'' is selected.'};


filterorder         = cfg_entry;
filterorder.name    = 'Filter order';
filterorder.tag     = 'filterorder';       
filterorder.strtype = 'r';
filterorder.num     = [1 Inf];
filterorder.val     = {4};
filterorder.help    = {'Butterworth filter order. Design a Butterworth digital filter and apply using a zero-phase digital filtering (filtfilt.m). '};

% Executable Branch
E_filter      = cfg_exbranch;
E_filter.name = 'Filter';
E_filter.tag  = 'bpfilt';
E_filter.val  = {NIRSmat DelPreviousData lowcutfreq highcutfreq filterorder paddingsymfilter interpolatebadfilter};
E_filter.prog = @nirs_run_filter;
E_filter.vout = @nirs_cfg_vout_filter;
E_filter.help = {'Bandpass filtering of the data. Bad intervals are interpolated prior to filtering.',...
    'Input: Normalized intensity or concentrations. Design a Butterworth digital filter and apply using a zero-phase digital filtering (filtfilt.m).'};

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
E_detrend.help = {'Use detrend on each block.',...
    'Subtract linear trends on the data. The trend is estimated using the beginning and end of the block.'};

function vout = nirs_cfg_vout_detrend(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end




E_prewhitening = cfg_exbranch;
E_prewhitening.name = 'PreWhitening';
E_prewhitening.tag  = 'E_prewhitening';
E_prewhitening.val  = {NIRSmat DelPreviousData }; 
E_prewhitening.prog = @nirs_run_E_prewhitening;
E_prewhitening.vout = @nirs_cfg_vout_prewhitening;
E_prewhitening.help = {'Use Prewhitening see ref: , TOBE TESTED RAE'};



function vout = nirs_cfg_vout_prewhitening(job);
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
b_ODtoHbOHbR_DPF1.help    = {'Adjust differential pathlength factor depending on the wavelength and age of the subject.',...
    'DOI: 10.1117/1.JBO.18.10.105004'};


%Choice
C_ODtoHbOHbR_DPF         = cfg_choice;
C_ODtoHbOHbR_DPF.tag     = 'C_ODtoHbOHbR_DPF';
C_ODtoHbOHbR_DPF.name    = 'DPF method';
C_ODtoHbOHbR_DPF.values  = {b_ODtoHbOHbR_DPF1,b_ODtoHbOHbR_DPF2,b_ODtoHbOHbR_DPF3 };
C_ODtoHbOHbR_DPF.val     = {b_ODtoHbOHbR_DPF1};
C_ODtoHbOHbR_DPF.help    = {'Adjust differential pathlenght factor (DPF) in Modified Beer Lambert Law calculation.'};


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
            'Use the wavelength extenction coefficient.',...
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
i_norm_dur.help    = {'Set amplitude variation threshold between 2 good subinterval group the interval or separate it'};




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
avg_datatype.help    = {'Boxy data give DC, AC and PH component. Select the one you want to average.'};


badintervalratio         = cfg_entry;
badintervalratio.tag     = 'badintervalratio';
badintervalratio.name    = 'Reject trial ratio';
badintervalratio.strtype = 'r';
badintervalratio.num     = [1 Inf];
badintervalratio.val     = {0.5};
badintervalratio.help    = {'Reject the trial if more than xx% of its duration is marked as a bad interval, set to 1 to keep trial.',...
    'Enter the percentage as a decimal number. Example: 0.5 for 50%'};

badchannelratio         = cfg_entry;
badchannelratio.tag     = 'badchannelratio';
badchannelratio.name    = 'Reject channel ratio';
badchannelratio.strtype = 'r';
badchannelratio.num     = [1 inf];
badchannelratio.val     = {0.5};
badchannelratio.help    = {'It Rejects the channel if less than xx% of the trial is rejected, set to 0 to keep all channels.',...
    'Enter the percentage as a decimal number. Example: 0.5 for 50%'};

helpmemoryprob          = cfg_menu;
helpmemoryprob.tag      = 'helpmemoryprob';
helpmemoryprob.name     = 'Do you need to solve help memory problem';
helpmemoryprob.labels   = {'No','No'};
helpmemoryprob.values   = {0,0};
helpmemoryprob.val      = {0};
helpmemoryprob.help     = {'Not available now, Help to solve help memory problem especially useful for large .nir data file, each segment trial will be read independently in spite of loading the whole file.'};

e_base_pretime         = cfg_entry;
e_base_pretime.name    = 'PreTime baseline correction';
e_base_pretime.tag     = 'e_base_pretime';       
e_base_pretime.strtype = 'r';
e_base_pretime.num     = [1 1];
e_base_pretime.val     = {0}; 
e_base_pretime.help    = {'Pretime for the baseline correction (Trigger reference = time 0).'}; 

e_base_posttime         = cfg_entry;
e_base_posttime.name    = 'PostTime baseline correction';
e_base_posttime.tag     = 'e_base_posttime';       
e_base_posttime.strtype = 'r';
e_base_posttime.num     = [1 1];
e_base_posttime.val     = {0}; 
e_base_posttime.help    = {'PostTime for the baseline correction (Trigger reference = time 0).'}; 

m_mean          = cfg_menu;
m_mean.tag      = 'm_mean';
m_mean.name     = 'Subtract value';
m_mean.labels   = {'Mean','Median'};
m_mean.values   = {0,1};
m_mean.val      = {0};
m_mean.help     = {'Subtract value using mean or median.'};

b_manualbaseline_corr        = cfg_branch;
b_manualbaseline_corr.tag    = 'b_manualbaseline_corr';
b_manualbaseline_corr.name   = 'Manual';
b_manualbaseline_corr.val    = {m_mean,e_base_pretime,e_base_posttime};
b_manualbaseline_corr.help   = {'Manual time baseline correction for each stim.'}; 

m_nobaseline_corr           = cfg_menu;
m_nobaseline_corr.tag       = 'm_nobaseline_corr';
m_nobaseline_corr.name      = 'No baseline correction';
m_nobaseline_corr.labels    = {'No'};
m_nobaseline_corr.values    = {0};
m_nobaseline_corr.val       = {0};
m_nobaseline_corr.help      = {'Keep the original data and average without applying baseline correction.'};

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
m_Tvalueoption.help   = {'Select one of the options:',...
    'Simple t-test against 0', 'Simple t-test against mean baseline value'}';

m_noreject_trial      = cfg_menu;
m_noreject_trial.tag  = 'm_noreject_trial';
m_noreject_trial.name = 'Keep trial';
m_noreject_trial.labels = {'Yes'};
m_noreject_trial.values = {1};
m_noreject_trial.val  = {1};
m_noreject_trial.help = {'Keep all trials, except if the trial is noisier than the reject trial ratio.'}; 


e_reject_outlier_threshold         = cfg_entry;
e_reject_outlier_threshold.name    = 'Threshold z-score';
e_reject_outlier_threshold.tag     = 'e_reject_outlier_threshold';       
e_reject_outlier_threshold.strtype = 'r';
e_reject_outlier_threshold.num     = [1 1];
e_reject_outlier_threshold.val     = {6};
e_reject_outlier_threshold.help    = {'Exclude trial outliers based on the z-score compute in comparison to the other trials. The z-score is computed individually for each channel by using the mean and standard deviation of all trials.'};

m_reject_outlier_printreport        = cfg_menu;
m_reject_outlier_printreport.tag    = 'm_reject_outlier_printreport';
m_reject_outlier_printreport.name   = 'Print report';
m_reject_outlier_printreport.labels = {'No','Yes'};
m_reject_outlier_printreport.values = {0,1};
m_reject_outlier_printreport.val    = {1};
m_reject_outlier_printreport.help   = {'Print a report of rejected trial.'}; 


b_reject_trial        = cfg_branch;
b_reject_trial.tag    = 'b_reject_trial';
b_reject_trial.name   = 'Reject outlier trial z-score';
b_reject_trial.val    = {e_reject_outlier_threshold,m_reject_outlier_printreport};
b_reject_trial.help   = {'Warning: Bad trials will be set as NaN, but no traces of which trials are rejected are written in the vmrk file. Further versions may include track of excluded channel in DisplayGui.'};


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
choiceave.val     = {trigger pretime posttime badintervalratio badchannelratio,avg_datatype,c_baseline_corr,c_rejecttrial,m_Tvalueoption};
choiceave.help    = {''}'; %avtype 

% Executable Branch
E_average      = cfg_exbranch;
E_average.name = 'Epoch averaging';
E_average.tag  = 'E_average';
E_average.val  = {NIRSmat DelPreviousData choiceave}; % savenirs
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
% PCA.help = {'� ENLEVER Apply PCA on selected time'};
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
NIRSname.help    = {'OutputName only used to name the concatenate file, else initial name will be kept.'};


NIRS_export_Concatenatefile         = cfg_branch;
NIRS_export_Concatenatefile.tag     = 'NIRS_export_Concatenatefile';
NIRS_export_Concatenatefile.name    = 'Concatenate file';
NIRS_export_Concatenatefile.val     = {NIRSname};
NIRS_export_Concatenatefile.help    = {'Export .nirs files with all segments one after the other.'}';

NIRS_export_Separatefile         = cfg_branch;
NIRS_export_Separatefile.tag     = 'NIRS_export_Separatefile';
NIRS_export_Separatefile.name    = 'Separate file';
NIRS_export_Separatefile.val     = {};
NIRS_export_Separatefile.help    = {'Export each file in a separated .nirs file. The output file is located in the same folder of the NIRS.mat file.'}';

%to be delete
NIRS_exportoption         = cfg_menu;
NIRS_exportoption.name    = 'Option';
NIRS_exportoption.tag     = 'NIRS_exportoption';      
NIRS_exportoption.labels  = {'Separate file','Concatenate file' };
NIRS_exportoption.values  = {2,1};
NIRS_exportoption.val     = {2};
NIRS_exportoption.help    = {'Export each file in a separated .nirs file or in one nirs file with all segments one after the other. '};


c_NIRS_exportoption         = cfg_choice;
c_NIRS_exportoption.tag     = 'c_NIRS_exportoption';
c_NIRS_exportoption.name    = 'Options';
c_NIRS_exportoption.values  = {NIRS_export_Separatefile,NIRS_export_Concatenatefile };
c_NIRS_exportoption.val     = {NIRS_export_Separatefile}; %Default option
c_NIRS_exportoption.help    = {''};

%Old version to be remove
% % Executable Branch
% E_writeNIRS      = cfg_exbranch;
% E_writeNIRS.name = 'Write NIRS';
% E_writeNIRS.tag  = 'WriteNIRS';
% E_writeNIRS.val  = {NIRSsession,NIRSmat,c_NIRS_exportoption};
% E_writeNIRS.prog = @nirs_run_writenirs;
% E_writeNIRS.vout = @nirs_cfg_vout_writenirs;
% E_writeNIRS.help = {'Write in .nirs format for Homer.'};
% 
% function vout = nirs_cfg_vout_writenirs(job)
%     vout = cfg_dep;                    
%     vout.sname      = 'NIRS.mat';       
%     vout.src_output = substruct('.','NIRSmat'); 
%     vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

f_writeNIRSdir        = cfg_files;
f_writeNIRSdir.name    = 'Select output folder'; 
f_writeNIRSdir.tag     = 'f_writeNIRSdir';       %file names
f_writeNIRSdir.filter  = {'dir'};
f_writeNIRSdir.ufilter = '.*';    
f_writeNIRSdir.num     = [1 1];      % Number of inputs required 
f_writeNIRSdir.help    = {'Select the output folder where to export data. Data will be placed in a subfolder with the subject name to facilitate the organization.'}; 

% Executable Branch
E_writeNIRSHomer      = cfg_exbranch;
E_writeNIRSHomer.name = 'Write .nirs homer';
E_writeNIRSHomer.tag  = 'E_writeNIRSHomer';
E_writeNIRSHomer.val  = {NIRSmat,f_writeNIRSdir};
E_writeNIRSHomer.prog = @nirs_run_writenirshomer;
E_writeNIRSHomer.vout = @nirs_cfg_vout_writenirshomer;
E_writeNIRSHomer.help = {'Write in .nirs format for Homer.'};

function vout = nirs_cfg_vout_writenirshomer(job)
    vout = cfg_dep;                    
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
SNIRFname.help    = {'Used to name the resulting output .snirf file.'};

SNIRF_outpath         = cfg_entry; %path
SNIRF_outpath.name    = 'Folder path for output .snirf file';
SNIRF_outpath.tag     = 'SNIRF_outpath';       
SNIRF_outpath.strtype = 's';
SNIRF_outpath.num     = [1 Inf];     
SNIRF_outpath.val     = {'C:\data\'}; 
SNIRF_outpath.help    = {'The path of the folder in which the .snirf file is going to be created.',...
    'Ex: C:\data\user1\SNIRFfolder\'}; 

% f_SNIRFfile         = cfg_files;
% f_SNIRFfile.name    = 'Select existing SNIRF file'; 
% f_SNIRFfile.tag     = 'f_SNIRFfile';       %file names
% f_SNIRFfile.filter  = 'snirf';
% f_SNIRFfile.ufilter = '.snirf$';    
% f_SNIRFfile.num     = [1 1];     % Number of inputs required 
% f_SNIRFfile.help    = {'Select the SNIRF file to which you want to append your data.'}; 
% 
% B_SNIRFAppend       = cfg_branch;
% B_SNIRFAppend.name  = 'Append to an existing SNIRF file';
% B_SNIRFAppend.tag   = 'B_SNIRFAppend';
% B_SNIRFAppend.val   = {f_SNIRFfile};
% B_SNIRFAppend.help   = {'To append the data to an existing file.'};
% 
% B_SNIRFCreate       = cfg_branch;
% B_SNIRFCreate.name  = 'Create a new SNIRF file';
% B_SNIRFCreate.tag   = 'B_SNIRFCreate';
% B_SNIRFCreate.val   = {SNIRFname SNIRF_outpath};
% B_SNIRFCreate.help  = {'Export the data to a new SNIRF file.'};
% 
% B_SNIRFExportLocation         = cfg_choice;
% B_SNIRFExportLocation.tag     = 'B_SNIRFExportLocation';
% B_SNIRFExportLocation.name    = 'Choose SNIRF export type';
% B_SNIRFExportLocation.val     = {B_SNIRFCreate};
% B_SNIRFExportLocation.values  = {B_SNIRFCreate  B_SNIRFAppend};
% B_SNIRFExportLocation.help    = {'Choose whether you want to export your data to a new SNIRF file or append it to an existing SNIRF file.'};

% Executable Branch
E_writeSNIRF      = cfg_exbranch;
E_writeSNIRF.name = 'Write .snirf';
E_writeSNIRF.tag  = 'E_writeSNIRF';
E_writeSNIRF.val  = {NIRSmat,f_writeNIRSdir};
E_writeSNIRF.prog = @nirs_run_writeSNIRF;
E_writeSNIRF.vout = @nirs_cfg_vout_writeSNIRF;
E_writeSNIRF.help = {'This toolbox enables you to export your data as .snirf files. A description of this format is available here https://github.com/fNIRS/snirf/blob/master/snirf_specification.md. LIONirs uses the class develop in Homer3 to save .snirf format please refer and install  https://github.com/fNIRS/snirf_homer3. A similar data structure to .nirs are used to save the .snirf format the export will contain : data: ''d''; probe coordinate: ''SD''; sample time: ''t''; stim trigger: ''s''; auxiliary export are not well-supported: ''aux''.'};

function vout = nirs_cfg_vout_writeSNIRF(job)
    vout = cfg_dep;                    
    vout.sname      = 'SNIRF.mat';       
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
    'It will cause an error if there is any epoch averaging in the NIRS.mat session.' };

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
NIR_FileIn.help    = {'Select a .nir file to extract a smaller time segment and create a new .nir file.'}; 


NIR_START_TIME         = cfg_entry;
NIR_START_TIME.name    = 'Time start';
NIR_START_TIME.tag     = 'NIR_START_TIME';       
NIR_START_TIME.strtype = 'r';
NIR_START_TIME.num     = [1 1];
NIR_START_TIME.val     = {0};
NIR_START_TIME.help    = {'The starting time of the selection.'};

NIR_STOP_TIME         = cfg_entry;
NIR_STOP_TIME.name    = 'Time stop';
NIR_STOP_TIME.tag     = 'NIR_STOP_TIME';       
NIR_STOP_TIME.strtype = 'r';
NIR_STOP_TIME.num     = [1 1];
NIR_STOP_TIME.val     = {0};
NIR_STOP_TIME.help    = {'The ending time of the selection.'};


% Executable Branch
E_NIR_segment      = cfg_exbranch;
E_NIR_segment.name = 'Write NIR indivual file segment';
E_NIR_segment.tag  = 'Write_NIR_segment';
E_NIR_segment.val  = {NIRSmat, NIR_FileIn, NIR_START_TIME,NIR_STOP_TIME};
E_NIR_segment.prog = @nirs_run_E_NIR_segment;
E_NIR_segment.vout = @nirs_run_vout_E_NIR_segment;
E_NIR_segment.help = {'Use a .nir file to write a smaller segment, can be used to specify only the baseline for the PCA artifact, for example.'};

function vout = nirs_run_vout_E_NIR_segment(job)
    vout            = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
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
HMR_datad1mat.help        = {'Select a matrix of data to display in the same order as the subject in the hmr file. Let empty if you want to select concentration data type'}; 

HMR_savenormfig          = cfg_menu;
HMR_savenormfig.tag      = 'HMR_savenormfig';
HMR_savenormfig.name     = 'Save Normalised Figure';  
HMR_savenormfig.labels   = {'No','Yes'}; 
HMR_savenormfig.values   = {0,1};
HMR_savenormfig.val      = {1};
HMR_savenormfig.help     = {'Save concentration and normalized concentration figures'};

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

NIRS_dataHemo           = cfg_menu;
NIRS_dataHemo.tag       = 'NIRS_dataHemo';
NIRS_dataHemo.name      = 'Concentration';
NIRS_dataHemo.labels    = { 'HbO and HbR','HbO', 'HbR',};
NIRS_dataHemo.values    = {1, 2, 3};
NIRS_dataHemo.val       = {1};
NIRS_dataHemo.help      = {'Display HbO and HbR or wavelength first and second if concentration were not compute'};


b_NIRSdata        = cfg_branch;
b_NIRSdata.tag    = 'b_NIRSdata';
b_NIRSdata.name   = 'NIRS epoch average data';
b_NIRSdata.val    = {NIRSmat NIRS_datatype HMR_savenormfig, NIRS_dataHemo};
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
vColordata.help    = {'Select a serie of topographic files as .img or structure TOPO.mat. ',...
    'If you have selected .img file note that they should have been done on the same skin or cortex reconstruction.',...
    'Only one .prj can be used then. '}; 

b_vColordata        = cfg_branch;
b_vColordata.tag    = 'b_vColordata';
b_vColordata.name   = 'Load topographic map';
b_vColordata.val    = {vColordata};
b_vcolordata.help   = {''};

typedata         =  cfg_choice;
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
m_prjfile_mode.help     = {'The project file (.prj) used to create the topography can be paired which all file selected above, file1 use prj1, file2 use prj2...',...
    'You can also use the same prj for all the topography',...
    'In case of the selection Load topographic map is selected, the same prj skin and cortex have been used for all files'};

%prjfile define above 
m_projection_mode          = cfg_menu;
m_projection_mode.tag      = 'm_projection_mode';
m_projection_mode.name     = 'Projection';  
m_projection_mode.labels   = {'Banana shape, radial projection', 'Inverse weighted distance'}; 
m_projection_mode.values   = {0,1};
m_projection_mode.val      = {0};
m_projection_mode.help     = {'Topographic projection',...
    'Banana shape, project color scale amplitude using the radial angle in spherical coordinate between source and detector.',...
    'Inverse weighted distance use the intersection between source and detector and project color scale using inverse weighted distance.'};


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
v_view.help     = {['Write view angle list as : ''head rotation angle,inclination angle;...''  as an example ''90,30;-90,30'' ',...
    'will produce a left view with 30 degrees of inclination and a right view with 30 degrees of inclination. ',...
    'Rotation angle : 0 front, 90 left, 180, back, -90 right view ',...
    'inclination angle : 0 no inclination, 90 top view ',....
    'If you enter nothing no figure or video will be produced. ']};

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
v_outpath.help    = {['Path for video output files: write video_cond1(omit backslashes) data will be saved in the subject folder like  ...subject\dataspm\video\ ',...
    'or write complete directory name as c:/data/video_cond1 to define manual the path']}; 

v_outtimeprecisionleft         = cfg_entry;
v_outtimeprecisionleft.name    = 'Number of left decimal output name time';
v_outtimeprecisionleft.tag     = 'v_outtimeprecisionleft';       
v_outtimeprecisionleft.strtype = 'r';
v_outtimeprecisionleft.num     = [1 Inf];     
v_outtimeprecisionleft.val     = {4}; 
v_outtimeprecisionleft.help    = {['Output file name will be written with this number of decimal example : 4=>0012.XX(s), 3 => 012.XXX(s)']}; 

v_outtimeprecision         = cfg_entry;
v_outtimeprecision.name    = 'Number of decimal output name time';
v_outtimeprecision.tag     = 'v_outtimeprecision';       
v_outtimeprecision.strtype = 'r';
v_outtimeprecision.num     = [1 Inf];     
v_outtimeprecision.val     = {2}; 
v_outtimeprecision.help    = {['Output file name will be written with this number of decimal example : 2=>0XX.23(s), 0 => 0XX(s)']}; 

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
v_showcover.val    = {0};
v_showcover.help   = {'Show uncovered region by the montage in grey'}';

b_option         = cfg_branch;
b_option.tag     = 'b_option';
b_option.name    = 'Define display option';
b_option.help    = {'Display option for HSJ helmet.'};
b_option.val     = {prjfile m_prjfile_mode m_projection_mode, v_cmin,v_cmax,v_ctr,v_view,v_start,v_stop,v_step,v_outpath,v_outtimeprecisionleft,v_outtimeprecision,v_showcover,v_skin};

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
m_extractcomponentfigure.val    = {0}; 
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
i_extractnoise_nbPCA.help    = {'Find number of components to identify as an artifact in PCA setting the minimal xx % of explained variance.'};

b_extractcomponent_PCA        = cfg_branch;
b_extractcomponent_PCA.tag    = 'b_extractcomponent_PCA';
b_extractcomponent_PCA.name   = 'Identify PCA for artifact period';
b_extractcomponent_PCA.val    = {NIRSmat,  m_extractcomponent,m_extractcomponentfigure,i_extract_pourcentagech, i_extractnoiseupto_nbPCA, i_extractnoise_nbPCA};
b_extractcomponent_PCA.help   = {'First, identify noisy intervals using artifact detection or a manual revision. This function runs a PCA decomposition (targetPCA) on each bad interval (yellow segment in the DisplayGUI). The decomposition is performed on identified channels during a continuous bad interval. PCA decomposition sort components according to the explained variance. The component(s) explaining the highest variance during the artifactual interval is assumed to be mainly related to the artifact event. They will be stored in the component list with the label MVTPCA. We recommend using the module �Subtract Components� that will subtract all the components identified with a specific label.'};

f_extractcomponent_physzone         = cfg_files;
f_extractcomponent_physzone.name    = 'Enter Regressor zone'; 
f_extractcomponent_physzone.tag     = 'f_extractcomponent_physzone';       %file names
f_extractcomponent_physzone.filter  = 'zone';
f_extractcomponent_physzone.ufilter = '.zone';    
f_extractcomponent_physzone.num     = [1 Inf];     % Number of inputs required regression 
f_extractcomponent_physzone.help    = {'Use to define channel use in the short distance regression: regressor zone1 will regress in zone1.'}; 

m_extractcomponent_physzone        = cfg_menu;
m_extractcomponent_physzone.tag    = 'm_extractcomponent_physzone';
m_extractcomponent_physzone.name   = 'Use ';
m_extractcomponent_physzone.labels = {'Mean' ,'PCA'};
m_extractcomponent_physzone.values = {0,1};
m_extractcomponent_physzone.val    = {0}; 
m_extractcomponent_physzone.help   = {'Use the mean or the first principal component of the regressor zone; useful when the regressor zone contains many channels.'};

b_extractcomponent_phys      = cfg_branch;
b_extractcomponent_phys.tag  = 'b_extractcomponent_phys';
b_extractcomponent_phys.name = 'Identify physiology regression';
b_extractcomponent_phys.val  = {NIRSmat,  f_extractcomponent_physzone ,m_extractcomponent_physzone};
b_extractcomponent_phys.help = {'Apply subtraction of the short distance physiology as describe in Saager and Berger 2008, https://doi.org/10.1117/1.2940587.',...
                                    'One signal was scaled to fit the other in a least-squares LS sense the scale estimation is saved as (SHORTGLM) component and must be subtracted.',...
                                    'When artefact period (yellow) are marked, the regression will be apply piecewise for each segment of data conserved.'};

m_replaceglmlist       = cfg_menu;
m_replaceglmlist.tag    = 'm_replaceglmlist';
m_replaceglmlist.name   = 'Replace current component ';
m_replaceglmlist.labels = {'Yes' ,'No'};
m_replaceglmlist.values = {1,2};
m_replaceglmlist.val    = {1}; 
m_replaceglmlist.help   = {'Avoid to run twice the extract GLM on the data, juste erase the present component if they have some in the actual data'};
                               
m_glmlist_autoexport_HbO       = cfg_menu;
m_glmlist_autoexport_HbO.tag    = 'm_glmlist_autoexport_HbO';
m_glmlist_autoexport_HbO.name   = 'Export chromophore';
m_glmlist_autoexport_HbO.labels = {'HbO' ,'HbR','HbO&HbR'};
m_glmlist_autoexport_HbO.values = {1,2,3};
m_glmlist_autoexport_HbO.val    = {1}; 
m_glmlist_autoexport_HbO.help   = {''};

i_glmlist_autoexport_Xi         = cfg_entry;
i_glmlist_autoexport_Xi.name    = 'Select which regressor to export X0=0, X1=1, X2=2...';
i_glmlist_autoexport_Xi.tag     = 'i_glmlist_autoexport_Xi';       
i_glmlist_autoexport_Xi.strtype = 'r';
i_glmlist_autoexport_Xi.num     = [1 Inf];
i_glmlist_autoexport_Xi.val     = {0}; 
i_glmlist_autoexport_Xi.help    = {'Associate to the regressor X0, X1, X2,... to export'};

i_glmlist_autoexport_labelbad_threshold         = cfg_entry;
i_glmlist_autoexport_labelbad_threshold.name    = 'Percentage of yellow interval to identify trial as rejected';
i_glmlist_autoexport_labelbad_threshold.tag     = 'i_glmlist_autoexport_labelbad_threshold';       
i_glmlist_autoexport_labelbad_threshold.strtype = 'r';
i_glmlist_autoexport_labelbad_threshold.num     = [1 Inf];
i_glmlist_autoexport_labelbad_threshold.val     = {50}; 
i_glmlist_autoexport_labelbad_threshold.help    = {'Associate to the export topo the suffix ''bad'' if  the associate interval is mark of more than 50% of Yellow',...
                                                    'If you wish also to force them as a missing value use Reject to Rejected trial instead of keep'};

i_glmlist_autoexport_nan_chrejected         = cfg_menu;
i_glmlist_autoexport_nan_chrejected.name    = 'Rejected channel set Beta as NAN';
i_glmlist_autoexport_nan_chrejected.tag     = 'i_glmlist_autoexport_nan_chrejected';       
i_glmlist_autoexport_nan_chrejected.labels  = {'Reject' ,'Keep'};
i_glmlist_autoexport_nan_chrejected.values  = {1,2};
i_glmlist_autoexport_nan_chrejected.val     = {1}; 
i_glmlist_autoexport_nan_chrejected.help    = {'Set as missing value (NAN) the ''Beta'' for rejected channel in dotted for the export '};

                        
i_glmlist_autoexport_nan_int_rejected         = cfg_menu;
i_glmlist_autoexport_nan_int_rejected.name    = 'Rejected trial set Beta as NAN';
i_glmlist_autoexport_nan_int_rejected.tag     = 'i_glmlist_autoexport_nan_int_rejected';       
i_glmlist_autoexport_nan_int_rejected.labels  = {'Reject' ,'Keep'};
i_glmlist_autoexport_nan_int_rejected.values  = {1,2};
i_glmlist_autoexport_nan_int_rejected.val     = {2}; 
i_glmlist_autoexport_nan_int_rejected.help    = {'Set as missing value (NAN) the ''Beta'' for rejected interval in yellow,'...
    'use the Percentage of yellow interval to identify bad trial',... 
    'use keep if they were corrected '};

    % Folder selector.
f_extractglmlist_autoexport_yes_ChannelList          = cfg_files; %path
f_extractglmlist_autoexport_yes_ChannelList.name     = 'Enter Channel List'; 
f_extractglmlist_autoexport_yes_ChannelList.tag      = 'f_extractglmlist_autoexport_yes_ChannelList';
f_extractglmlist_autoexport_yes_ChannelList.filter   = {'txt'};
f_extractglmlist_autoexport_yes_ChannelList.ufilter  = '.txt';    %
f_extractglmlist_autoexport_yes_ChannelList.num      = [0 1];     % Number of inputs required 
f_extractglmlist_autoexport_yes_ChannelList.help     = {'Channel list is important to ensure channel order is repected if different project configuration are used (source and detector position). If you let this field empty a default channel list for this data file will be create'};    

f_extractglmlist_autoexport_yes          = cfg_files; %path
f_extractglmlist_autoexport_yes.name     = 'Result folder';
f_extractglmlist_autoexport_yes.tag      = 'f_extractglmlist_autoexport_yes';
f_extractglmlist_autoexport_yes.filter   = {'dir'};
f_extractglmlist_autoexport_yes.ufilter  = '.*';    %
f_extractglmlist_autoexport_yes.num      = [1 1];     % Number of inputs required 
f_extractglmlist_autoexport_yes.help     = {'Result of component will be saved in this folder.'};
   
b_extractglmlist_autoexport_yes         = cfg_branch;
b_extractglmlist_autoexport_yes.tag      = 'b_extractglmlist_autoexport_yes';
b_extractglmlist_autoexport_yes.name     = 'Yes';
b_extractglmlist_autoexport_yes.val      = {f_extractglmlist_autoexport_yes, f_extractglmlist_autoexport_yes_ChannelList, m_glmlist_autoexport_HbO, i_glmlist_autoexport_Xi, i_glmlist_autoexport_labelbad_threshold,i_glmlist_autoexport_nan_chrejected,i_glmlist_autoexport_nan_int_rejected,m_replaceglmlist};
b_extractglmlist_autoexport_yes.help     = {'Result of component will be saved in this folder.'};

e_extractglmlist_autoexport_no        = cfg_entry;
e_extractglmlist_autoexport_no.name    = 'No';
e_extractglmlist_autoexport_no.tag     = 'c_extractglmlist_autoexport_no';       
e_extractglmlist_autoexport_no.strtype = 's';
e_extractglmlist_autoexport_no.num     = [1 Inf];
e_extractglmlist_autoexport_no.val     = {'No'}; 
e_extractglmlist_autoexport_no.help    = {'Do not save extract'};
                 
                                
                                
c_extractglmlist_autoexport         = cfg_choice; 
c_extractglmlist_autoexport.tag     = 'c_extractglmlist_autoexport';
c_extractglmlist_autoexport.name    = 'Export in folder';
c_extractglmlist_autoexport.values  = {e_extractglmlist_autoexport_no,b_extractglmlist_autoexport_yes};
c_extractglmlist_autoexport.val     = {e_extractglmlist_autoexport_no}; %Default option
c_extractglmlist_autoexport.help    = {'Export components directly in a specify folder.'};

                                
f_extractcomponent_glmlist         = cfg_files;
f_extractcomponent_glmlist.name    = 'List GLM to identify (xls)'; 
f_extractcomponent_glmlist.tag     = 'f_extractcomponent_glmlist';       %file names
f_extractcomponent_glmlist.filter  = {'xlsx','xls','txt'};
f_extractcomponent_glmlist.ufilter = '.*';
f_extractcomponent_glmlist.num     = [0 Inf];     % Number of inputs required 
f_extractcomponent_glmlist.help    = {'Enter the list xls to extract GLM, Let the field empty if you want to used the default configuration create by CreateAUX using trigger',...
    'the list must include the following columns:',...
    '''NIRS.mat folder'': directory of the NIRS.mat to use as: C:\data\Analyze\C01;',...
    '''file'': number to identify the file to use as: 1;',...
    '''tStart'': starting time in seconds to define the beginning of the period where the GLM will be applied as; 10',...
    '''tStop'': time stop in seconds to define the end of the period where the GLM will be applied as; 40',...
    '''Label'': Label identification to write in the event (useful to manage the export)as HRF;',...
    'Regressor(s) (''X0'', ''X1'',''X2'',''X3'',''X4''...): the regressor could be a name in the aux list, a channel zone (regressor zone 1 apply to zone 1) be realy picky on the spelling as; HRF task ',...
    'If you wish to use the one create in the ''Create AUX'' previous operation call ExtractHRF.xlsx in the current directory let the option empty'}; 




b_extractcomponent_glm          = cfg_branch;
b_extractcomponent_glm.tag      = 'b_extractcomponent_glm';
b_extractcomponent_glm.name     = 'Identify GLM';
b_extractcomponent_glm.val      = {NIRSmat, f_extractcomponent_glmlist, c_extractglmlist_autoexport};
b_extractcomponent_glm.help     = {'Apply multiple linear regressions (regress.m).'};


f_extractcomponent_PARAFAClist         = cfg_files;
f_extractcomponent_PARAFAClist.name    = 'List Parafac to identify (xls)'; 
f_extractcomponent_PARAFAClist.tag     = 'f_component_PARAFAClist';       %file names
f_extractcomponent_PARAFAClist.filter  = {'xlsx','xls','txt'};
f_extractcomponent_PARAFAClist.ufilter = '.*';    
f_extractcomponent_PARAFAClist.num     = [1 Inf];     % Number of inputs required 
f_extractcomponent_PARAFAClist.help    = {'Enter the list xls to extract PARAFAC, the list must include the following columns:',...
    '''NIRS.mat folder'': directory of the NIRS.mat to use;',...
    '''file'': number to identify the file to use;',...
    '''tStart'': time start in seconds to define the beginning of the period where the GLM will be applies;',...
    '''tStop'': time stop in seconds to define the end of the period where the GLM will be applied;',...
    '''Label'': Label identification to write in the event (useful manage the export)'};

b_extractcomponent_PARAFAC          = cfg_branch;
b_extractcomponent_PARAFAC.tag      = 'b_extractcomponent_PARAFAC';
b_extractcomponent_PARAFAC.name     = 'Identify PARAFAC';
b_extractcomponent_PARAFAC.val      = {f_extractcomponent_PARAFAClist};
b_extractcomponent_PARAFAC.help     = {'Display option for HSJ helmet.'};



f_extractAVGlist_autoexport_yes          = cfg_files; %path
f_extractAVGlist_autoexport_yes.name     = 'Result folder';
f_extractAVGlist_autoexport_yes.tag      = 'f_extractAVGlist_autoexport_yes';
f_extractAVGlist_autoexport_yes.filter   = {'dir'};
f_extractAVGlist_autoexport_yes.ufilter  = '.*';    %
f_extractAVGlist_autoexport_yes.num      = [1 1];     % Number of inputs required 
f_extractAVGlist_autoexport_yes.help     = {'Result of component will be saved in this folder.'};
   
b_extractAVGlist_autoexport_yes         = cfg_branch;
b_extractAVGlist_autoexport_yes.tag      = 'b_extractAVGlist_autoexport_yes';
b_extractAVGlist_autoexport_yes.name     = 'Yes';
b_extractAVGlist_autoexport_yes.val      = {f_extractAVGlist_autoexport_yes};
b_extractAVGlist_autoexport_yes.help     = {'Result of component will be saved in this folder.'};

e_extractAVGlist_autoexport_no        = cfg_entry;
e_extractAVGlist_autoexport_no.name    = 'No';
e_extractAVGlist_autoexport_no.tag     = 'c_extractAVGlist_autoexport_no';       
e_extractAVGlist_autoexport_no.strtype = 's';
e_extractAVGlist_autoexport_no.num     = [1 Inf];
e_extractAVGlist_autoexport_no.val     = {'No'}; 
e_extractAVGlist_autoexport_no.help    = {'Do not save extract'};

c_extractAVGlist_autoexport         = cfg_choice; 
c_extractAVGlist_autoexport.tag     = 'c_extractAVGlist_autoexport';
c_extractAVGlist_autoexport.name    = 'Export in folder';
c_extractAVGlist_autoexport.values  = {e_extractAVGlist_autoexport_no,b_extractAVGlist_autoexport_yes};
c_extractAVGlist_autoexport.val     = {e_extractAVGlist_autoexport_no}; %Default option
c_extractAVGlist_autoexport.help    = {'Export components directly in a specify folder.'};

f_extractcomponent_AVGlist         = cfg_files;
f_extractcomponent_AVGlist.name    = 'List AVG to identify (xls)'; 
f_extractcomponent_AVGlist.tag     = 'f_component_AVGlist';       %file names
f_extractcomponent_AVGlist.filter  = {'xlsx','xls','txt'};
f_extractcomponent_AVGlist.ufilter = '.*';    
f_extractcomponent_AVGlist.num     = [1 Inf];     % Number of inputs required 
f_extractcomponent_AVGlist.help    = {'NIRS.mat folder: Directory to locate the data to extract.',... 
    '''File'': block in the NIRS.mat data file;',... 
    '''tStart'': to get the curve start point;',... 
    '''tStop: to get the curve stop point;',... 
    '''tStartavg'': to get the average starting point;',... 
    '''tStopavg'': to get the average stopping point;',... 
    '''Label'': Write label in the component name;',... 
    '''ZoneDisplay'': use the first zone channel to plot the average. Or let it empty Keep the zone file in the same folder as the excel ExtractAVG setting.'}; 

b_extractcomponent_AVG          = cfg_branch;
b_extractcomponent_AVG.tag      = 'b_extractcomponent_AVG';
b_extractcomponent_AVG.name     = 'Identify AVG';
b_extractcomponent_AVG.val      = {f_extractcomponent_AVGlist, c_extractAVGlist_autoexport};
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
i_extractnoise_nbPARAFAC.help    = {'Find the optimal number of components to decompose PARAFAC. Try up to x components to find best Concordia and minimise the error .'};

b_extractnoise_PARAFAC          = cfg_branch;
b_extractnoise_PARAFAC.tag      = 'b_extractnoise_PARAFAC';
b_extractnoise_PARAFAC.name     = 'Identify PARAFAC for artifact period';
b_extractnoise_PARAFAC.val      = {NIRSmat,  i_extractnoise_labelPARAFAC,m_extractcomponentfigure,i_extract_pourcentagech,i_extractnoise_nbPARAFAC};
b_extractnoise_PARAFAC.help     = {'First, identify noisy intervals. This function runs a PARAFAC decomposition on each pre-identified bad intervals (yellow segment in the DisplayGUI). The decomposition is performed on the noisy channel during this interval as a Target PARAFAC. The decomposition explaining the most of variance in the data will save as the one representing the artifact  will be stored in the component list with the label MVTPARAFAC. The component list is stored in �SelectedFactors.mat� in the same folder as the NIRS.mat. You could visualize them in the DisplayGUI before subtracting them. We recommend using the module �Subtract Components� that will subtract all the components identified with the same label.'};

c_extractcomponent         = cfg_choice;
c_extractcomponent.tag     = 'c_extractcomponent';
c_extractcomponent.name    = 'Option';
c_extractcomponent.val     = {b_extractcomponent_PCA};
c_extractcomponent.values  = {b_extractcomponent_PCA,b_extractnoise_PARAFAC, b_extractcomponent_phys, b_extractcomponent_glm,b_extractcomponent_PARAFAC,b_extractcomponent_AVG };
c_extractcomponent.help    = {''};


% Executable Branch
E_extractcomponent      = cfg_exbranch;
E_extractcomponent.name = 'Extract';
E_extractcomponent.tag  = 'E_extractcomponent';
E_extractcomponent.val  = {c_extractcomponent};
E_extractcomponent.prog = @nirs_run_E_extractcomponent;
E_extractcomponent.vout = @nirs_cfg_vout_E_extractcomponent;
E_extractcomponent.help = {'Identify data components using several data decomposition methods on target data. The component will be added to a list that could be visualized in the DisplayGUI (extract), subtract in case of an artifactual component using �Subtract component� or export as a relevant activity for further statistics using �Export component�. '};
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
m_substract_and_offsetcorrection.help   = {'Force offset adjustment after a substracted artifact.'};

i_substractcomponent_label         = cfg_entry;
i_substractcomponent_label.name    = 'Component identification';
i_substractcomponent_label.tag     = 'i_substractcomponent_label';       
i_substractcomponent_label.strtype = 's';
i_substractcomponent_label.num     = [1 inf];
i_substractcomponent_label.val     = {'MVTPARAFAC'}; 
i_substractcomponent_label.help    = {'Identify the label of the component you want to subtract, as an example MVTPARAFAC MVTPCA for the automatic detection.'};

% Executable Branch
E_substractcomponent      = cfg_exbranch;
E_substractcomponent.name = 'Subtract';
E_substractcomponent.tag  = 'E_substractcomponent';
E_substractcomponent.val  = {NIRSmat,i_substractcomponent_label,m_substract_and_offsetcorrection};
E_substractcomponent.prog = @nirs_run_E_substractcomponent;
E_substractcomponent.vout = @nirs_cfg_vout_E_substractcomponent;
E_substractcomponent.help = {'Subtract from the data the component identify by the label.'};

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
E_exportcomponent.help = {'Subtract from the data the component identify by the label'};

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
f_component_list.name    = 'Enter list of components to export (xls)'; 
f_component_list.tag     = 'f_component_list';       %file names
f_component_list.filter  = {'xlsx','xls','txt'};
f_component_list.ufilter = '.*';    
f_component_list.num     = [1 Inf];     % Number of inputs required 
f_component_list.help    = {'The export will be organized using the channel list order.  Each subject must have a close localization on the head to be compared. Enter the list to export columns such as:  �NIRS.mat folder�, �Type� (GLM, PARAFAC,PCA), �Label�, Component name to filter, �Channel List� full file or .txt same folder as xls file and �Name� use as the output name of the export.'};

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
f_component_zone.name    = 'Enter list of components to export and there zone (xls)'; 
f_component_zone.tag     = 'f_component_zone';       %file names
f_component_zone.filter  = {'xlsx','xls','txt'};
f_component_zone.ufilter = '.*';    
f_component_zone.num     = [1 Inf];     % Number of inputs required 
f_component_zone.help    = {'Enter the list to export columns such as: �NIRS.mat folder�, �Type� (GLM, PARAFAC,PCA),�Label� ,Component name to filter, �Zone List� full file or .txt same folder as xls file and �Name� the output name of the export .'}; 

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

m_statcomponent_saveoption        = cfg_menu;
m_statcomponent_saveoption.name   = 'Save info.mat data and parameters';
m_statcomponent_saveoption.tag    = 'm_statcomponent_saveoption';       
m_statcomponent_saveoption.labels = {'Yes','No'};
m_statcomponent_saveoption.values = {0,1};
m_statcomponent_saveoption.val    = {0};
m_statcomponent_saveoption.help   = {'Save in the result file a file with the data and parameters used to perform the test. Output file info.mat. These data could be used in external statistical software.'};

e_statcomponent_alpha      = cfg_entry;
e_statcomponent_alpha.name    = 'Alpha threshold p value < ';
e_statcomponent_alpha.tag     = 'e_statcomponent_alpha';       
e_statcomponent_alpha.strtype = 'r';
e_statcomponent_alpha.num     = [1 inf];
e_statcomponent_alpha.val     = {0.05}; 
e_statcomponent_alpha.help    = {'Define the alpha threshold of significance to mask result map.'};

% Folder selector.
e_STATCOMPPath          = cfg_files; %path
e_STATCOMPPath.name     = 'Result folder';
e_STATCOMPPath.tag      = 'e_STATCOMPPath';
e_STATCOMPPath.filter   = {'dir'};
e_STATCOMPPath.ufilter  = '.*';    %
e_STATCOMPPath.num      = [1 1];     % Number of inputs required 
e_STATCOMPPath.help     = {'Result of the statistic will be saved in this folder.'};



f_componentG1.name    = 'Group 1 (components export)'; 
f_componentG1.tag     = 'f_componentG1';       %file names
f_componentG1.filter  = 'mat';
f_componentG1.ufilter = '.mat';    
f_componentG1.num     = [1 Inf];     % Number of inputs required 
f_componentG1.help    = {'Enter the list of components to test statistically.'}; 



f_MCT_ZoneBased        = cfg_files;
f_MCT_ZoneBased.tag    = 'b_MCT_ZoneBased';
f_MCT_ZoneBased.name   = 'Zone' ;
f_MCT_ZoneBased.filter  = 'zone';
f_MCT_ZoneBased.ufilter = '.zone';    
f_MCT_ZoneBased.num     = [1 Inf];     % Number of inputs required 
f_MCT_ZoneBased.help   = {'The cluster need the distance between channel to define neighbord. ',...
    'The zone is use to define channel position and to define a subset of channel if necessary.',...
    'Zone must be define in the display GUI or a zone with all channel is automaticly save at the data import step Global.zone'};

%THIS USED in component stat

e_neighbourdist         =  cfg_entry; %component stat
e_neighbourdist.name    = 'Neighbor distance';
e_neighbourdist.tag     = 'e_neighbourdist';       
e_neighbourdist.strtype = 's';
e_neighbourdist.num     = [0 inf];
e_neighbourdist.val     = {'4'};
e_neighbourdist.help    = { 'Write a number to define the distance between optode to be include in a neighboring to define a cluster. The unit is the same as the channel position. '};

b_MCT_ClusterBased         = cfg_branch;
b_MCT_ClusterBased.tag    = 'b_MCT_ClusterBased';
b_MCT_ClusterBased.name   = 'Cluster Based' ;
b_MCT_ClusterBased.val    = {e_neighbourdist, f_MCT_ZoneBased};
b_MCT_ClusterBased.help   = {'Cluster based permutation'};

b_MCT_Max         = cfg_branch;
b_MCT_Max.tag    = 'b_MCT_Max';
b_MCT_Max.name   = 'No' ;
b_MCT_Max.val    = {};
b_MCT_Max.help   = {'Use the maximal univariate comparaison to create null distribution; Troendle, J.F., 1995. A Stepwise Resampling Method of Multiple Hypothesis Testing. Journal of the American Statistical Association 90, 370�378. https://doi.org/10.1080/01621459.1995.10476522'};

 

c_statMultipleComparaisonTesting        = cfg_choice;
c_statMultipleComparaisonTesting.tag     = 'c_statMultipleComparaisonTesting';
c_statMultipleComparaisonTesting.name    = 'Multiple Comparaison';
c_statMultipleComparaisonTesting.values  = {b_MCT_Max, b_MCT_ClusterBased};
c_statMultipleComparaisonTesting.val     = {b_MCT_ClusterBased}; %Default option
c_statMultipleComparaisonTesting.help    = {'Cluster-based permutation adapt clusterstat function from fieldtrip function and see for reference Maris, E., Oostenveld, R., 2007. Nonparametric statistical testing of EEG- and MEG-data. Journal of Neuroscience Methods 164, 177�190. https://doi.org/10.1016/j.jneumeth.2007.03.024'};


e_npermutation         =  cfg_entry;
e_npermutation.name    = 'Nb permutations';
e_npermutation.tag     = 'e_npermutation';       
e_npermutation.strtype = 's';
e_npermutation.num     = [0 inf];
e_npermutation.val     = {'500'};
e_npermutation.help    = {'Define the number of non-parametric permutations to use.'};


b_permutation        = cfg_branch;
b_permutation.tag    = 'b_permutation';
b_permutation.name   = 'Yes' ;
b_permutation.val    = {e_npermutation, c_statMultipleComparaisonTesting};
b_permutation.help   = {'Apply permutation test'};


b_Nopermutation        = cfg_branch;
b_Nopermutation.tag    = 'b_Nopermutation';
b_Nopermutation.name   = 'No' ;
b_Nopermutation.val    = {};
b_Nopermutation.help   = {'Do not apply permutation test'};

c_statcomppermutation        = cfg_choice;
c_statcomppermutation.tag     = 'c_statpermutation';
c_statcomppermutation.name    = 'Univariate permutations';
c_statcomppermutation.values  = {b_permutation,b_Nopermutation };
c_statcomppermutation.val     = {b_Nopermutation}; %Default option
c_statcomppermutation.help    = {'Use permutation to create the empirical null distribution.'};

m_TtestOneSample        = cfg_menu;
m_TtestOneSample.tag    = 'm_TtestOneSample';
m_TtestOneSample.name   = 'Use ';
m_TtestOneSample.labels = {'2-tailed','1-tailed left (neg) ' ,'1-tailed right (pos) ' };
m_TtestOneSample.values = {1,2,3};
m_TtestOneSample.val    = {1}; 
m_TtestOneSample.help   = {'Use 2-tailed or one-tailed assumption.'};


f_component         = cfg_files;
f_component.name    = 'Enter list of components'; 
f_component.tag     = 'f_component';       %file names
f_component.filter  = 'mat';
f_component.ufilter = '.mat';    
f_component.num     = [1 Inf];     % Number of inputs required 
f_component.help    = {'Enter the list of components to test statistically.'}; 


b_TtestOneSample        = cfg_branch;
b_TtestOneSample.tag    = 'b_TtestOneSample';
b_TtestOneSample.name   = 'One-sample t-test' ;
b_TtestOneSample.val    = {f_component,m_TtestOneSample};
b_TtestOneSample.help   = {'The one-sample t-test is a parametric test of the location parameter when the population standard deviation is unknown.',...
    'T = mean(X)/(STD(x)*sqrt(n))'};


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
    'D1matrix or mean by subjects from the export list.'}; 


b_TtestUnpaired        = cfg_branch;
b_TtestUnpaired.tag    = 'b_TtestUnpaired';
b_TtestUnpaired.name   = 'Unpaired t-test' ;
b_TtestUnpaired.val    = {f_componentG1,f_componentG2,m_TtestOneSample, c_statcomppermutation};
b_TtestUnpaired.help   = {'The implementation an unpaired t-test.'};

b_Ttestpaired        = cfg_branch;
b_Ttestpaired.tag    = 'b_Ttestpaired';
b_Ttestpaired.name   = 'Paired t-test' ;
b_Ttestpaired.val    = {f_componentG1,f_componentG2,m_TtestOneSample};
b_Ttestpaired.help   = {'The implementation a paired t-test. To make sure that the test is relevant, the two groups need to be paired data, which means the data can be naturally matched between the two groups.'};




f_anovan         = cfg_files;
f_anovan.name    = 'Group identification'; 
f_anovan.tag     = 'f_anovan';       %file names
f_anovan.filter  = {'xlsx','xls','txt'};
f_anovan.ufilter =  '.*';    
f_anovan.num     = [1 Inf];     % Number of inputs required 
f_anovan.help    = {'First row: dir',...
    'Second row: component file name',...
    'Third row and following for factor identification',... 
    'In the factor identification used 0 if the subject need to be rejected from the stat.'}; 

 

b_ANOVAN        = cfg_branch;
b_ANOVAN.tag    = 'b_ANOVAN';
b_ANOVAN.name   = 'By channel' ;
b_ANOVAN.val    = {f_anovan};
b_ANOVAN.help   = {'Apply the anonan on each channel of the xls list. Enter the observation to compute the anovan.',...
    'First row: dir',...
    'Second row: component file name',...
    'Third row and following for factor identification',... 
    'In the factor identification used 0 if the subject need to be rejected from the stat.'};


b_ANOVANzone        = cfg_branch;
b_ANOVANzone.tag    = 'b_ANOVANzone';
b_ANOVANzone.name   = 'By zone' ;
b_ANOVANzone.val    = {f_anovan};
b_ANOVANzone.help   = {'Apply the anonan on each specific region define by the xls file. Enter the observation to compute the anovan.',...
    'First row: dir',...
    'Second row: component file name',...
    'Third row for the zone file',...
    'Fourth row for the region label to use',...
    'Fifth and following for factor identification '};

c_ANOVAN         = cfg_choice;
c_ANOVAN.tag     = 'c_ANOVAN';
c_ANOVAN.name    = 'Anovan' ;
c_ANOVAN.values  = {b_ANOVAN,b_ANOVANzone};
c_ANOVAN.val     = {b_ANOVAN};
c_ANOVAN.help    = {'N-way analysis of variance, the anova could be performed by channel (channelist) or by zone.'};


f_componentexportzone         = cfg_files;
f_componentexportzone.name    = 'Component to export by zone'; 
f_componentexportzone.tag     = 'f_componentexportzone';       %file names
f_componentexportzone.filter  = {'xlsx','xls','txt'};
f_componentexportzone.ufilter =  '.*';    
f_componentexportzone.num     = [1 Inf];     % Number of inputs required 
f_componentexportzone.help    = {'First row: dir',...
    'Second row: component file name .mat file with each channel use for the topo',...
    'Third row and following for factor identification for further stat'}; 


f_componentexportch         = cfg_files;
f_componentexportch.name    = 'Component to export by channel'; 
f_componentexportch.tag     = 'f_componentexportch';       %file names
f_componentexportch.filter  = {'xlsx','xls','txt'};
f_componentexportch.ufilter =  '.*';    
f_componentexportch.num     = [1 Inf];     % Number of inputs required 
f_componentexportch.help    = {'First row: dir',...
    'Second row: component file name .mat file with each channel use for the topo',...
    'Third row channelList',...
    'Additional rows factor of identification for further stat'}; 

c_StatcomponentExport        = cfg_choice;
c_StatcomponentExport.tag     = 'c_StatcomponentExport';
c_StatcomponentExport.name    = 'Export' ;
c_StatcomponentExport.values  = {f_componentexportch, f_componentexportzone};
c_StatcomponentExport.val     = {f_componentexportch};
c_StatcomponentExport.help    = {'Create a mat file with all component for stat in external program'};


c_statcomponent         = cfg_choice; 
c_statcomponent.tag     = 'c_statcomponent';
c_statcomponent.name    = 'Choose the statistical test';
c_statcomponent.values  = {b_TtestOneSample, b_Ttestpaired,b_TtestUnpaired,c_ANOVAN,c_StatcomponentExport};
c_statcomponent.val     = {b_TtestOneSample}; %Default option
c_statcomponent.help    = {'Statistical tests are available: One sample t-test, Paired t-test, Unpaired t-test, Anovan (use matlab Statistics and Machine Learning Toolbox). The first step is to choose which test you want to perform.'};

E_statcomponent      = cfg_exbranch;
E_statcomponent.name = 'Stats components';
E_statcomponent.tag  = 'E_statcomponent';
E_statcomponent.val  = {c_statcomponent,m_statcomponent_saveoption,e_statcomponent_alpha, e_STATCOMPPath};
E_statcomponent.prog = @nirs_run_E_statcomponent;
E_statcomponent.vout = @nirs_cfg_vout_E_statcomponent;
E_statcomponent.help = {'Apply basic statistics on exported components.'};

function vout = nirs_cfg_vout_E_statcomponent(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end




% Executable Branch
E_MarkCorrectionAsNoise      = cfg_exbranch;
E_MarkCorrectionAsNoise.name = 'Mark Correction As Noise';
E_MarkCorrectionAsNoise.tag  = 'E_MarkCorrectionAsNoise';
E_MarkCorrectionAsNoise.val  = {NIRSmat, DelPreviousData};
E_MarkCorrectionAsNoise.prog = @nirs_run_E_MarkCorrectionAsNoise;
E_MarkCorrectionAsNoise.vout = @nirs_cfg_vout_E_MarkCorrectionAsNoise;
E_MarkCorrectionAsNoise.help = {'Identify intervall were a correction were apply and mark as noise (yellow) to show were the modification were don '};
%make NIRS.mat available as a dependency
function vout = nirs_cfg_vout_E_MarkCorrectionAsNoise(job)
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
    'If it is used before the normalization step, the block will be cut and associate to the trig in nirs data.',...
    'If it is used after the normalization step, you will need to match the blocks to the nirs data manually.'};


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
e_readEEGMarker_Marker.help    = {'Write the markers to import; use the exact label.',...
    'As an example: MARKER1:1,MARKER2=2,MARKER3=3'}; 


E_readEEGMarker      = cfg_exbranch;      
E_readEEGMarker.name = 'Import EEG marker' ;            
E_readEEGMarker.tag  = 'E_readEEGMarker'; 
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
Video_files.val  = {''};
Video_files.help    = {'Open Video or Audio files.',...
'The issue with 32 bit codec is that it needs 32 bits Matlab version.'}; % help text displayed


video_files_sync         = cfg_menu;
video_files_sync.name    = 'Ref Synchronisation trigger';
video_files_sync.tag     = 'video_files_sync';       
video_files_sync.labels  = {'EEG ','NIRS ','AUX '};
video_files_sync.values  = {0,1,2};
video_files_sync.val     = {0};
video_files_sync.help    = {'Video is recording simultaneously to EEG or NIRS data recording, if there is a delay please edit the video to removing exceeding part at the beginning.'};


b_videooffset_no         = cfg_branch;
b_videooffset_no.tag     = 'b_videooffset_no';
b_videooffset_no.name    = 'No';
b_videooffset_no.val     = {};
b_videooffset_no.help    = {'The start of video is synchronized with the triggers.'}';

i_videooffset         = cfg_entry;
i_videooffset.name    = 'Lag(second) = Video - Ref';
i_videooffset.tag     = 'i_videooffset';       
i_videooffset.strtype = 'r';
i_videooffset.num     = [1 inf];
i_videooffset.val     = {0}; 
i_videooffset.help    = {'As an example if the in video the event occur at 5 seconds and the trig in EEG happen at 10 you should adjust the lag to -5 as AUDIO(5)- EEG(10)=LAG(-5)'};

b_videooffset_yes         = cfg_branch;
b_videooffset_yes.tag     = 'b_videooffset_yes';
b_videooffset_yes.name    = 'Yes';
b_videooffset_yes.val     = {i_videooffset};
b_videooffset_yes.help    = {'There is a delay between the video and the trigger. If the video starts after, use a positive value, if the video starts before, use a negative value'}';

i_videooffsetwithtrig         = cfg_entry;
i_videooffsetwithtrig.name    = 'trig';
i_videooffsetwithtrig.tag     = 'i_videooffsetwithtrig';       
i_videooffsetwithtrig.strtype = 'r';
i_videooffsetwithtrig.num     = [1 inf];
i_videooffsetwithtrig.val     = {11}; 
i_videooffsetwithtrig.help    = {'Enter the trig that mark the lag identify for the recording'};


b_videooffsetwithtrig_yes         = cfg_branch;
b_videooffsetwithtrig_yes.tag     = 'b_videooffsetwithtrig_yes';
b_videooffsetwithtrig_yes.name    = 'Yes sync with a specific trigger';
b_videooffsetwithtrig_yes.val     = {i_videooffset, i_videooffsetwithtrig };
b_videooffsetwithtrig_yes.help    = {'There is a delay between the video and the trigger. If the video starts after, use a positive value, if the video starts before, use a negative value'}';


c_videooffset        = cfg_choice;
c_videooffset.tag    = 'c_videooffset';
c_videooffset.name   = 'Offset';
c_videooffset.values = {b_videooffset_yes, b_videooffset_no,b_videooffsetwithtrig_yes};
c_videooffset.val    = {b_videooffset_no};
c_videooffset.help   = {'If there is lag between Video and Ref Synchronisation trigger, select ''yes''. Otherwise, select ''No''.'}';

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
Audio_files_sync.help    = {'Audio is recorded simultaneously to EEG or NIRS data recording, if they have a delay please edit offset.'};


b_Audiooffset_no         = cfg_branch;
b_Audiooffset_no.tag     = 'b_Audiooffset_no';
b_Audiooffset_no.name    = 'No';
b_Audiooffset_no.val     = {};
b_Audiooffset_no.help    = {'The Audio start is synchronise with the triggers'}';

i_Audiooffset         = cfg_entry;
i_Audiooffset.name    = 'Lag(second) = Audio - Ref';
i_Audiooffset.tag     = 'i_Audiooffset';       
i_Audiooffset.strtype = 'r';
i_Audiooffset.num     = [1 inf];
i_Audiooffset.val     = {0}; 
i_Audiooffset.help    = {'As example if in audio the event occurs at 5 seconds and the trig in EEG happen at 10 you should adjust the lag to -5 as AUDIO(5)- EEG(10)=LAG(-5)'};

b_Audiooffset_yes         = cfg_branch;
b_Audiooffset_yes.tag     = 'b_Audiooffset_yes';
b_Audiooffset_yes.name    = 'Yes';
b_Audiooffset_yes.val     = {i_Audiooffset};
b_Audiooffset_yes.help    = {'They have a delay between the audio and the trigger. If the video starts after, use a positive value, if the audio starts before, use a negative value.'}';

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
E_readAudio.help = {'Read simultaneous audio data. The synchronization could be done using EEG AUX or NIRS trig, off'};

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


AUX_files_paired_NIRS         = cfg_menu;
AUX_files_paired_NIRS.name    = 'Sync option';
AUX_files_paired_NIRS.tag     = 'AUX_files_paired_NIRS';       
AUX_files_paired_NIRS.labels  = {'Defaults ','Paired one to one with NIRS file'};
AUX_files_paired_NIRS.values  = {0,1};
AUX_files_paired_NIRS.val     = {0};
AUX_files_paired_NIRS.help    = {'Usualy we do expect one AUX file for one nirs file but for some particular the sessions could be cut in two recording use the option paired one to one with nirs file if you enconter this cas, in that case you are limited to one AUX file by nirs file.'};



% Executable Branch
E_readAUX      = cfg_exbranch;      
E_readAUX.name = 'Read AUX' ;            
E_readAUX.tag  = 'E_readAUX'; 
E_readAUX.val  = {NIRSmat,AUX_files, AUX_files_paired_NIRS};   
E_readAUX.prog = @nirs_run_readAUX;  
E_readAUX.vout = @nirs_cfg_vout_readAUX;
E_readAUX.help = {'Read simultaneous AUX data; Note trig must be synchronized with the NIRS',... 
    'If it is used before the normalization step, the block will be cut and associate to the trig in nirs data.',...
    'If it is used after the normalization step, you will need to match the blocks to the nirs data manually.'};



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
m_replaceaux.help     = {'Replace aux will erase any reference to previous aux file, add aux will add the new hrf as an additional auxiliary file.',...
    'Please do not use the same name if you want to add a new one.'};

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
e_AUXdir.num = [0 Inf];  % Number of inputs required 
e_AUXdir.val  = {''};
%e_AUXdir.check = [];
e_AUXdir.help    = {'Auxiliary output folder. If empty the output folder will be the same as data nirs'}; % help text displayed


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

e_HRFpretime        = cfg_entry; 
e_HRFpretime.name    = 'Pretime (s)'; 
e_HRFpretime.tag     = 'e_HRFpretime';      
e_HRFpretime.strtype = 'r';       
e_HRFpretime.num     = [0 inf];     
e_HRFpretime.val     = {['']};
e_HRFpretime.help    = {'Pretime interval, let empty to use same as normalization pretime'};

e_HRFposttime        = cfg_entry; 
e_HRFposttime.name    = 'Posttime (s)'; 
e_HRFposttime.tag     = 'e_HRFposttime';      
e_HRFposttime.strtype = 'r';       
e_HRFposttime.num     = [0 inf];     
e_HRFposttime.val     = {['']};
e_HRFposttime.help    = {'Posttime interval, let empty to use same as normalization posttime'};

e_HRF_SDmodel        = cfg_entry;
e_HRF_SDmodel.tag    = 'e_HRF_SDmodel';
e_HRF_SDmodel.name   = 'Add SD zone to add physiological regressor in the model';
e_HRF_SDmodel.strtype = 's';      
e_HRF_SDmodel.num     = [0 Inf];     % Number of inputs required 
e_HRF_SDmodel.val  = {'No'};
e_HRF_SDmodel.help   = {['Add short distance in the GLM model, ' ...
    'If you write ''Global'', the zone Global.zone will be used as defaul regressor in your folder subject for the physiology using the average of all channel average.,...' ...
    'If you write a ''AllShortDistance.zone'' the default short distance file will be used in your folder subject ' ...
    'Use a specific file name if you want to specify you choice at a custom location' ...
    ' If you enter No, or let the field empty no physiological regressor will be add in the model']};


e_HRF_AUXmodel        = cfg_entry;
e_HRF_AUXmodel.tag    = 'e_HRF_AUXmodel';
e_HRF_AUXmodel.name   = 'Add AUX in the model';
e_HRF_AUXmodel.strtype = 's';   
e_HRF_AUXmodel.num     = [0 Inf];     % Number of inputs required 
e_HRF_AUXmodel.val  = {''};
e_HRF_AUXmodel.help   = {'Add auxilairy field in the model'};


b_HRFtriggeronset        = cfg_branch;
b_HRFtriggeronset.tag    = 'b_HRFtriggeronset';
b_HRFtriggeronset.name   = 'HRF trigger onset';
b_HRFtriggeronset.val    = {e_HRFtrigger, e_HRFduration, e_HRFpretime, e_HRFposttime,e_HRFlabel e_TimetoPeak1, e_FWHM1,e_TimetoPeak2,e_FWHM2,e_DIP,e_AUXdir,e_HRF_SDmodel}; %,,e_HRF_AUXmodel 
b_HRFtriggeronset.help   = {'Model HRF response using trigger onset and the user defines fix duration.'};


e_HRFxlsfiles         = cfg_files;  
e_HRFxlsfiles.name    = 'XLS onset files'; % The displayed name
e_HRFxlsfiles.tag     = 'e_HRFxlsfiles';          
e_HRFxlsfiles.num     = [0 Inf];     % Number of inputs required 
e_HRFxlsfiles.val{1}  = {''};
e_HRFxlsfiles.help    = {'Open excel or text file with HRF table onset information. This file must contain 3 columns, first column = onset, second column = duration and third column=weight',...
                        'fourth column pretime in second and fifth column postime in second', }; % help text displayed



b_HRFxlsonset        = cfg_branch;
b_HRFxlsonset.tag    = 'b_HRFxlsonset';
b_HRFxlsonset.name   = 'HRF xls onset';
b_HRFxlsonset.val    = {e_HRFxlsfiles  e_HRFlabel e_TimetoPeak1, e_FWHM1,e_TimetoPeak2,e_FWHM2,e_DIP,e_AUXdir e_HRF_SDmodel  };
b_HRFxlsonset.help   = {'Model HRF response using onset write the table (.xlsx or .txt file) using 3 columns first the onset time in seconds, second the duration in seconds and third the weight 1 or other'};




c_createAUXauto         = cfg_choice;
c_createAUXauto.tag     = 'c_createAUXauto';
c_createAUXauto.name    = 'Choose AUX to generate';
c_createAUXauto.values  = {b_HRFtriggeronset b_HRFxlsonset};
c_createAUXauto.val     = {b_HRFtriggeronset}; %Default option
c_createAUXauto.help    = {''};

% Executable Branch
E_createAUXauto      = cfg_exbranch;      
E_createAUXauto.name = 'Create AUX' ;            
E_createAUXauto.tag  = 'E_createAUXauto'; 
E_createAUXauto.val  = {NIRSmat,m_replaceaux,c_createAUXauto};   
E_createAUXauto.prog = @nirs_run_createAUXauto;  
E_createAUXauto.vout = @nirs_cfg_vout_createAUXauto;
E_createAUXauto.help = {'Read simultaneous AUX data; Note trig must be synchronized with the NIRS.',... 
    'If it is used before the normalization step, the block will be cut and associate to the trig in nirs data.',...
    'If it is used after the normalization step, you will need to match the blocks to the nirs data manually.'};



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
    'To create this list automatically, use the module Utility /Create channel list.'};

b_nodelist         =  cfg_choice; %Select NIRS.mat for this subject 
b_nodelist.name    = 'Node list'; % The displayed name
b_nodelist.tag     = 'b_nodelist';       %file names
b_nodelist.val     = {I_chcorrlist };
b_nodelist.values  = {I_chcorrlist I_zonecorrlist};  
b_nodelist.help    = {'Nodes could be defined by a region of interest (zones) or each channel.',...
    'Technical note: to allow subject comparison to ensure channels or zones spatial localization analogous.'};



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
i_RandomSample_crossspectrum.help    = {'Number of random samples to compute the cross spectrum, they are taken using pseudorandom segment in the whole block.'}; 

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
b_crossspectrum.help   = {'Coherence is a statistic representing the relationship between two signals and is also an extension of correlation to the frequency domain (Kida, 2016). Coherence is known as magnitude squared coherence is defined as the complex conjugate product of the Fourier transforms data X(f)* Y*T(f). x(t) and y(t) are two time series, Gxy(f) is the cross-spectral density between x and y, and Gxx(f) and Gyy(f) are the auto spectral densities of x and y, respectively. The coherence is implemented to use one long continuous segment of the recording. In case you record multiple sessions, you may join them using the concatenate module. A large number of segments (Number of random samples) of duration (Length of the segment) will be picked randomly (circular bootstrap). Any segments that belong to a specific artifact period will be excluded from the coherence calculation. The segment will be randomly segmented to calculate coherence based on many segments. An FFT is computed on each random segment and the coherence is measured based on the specified frequency range (The frequency range to obtain Cxy(f)) to obtain a connectivity matrix representative of the whole recording. '};


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
edownsample_Granger.help    = {['Downsample factor useful if you want to do granger causality efficiently on slow frequency signal ']}; 


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
                
RespirationBBM           = cfg_menu;
RespirationBBM.tag       = 'RespirationBBM';
RespirationBBM.name      = 'Respiration beat by minute';
RespirationBBM.labels    = {'Yes','No'};
RespirationBBM.values    = {1, 0};
RespirationBBM.val       = {1};
RespirationBBM.help      = {'Mesure in the auxilairy (Resp) channel to save in the structure the rate of respiration by minute. It could be use later to help to determine the state of the participant. '};
         
                     
b_PearsonBootstrap         = cfg_branch;
b_PearsonBootstrap.tag     = 'b_PearsonBootstrap';
b_PearsonBootstrap.name    = 'Circular bootstrap';
b_PearsonBootstrap.val     = {i_TrialLenght_crossspectrum,i_RandomSample_crossspectrum,i_OutlierControl_crossspectrum, RespirationBBM};
b_PearsonBootstrap.help    = {'Use circular bootstrap to compute cross-correlation analysis.'};
                
m_Pearson           = cfg_menu;
m_Pearson.tag       = 'm_Pearson';
m_Pearson.name      = 'Segments';
m_Pearson.labels    = {'Zero lag Cross-Correlation (Pearson)'};
m_Pearson.values    = {1};
m_Pearson.val       = {1};
m_Pearson.help      = {'Apply the zero lag Cross-Correlation on each segmented file of the data set.'};


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
b_Pearson.help   = {'Compute zeros lag cross-correlation between channels or region of interest.'};                  



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
b_HilbertBootstrap.help    = {'Use circular bootstrap to segment randomly multiple segments in the fNIRS files to compute the Hilbert transform.'};



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
b_Hilbert.help   = {'See B. Molavi et al., �Analyzing the resting state functional connectivity in the human language system using near infrared spectroscopy.�, Front. Hum. Neurosci. 7(January), 921 (2013) [doi:10.3389/fnhum.2013.00921].',...
                    'Use function from the CircStat2012a toolbox.' };



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
b_waveletcluster.name   = 'Wavelet Coherence NOT TESTED';
b_waveletcluster.val    = {e_startwaveletcluster,e_stopwaveletcluster,i_Freq_crossspectrum, e_widthwaveletcluster };
b_waveletcluster.help   = {''};


I_chcorrlist_type         = cfg_choice;
I_chcorrlist_type.tag     = 'I_chcorrlist_type';
I_chcorrlist_type.name    = 'Connectivity to use';
I_chcorrlist_type.val     = {b_crossspectrum};
I_chcorrlist_type.help    = {''};
I_chcorrlist_type.values  = {b_Pearson, b_Hilbert,b_crossspectrum, b_waveletcluster};%,b_Granger,b_Phase,b_waveletcluster};%,b_Hilbert,b_Granger, b_Phase, b_crossspectrum
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
I_ConnectivityMATName.help    = {'OutputName for the matrix. Use an identifier as subject number.'};


E_chcorrMat      = cfg_exbranch;      
E_chcorrMat.name = 'Connectivity matrix' ;            
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
zonemat.help    = {'Select zone to create channel list for verification, zone label as well as detector source combination will be listed.'}; % help text displayed

% Executable Branch
E_zone2channellist      = cfg_exbranch;      
E_zone2channellist.name = 'Transfer .zone into a .txt list' ;            
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
zonetxt.help    = {'Select .txt list that contains zonelabel  to create channel list for verification, zone label as well as detector source combination will be listed.'}; % help text displayed

% Executable Branch
E_channellist2zone      = cfg_exbranch;      
E_channellist2zone.name = 'Transfer a .txt list into .zone' ;            
E_channellist2zone.tag  = 'E_channellist2zone'; 
E_channellist2zone.val  = {NIRSmat,zonetxt, m_zonematformat};   
E_channellist2zone.prog = @nirs_run_E_channellist2zone;  
E_channellist2zone.vout = @nirs_cfg_vout_channellist2zone;
E_channellist2zone.help = {'Create a zone using list of detectors and sources belonging to this zone, use Label: zonename ;  RGBcolor: 255 0 0 ;  then a line for each Detector Source combination D01 E1.',...
    'The zone will be created at the NIRS.mat location to be associate to this subject.'};

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
m_fishertransform.labels    = {'Yes', 'No','Yes & valeur absolu'};
m_fishertransform.values    = {1,2,3};
m_fishertransform.val       = {1};
m_fishertransform.help      = {'Use the fisher transform 1/2 ln((1+p)/(1-p)), when the transformation is applied to the sample correlation coefficient, the sampling distribution of the resulting variable is approximately normal, with a variance that is stable over different values of the underlying true correlation.'};


m_nodeunit           = cfg_menu;
m_nodeunit.tag       = 'm_nodeunit';
m_nodeunit.name      = 'Nodes';
m_nodeunit.labels    = {'Channels', 'Zones'};
m_nodeunit.values    = {1,2};
m_nodeunit.val       = {1};
m_nodeunit.help      = {'Apply the statistics on each node. Define nodes as each channel or each average of zone channels.'};

%use for stat matrice
e_minnbchan         =  cfg_entry;
e_minnbchan.name    = 'minnbchan';
e_minnbchan.tag     = 'e_minnbchan';       
e_minnbchan.strtype = 's';
e_minnbchan.num     = [0 inf];
e_minnbchan.val     = {'4'};
e_minnbchan.help    = { 'Minnbchan number is a parameter to define minimal number of commun neighbord to be include in the neighboring.' };

e_neighbourdist         =  cfg_entry; %use for matrix stat
e_neighbourdist.name    = 'Neighbor';
e_neighbourdist.tag     = 'e_neighbourdist';       
e_neighbourdist.strtype = 's';
e_neighbourdist.num     = [0 inf];
e_neighbourdist.val     = {'3'};
e_neighbourdist.help    = {'Define how neighbor will be define:',...
    'Write ''link'' to consider optode with commun link together.',...
    'Write a number such as ''5'' to define the distance between optode to be include in a neighboring to define a cluster.',... 
    'All optode with a distance lower will be consider as neighbors.',...
    'Write ''all'' to use all possible other link as neighbors, that will result to one big cluster'};

b_MCT_ClusterBased         = cfg_branch;
b_MCT_ClusterBased.tag    = 'b_MCT_ClusterBased';
b_MCT_ClusterBased.name   = 'Cluster Based' ;
b_MCT_ClusterBased.val    = {e_neighbourdist, e_minnbchan};
b_MCT_ClusterBased.help   = {'Cluster based permutation'};

b_MCT_Max         = cfg_branch;
b_MCT_Max.tag    = 'b_MCT_Max';
b_MCT_Max.name   = 'No' ;
b_MCT_Max.val    = {};
b_MCT_Max.help   = {'Use the maximal univariate comparaison to create null distribution; Troendle, J.F., 1995. A Stepwise Resampling Method of Multiple Hypothesis Testing. Journal of the American Statistical Association 90, 370�378. https://doi.org/10.1080/01621459.1995.10476522'};

 

c_statMultipleComparaisonTesting        = cfg_choice;
c_statMultipleComparaisonTesting.tag     = 'c_statMultipleComparaisonTesting';
c_statMultipleComparaisonTesting.name    = 'Multiple Comparaison';
c_statMultipleComparaisonTesting.values  = {b_MCT_Max, b_MCT_ClusterBased};
c_statMultipleComparaisonTesting.val     = {b_MCT_ClusterBased}; %Default option
c_statMultipleComparaisonTesting.help    = {'Cluster-based permutation adapt clusterstat function from fieldtrip function and see for reference Maris, E., Oostenveld, R., 2007. Nonparametric statistical testing of EEG- and MEG-data. Journal of Neuroscience Methods 164, 177�190. https://doi.org/10.1016/j.jneumeth.2007.03.024'};



e_npermutation         =  cfg_entry;
e_npermutation.name    = 'Nb permutations';
e_npermutation.tag     = 'e_npermutation';       
e_npermutation.strtype = 's';
e_npermutation.num     = [0 inf];
e_npermutation.val     = {'500'};
e_npermutation.help    = {'Define the number of non-parametric permutations to use.'};

b_permutation        = cfg_branch;
b_permutation.tag    = 'b_permutation';
b_permutation.name   = 'Yes' ;
b_permutation.val    = {e_npermutation, c_statMultipleComparaisonTesting};
b_permutation.help   = {'Apply permutation test'};


b_Nopermutation        = cfg_branch;
b_Nopermutation.tag    = 'b_Nopermutation';
b_Nopermutation.name   = 'No' ;
b_Nopermutation.val    = {};
b_Nopermutation.help   = {'Do not apply permutation test'};


c_statpermutation        = cfg_choice;
c_statpermutation.tag     = 'c_statpermutation';
c_statpermutation.name    = 'Univariate permutations';
c_statpermutation.values  = {b_permutation,b_Nopermutation };
c_statpermutation.val     = {b_Nopermutation}; %Default option
c_statpermutation.help    = {'Use permutation to create the empirical null distribution to calculate.',...
    'Compared each link perm to its permuted distribution and also the max link value among the comparaison multiple comparaison to be more restrictive. See: Lage-Castellanos, A., Mart�nez-Montes, E., Hern�ndez-Cabrera, J.A., Gal�n, L., 2010. False discovery rate and permutation test: An evaluation in ERP data analysis. Statistics in Medicine 29, 63�74. https://doi.org/10.1002/sim.3784'};

m_export_matrix        = cfg_menu;
m_export_matrix.tag    = 'm_export_matrix';
m_export_matrix.name   = 'Export';
m_export_matrix.labels = {'.mat format' };
m_export_matrix.values = {1};
m_export_matrix.val    = {1}; 
m_export_matrix.help   = {'Export all observation in one file for external statistique.'};

m_TtestOneSample_matrix        = cfg_menu;
m_TtestOneSample_matrix.tag    = 'm_TtestOneSample_matrix';
m_TtestOneSample_matrix.name   = 'Tailed ';
m_TtestOneSample_matrix.labels = {'2-tailed','1-tailed left (neg)' ,'1-tailed right (pos)' };
m_TtestOneSample_matrix.values = {1,2,3};
m_TtestOneSample_matrix.val    = {1}; 
m_TtestOneSample_matrix.help   = {'Use one-tailed or two-tailed t-test.'};

f_matrix         = cfg_files;
f_matrix.name    = 'Enter list connectivity matrix'; 
f_matrix.tag     = 'f_matrix';       %file namesH
f_matrix.filter  = {'xlsx','xls','txt'};
f_matrix.ufilter = '.*';    
f_matrix.num     = [1 Inf];     % Number of inputs required 
f_matrix.help    = {'Enter the list of connectivity matrix to test statistically.',...
    'column 1 dir, column 2 name matrices,  column 3 zone id , column 4 group information label by number 1,2,3.... use 0 to exclude participant, following column be use any other covariables. '}; 

e_TtestOneSample_meanvalue         = cfg_entry;
e_TtestOneSample_meanvalue.name    = 'Hypothetical mean value';
e_TtestOneSample_meanvalue.tag     = 'e_TtestOneSample_meanvalue';       
e_TtestOneSample_meanvalue.strtype = 'r';
e_TtestOneSample_meanvalue.num     = [1 Inf];     
e_TtestOneSample_meanvalue.val     = {0}; 
e_TtestOneSample_meanvalue.help    = {'Enter the hypothetical mean value'};



e_TtestOneSampleGR         = cfg_entry; %path
e_TtestOneSampleGR.name    = 'Group identification';
e_TtestOneSampleGR.tag     = 'e_TtestOneSampleGR';       
e_TtestOneSampleGR.strtype = 'r';
e_TtestOneSampleGR.num     = [1 Inf];     
e_TtestOneSampleGR.val     = {1}; 
e_TtestOneSampleGR.help    = {['Enter the group to apply the one sample t-test (refer to group value in the xls file).']}; 

b_TtestOneSamplematrix        = cfg_branch;
b_TtestOneSamplematrix.tag    = 'b_TtestOneSamplematrix';
b_TtestOneSamplematrix.name   = 'One-sample t-test' ;
b_TtestOneSamplematrix.val    = {m_TtestOneSample_matrix,e_TtestOneSample_meanvalue,e_TtestOneSampleGR};
b_TtestOneSamplematrix.help   = {'The one-sample t-test is a parametric test of the location parameter when the population standard deviation is unknown.',
    'T = mean(X)/(STD(x)*sqrt(n))'};

e_TtestOneSampleGR2         = cfg_entry; %path
e_TtestOneSampleGR2.name    = 'Group 2 identification';
e_TtestOneSampleGR2.tag     = 'e_TtestOneSampleGR2';       
e_TtestOneSampleGR2.strtype = 'r';
e_TtestOneSampleGR2.num     = [1 Inf];     
e_TtestOneSampleGR2.val     = {2}; 
e_TtestOneSampleGR2.help    = {['Enter the second group compared (refer to group value in the xls file).']}; 


b_PermutationTest         = cfg_branch;
b_PermutationTest.tag    = 'b_PermutationTest';
b_PermutationTest.name   = 'Unpaired permutation t-test' ;
b_PermutationTest.val    = {e_npermutation, e_TtestOneSampleGR, e_TtestOneSampleGR2};
b_PermutationTest.help   = {'Compared 2 groups using permutation '};

b_UnpairedTtest        = cfg_branch;
b_UnpairedTtest.tag    = 'b_UnpairedTtest';
b_UnpairedTtest.name   = 'Unpaired t-test' ;
b_UnpairedTtest.val    = { e_TtestOneSampleGR, e_TtestOneSampleGR2, c_statpermutation};%m_TtestOneSample_matrix,
b_UnpairedTtest.help   = {'Compute Unpaired Ttest using using group in the xls. Use or not permutation test. Without permutation Cohen d size effet and tvalue will be computed'};

b_PairedTtest        = cfg_branch;
b_PairedTtest.tag    = 'b_PairedTtest';
b_PairedTtest.name   = 'Paired t-test' ;
b_PairedTtest.val    = { m_TtestOneSample_matrix, e_TtestOneSample_meanvalue, e_TtestOneSampleGR, e_TtestOneSampleGR2, c_statpermutation};
b_PairedTtest.help   = {'Compute Paired Ttest using repeated measures identify subject in the xls.  The subject list must be placed in order. The first subject in group 1 will be paired with the first subject in group 2 in the subject list.  The subtraction of each subject will be saved in the results folder.'};


b_exportNBSformat        = cfg_branch;
b_exportNBSformat.tag    = 'b_exportNBSformat';
b_exportNBSformat.name   = 'Export NBS format' ;
b_exportNBSformat.val    = {};
b_exportNBSformat.help   = {'Export to NBS network based statistic format. A Zalesky(2010) [doi: 10.1016/j.neuroimage.2010.06.041]. '};

e_GLMGR         = cfg_entry; %path
e_GLMGR.name    = 'Group to use';
e_GLMGR.tag     = 'e_GLMGR';       
e_GLMGR.strtype = 's';
e_GLMGR.num     = [1 Inf];     
e_GLMGR.val     = {'all'}; 
e_GLMGR.help    = {['Enter the group of the subject to use, subject no in this list will just be exclude keep all to use all subject.']}; 


b_Covariable_Mat         =  cfg_entry;
b_Covariable_Mat.name    = 'Covariable';
b_Covariable_Mat.tag     = 'b_Covariable_Mat';       
b_Covariable_Mat.strtype = 's';
b_Covariable_Mat.num     = [0 inf];
b_Covariable_Mat.val     = {'Name column'};
b_Covariable_Mat.help    = {'Example: Constant, Age. Use the exact column title to reconize which covariable to correlate with the connectivity scores. Separate by a comma when they are many covariables to explore. Add a constant in your model. '};


b_PearsonCorr_Mat        = cfg_branch;
b_PearsonCorr_Mat.tag    = 'b_PearsonCorr_Mat';
b_PearsonCorr_Mat.name   = 'Pearson correlation' ;
b_PearsonCorr_Mat.val    = {b_Covariable_Mat};
b_PearsonCorr_Mat.help   = {'Compute the pearson correlation coefficient between connectivity score and one covariable in the excel file.'};

b_substractidCovariable_Mat         =  cfg_entry;
b_substractidCovariable_Mat.name    = 'Covariable to subtract';
b_substractidCovariable_Mat.tag     = 'b_substractidCovariable_Mat';       
b_substractidCovariable_Mat.strtype = 's';
b_substractidCovariable_Mat.num     = [0 inf];
b_substractidCovariable_Mat.val     = {'Name column'};
b_substractidCovariable_Mat.help    = {'Subtract one covariable and save residual.'};


b_GLM_Mat        = cfg_branch;
b_GLM_Mat.tag    = 'b_GLM_Mat';
b_GLM_Mat.name   = 'GLM' ;
b_GLM_Mat.val    = {e_GLMGR,b_Covariable_Mat b_substractidCovariable_Mat,c_statpermutation};
b_GLM_Mat.help   = {'Apply the general linear model, Specify covariable as regressor y=b1*x1+b2*x2+...+c, do not forget to include a covariable for the constant. Use the function regress.m in matlab'};

b_LME_formula         =  cfg_entry;
b_LME_formula.name    = 'Formula';
b_LME_formula.tag     = 'b_LME_formula';       
b_LME_formula.strtype = 's';
b_LME_formula.num     = [0 inf];
b_LME_formula.val     = {'columnlabel~MAT+(1+id)'};
b_LME_formula.help    = {'Apply LME formula using wilkinson notation, please look on matlab documentation fitlme to see the formula example',...
                        'avoid space in column label to enter column label directly in the formula definition as example ',... 
                        'y ~ X1 + X2       => Fixed effects for the intercept, X1 (column label) and X2 (column label). This is equivalent to y ~ 1 + X1 + X2.',...
                        'y ~ -1 + X1 + X2  => No intercept and fixed effects for X1 and X2. The implicit intercept term is suppressed by including -1',... 
                        'y ~ 1 + (1 | g1)  => Fixed effects for the intercept plus random effect for the intercept for each level of the grouping variable g1.',...
                        'y ~ X1 + (X1 | g1)=> Random intercept and slope, with possible correlation between them. This is equivalent to y ~ 1 + X1 + (1 + X1|g1)',...
                        'y ~ X1 + (1 | g1) + (-1 + X1 | g1) => Independent random effects terms for intercept and slope.',...
                        'y ~ 1 + (1 | g1) + (1 | g2) + (1 | g1:g2) =>Random intercept model with independent main effects for g1 and g2, plus an independent interaction effect.' };
                       
b_LME_Mat        = cfg_branch;
b_LME_Mat.tag    = 'b_LME_Mat';
b_LME_Mat.name   = 'LME' ;
b_LME_Mat.val    = {e_GLMGR, b_LME_formula,c_statpermutation};
b_LME_Mat.help   = {'Fit linear mixed-effects model function FitLME.m '};

b_fitLM_Mat        = cfg_branch;
b_fitLM_Mat.tag    = 'b_LM_Mat';
b_fitLM_Mat.name   = 'LM' ;
b_fitLM_Mat.val    = {e_GLMGR, b_LME_formula,c_statpermutation};
b_fitLM_Mat.help   = {'Fit linear model function FitLM.m '};

e_GRzscore         = cfg_entry; %path
e_GRzscore.name    = 'Subject to apply zscore (group)';
e_GRzscore.tag     = 'e_GRzscore';       
e_GRzscore.strtype = 'r';
e_GRzscore.num     = [1 Inf];     
e_GRzscore.val     = {2}; 
e_GRzscore.help    = {['Enter the group identification to compute the zscore for each subject according to the group distribution ']};

b_zscore_Mat        = cfg_branch;
b_zscore_Mat.tag    = 'b_zscore_Mat';
b_zscore_Mat.name   = 'Zscore' ;
b_zscore_Mat.val    = {e_TtestOneSampleGR e_GRzscore};
b_zscore_Mat.help   = {'Find the normal distribution group 1 and compute zscore matrice for each individual belonging to group 2o a groupe 1'};



e_Anova1GR         = cfg_entry; %path
e_Anova1GR.name    = 'Group identification';
e_Anova1GR.tag     = 'e_Anova1GR';       
e_Anova1GR.strtype = 's';
e_Anova1GR.num     = [1 Inf];     
e_Anova1GR.val     = {'1, 2, 3'}; 
e_Anova1GR.help    = {['Enter the group to include in the ANOVA (refer to group value in the xls file).']};


b_anova1_Mat        = cfg_branch;
b_anova1_Mat.tag    = 'b_anova1_Mat';
b_anova1_Mat.name   = 'Anova one way' ;
b_anova1_Mat.val    = {e_Anova1GR, c_statpermutation};
b_anova1_Mat.help   = {'Find One-way analysis of variance anova1 list the group to evaluate in the anova fdr correction. Eta2 size effect are compute using Harald Hentschke (2022). hhentschke/measures-of-effect-size-toolbox (https://github.com/hhentschke/measures-of-effect-size-toolbox), GitHub. Retrieved November 16, 2022.  '};

b_ANCOVA_Covariable         =  cfg_entry;
b_ANCOVA_Covariable.name    = 'Covariable (only one)';
b_ANCOVA_Covariable.tag     = 'b_Covariable_Mat';       
b_ANCOVA_Covariable.strtype = 's';
b_ANCOVA_Covariable.num     = [0 inf];
b_ANCOVA_Covariable.val     = {'Name column'};
b_ANCOVA_Covariable.help    = {'Use the exact column title to recognize which covariable.'};

b_ANCOVA_Mat        = cfg_branch;
b_ANCOVA_Mat.tag    = 'b_ANCOVA_Mat';
b_ANCOVA_Mat.name   = 'ANCOVA (AOCtool)' ;
b_ANCOVA_Mat.val    = {e_Anova1GR,b_ANCOVA_Covariable};
b_ANCOVA_Mat.help   = {'ANCOVA (Analysis of covariance) use the matlab function aoctool' };


e_fitMANCOVAN_GR         = cfg_entry; %path
e_fitMANCOVAN_GR.name    = 'Group identification';
e_fitMANCOVAN_GR.tag     = 'e_fitMANCOVAN_GR';       
e_fitMANCOVAN_GR.strtype = 's';
e_fitMANCOVAN_GR.num     = [1 Inf];     
e_fitMANCOVAN_GR.val     = {'1, 2, 3'}; 
e_fitMANCOVAN_GR.help    = {['Enter the group to include in the ANOVA (refer to group value in the xls file).']};

b_fitMANCOVAN_Covariable         =  cfg_entry;
b_fitMANCOVAN_Covariable.name    = 'Covariable';
b_fitMANCOVAN_Covariable.tag     = 'b_fitMANCOVAN_Covariable';       
b_fitMANCOVAN_Covariable.strtype = 's';
b_fitMANCOVAN_Covariable.num     = [0 inf];
b_fitMANCOVAN_Covariable.val     = {'Name column'};
b_fitMANCOVAN_Covariable.help    = {'Use the exact column title to recognize which covariable.'};



b_fitMANCOVAN_Mat        = cfg_branch;
b_fitMANCOVAN_Mat.tag    = 'b_fitMANCOVAN_Mat';
b_fitMANCOVAN_Mat.name   = 'ANCOVA (analysis of covariance with covariate)' ;
b_fitMANCOVAN_Mat.val    = {e_fitMANCOVAN_GR,b_fitMANCOVAN_Covariable, c_statpermutation};
b_fitMANCOVAN_Mat.help   = {'ANCOVA : Analysis of covariance (ANCOVA) is a general linear model which blends ANOVA and regression.',...
    'ANCOVA evaluates whether the means of a dependent variable (DV) are equal across levels of a categorical independent variable (IV) groups',...
    'while statistically controlling for the effects of other continuous variables that are not of primary interest, known as covariates (CV) or nuisance variables.',... 
    'equivalent to SPSS univariate linear general model with fix factor of groupe and selected covariates.',...
    'Use: William Gruner (2022). MANCOVAN (https://www.mathworks.com/matlabcentral/fileexchange/27014-mancovan), MATLAB Central File Exchange. Retrieved October 24, 2022. '};


%%%%%%%%%%%%%%%%%%%%%%%%%
e_Anovarep_withinsubjectGR       =  cfg_entry;
e_Anovarep_withinsubjectGR.name    = 'Repeated measures (group G)';
e_Anovarep_withinsubjectGR.tag     = 'e_Anovarep_withinsubjectGR';       
e_Anovarep_withinsubjectGR.strtype = 's';
e_Anovarep_withinsubjectGR.num     = [1 inf];
e_Anovarep_withinsubjectGR.val     = {'1,2'};
e_Anovarep_withinsubjectGR.help    = {'Matrice will be repeated measure, use xls table as in the paired t-test ensure that Groupe 1',...
    'is the first observation for the subject while Groupe 2 is the second observation and follow'};

e_Anovarep_predictor         =  cfg_entry;
e_Anovarep_predictor.name    = 'Predictor measures (X1, X2)';
e_Anovarep_predictor.tag     = 'e_Anovarep_predictor';       
e_Anovarep_predictor.strtype = 's';
e_Anovarep_predictor.num     = [0 inf];
e_Anovarep_predictor.val     = {'Name column'};
e_Anovarep_predictor.help    = {'Use the exact column title to recognize which covariable, Such as AGE or Gender'};



e_anovarep_model        =  cfg_entry;
e_anovarep_model.name    = 'Model';
e_anovarep_model.tag     = 'e_anovarep_model';       
e_anovarep_model.strtype = 's';
e_anovarep_model.num     = [0 inf];
e_anovarep_model.val     = {'G1,G2 ~ X1'};
e_anovarep_model.help    = {'Use Wilkinson notation to define your model see help in matlab for the function fitrm or see :',...
    'Wilkinson, G.N., Rogers, C.E., 1973. Symbolic Description of Factorial Models for Analysis of Variance.',...
    'Journal of the Royal Statistical Society. Series C (Applied Statistics) 22, 392�399.',...
    'https://doi.org/10.2307/2346786',...
    'The terms used must be column name in the xls table to specify your design.Example G1,G2 ~ AGE,  '};

b_anovarep_Mat        = cfg_branch;
b_anovarep_Mat.tag    = 'b_anovarep_Mat';
b_anovarep_Mat.name   = 'Anova repeted measure' ;
b_anovarep_Mat.val    = {e_Anovarep_withinsubjectGR,e_Anovarep_predictor, e_anovarep_model};
b_anovarep_Mat.help   = {'Apply Anova repeted measure model, fonction fitrm. ',...
    'You must paired subject as repeted measure in the groupe identification,',...
    'first observation of group 1 will be paired with first observation of group 2. '};

b_anovarep_Mat        = cfg_branch;
b_anovarep_Mat.tag    = 'b_anovarep_Mat';
b_anovarep_Mat.name   = 'Anova repeted measure' ;
b_anovarep_Mat.val    = {e_Anovarep_withinsubjectGR,e_Anovarep_predictor, e_anovarep_model};
b_anovarep_Mat.help   = {'Apply Anova repeted measure model, fonction fitrm. ',...
    'You must paired subject as repeted measure in the groupe identification,',...
    'first observation of group 1 will be paired with first observation of group 2. '};


%%%%%%%%%%%%%%%%%%%%%%%%%
e_NWayANOVA_model_GR       =  cfg_entry;
e_NWayANOVA_model_GR.name    = 'Group to be include in the model';
e_NWayANOVA_model_GR.tag     = 'e_NWayANOVA_model_GR';       
e_NWayANOVA_model_GR.strtype = 's';
e_NWayANOVA_model_GR.num     = [1 inf];
e_NWayANOVA_model_GR.val     = {'1,2'};
e_NWayANOVA_model_GR.help    = {'Write the group number in the column 4 of the matrix file identification to be consider in the model'};



e_NWayANOVA_model        =  cfg_entry;
e_NWayANOVA_model.name    = 'Model';
e_NWayANOVA_model.tag     = 'e_NWayANOVA_model';       
e_NWayANOVA_model.strtype = 's';
e_NWayANOVA_model.num     = [0 inf];
e_NWayANOVA_model.val     = {['{''PREMA'',''VERBAL'',''RARE''}']};
e_NWayANOVA_model.help    = {'Write the variable to include in the Nway anova model, it must be exact the column name. See anovan help in matlab'};



b_NWayANOVA_Mat        = cfg_branch;
b_NWayANOVA_Mat.tag    = 'b_NWayANOVA_Mat';
b_NWayANOVA_Mat.name   = 'NWAY Anova' ;
b_NWayANOVA_Mat.val    = {e_NWayANOVA_model_GR  e_NWayANOVA_model };
b_NWayANOVA_Mat.help   = {'Apply NWAY ANOVA'};

b_kruskalwallis_Mat        = cfg_branch;
b_kruskalwallis_Mat.tag    = 'b_kruskalwallis_Mat';
b_kruskalwallis_Mat.name   = 'Kruskal Wallis' ;
b_kruskalwallis_Mat.val    = {e_Anova1GR, c_statpermutation };
b_kruskalwallis_Mat.help   = {'Apply kruskalwallis, For permutation: Odiase, J.I., Ogbonmwan, S.M., 2005. JMASM20: Exact Permutation Critical Values For The Kruskal-Wallis One-Way ANOVA. J. Mod. App. Stat. Meth. 4, 609�620. https://doi.org/10.22237/jmasm/1130804820'};


b_manova1_Covariable         =  cfg_entry;
b_manova1_Covariable.name    = 'Covariable ';
b_manova1_Covariable.tag     = 'b_Covariable_Mat';       
b_manova1_Covariable.strtype = 's';
b_manova1_Covariable.num     = [0 inf];
b_manova1_Covariable.val     = {'Name column'};
b_manova1_Covariable.help    = {'Use the exact column title to recognize which covariable. Need toolbox '};

b_manova1_Mat        = cfg_branch;
b_manova1_Mat.tag    = 'b_manova1_Mat';
b_manova1_Mat.name   = 'manova (manova1)' ;
b_manova1_Mat.val    = {e_Anova1GR, b_manova1_Covariable};
b_manova1_Mat.help   = {'Apply manova on specific groups.'};

c_statmatrix         = cfg_choice;
c_statmatrix.tag     = 'c_statmatrix';
c_statmatrix.name    = 'Choose the statistical test';
c_statmatrix.values  = {m_export_matrix, b_TtestOneSamplematrix,b_UnpairedTtest, b_PairedTtest, b_PearsonCorr_Mat, b_GLM_Mat,b_zscore_Mat,b_anova1_Mat,b_anovarep_Mat,b_NWayANOVA_Mat, b_kruskalwallis_Mat, b_fitMANCOVAN_Mat, b_exportNBSformat,b_PermutationTest, b_LME_Mat, b_fitLM_Mat};%b_ANCOVA_Mat aoctool remove,
c_statmatrix.val     = {b_TtestOneSamplematrix}; %Default option
c_statmatrix.help    = {'Select one of the statistical tests.'};

% Folder selector
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
E_statmatrix.val  = {f_matrix,m_fishertransform,e_statcomponent_alpha,m_nodeunit,c_statmatrix, e_statmatrixPath};
E_statmatrix.prog = @nirs_run_E_statmatrix;
E_statmatrix.vout = @nirs_cfg_vout_E_statmatrix;
E_statmatrix.help = {'Apply basic statistics on exported connectivity matrices.'};

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
M_readMultimodal.help   = {'Read EEG or Auxiliary format to do visualized multimodal NIRS in the GUI'};

%Module Read data 
M_readNIRS        = cfg_choice; 
M_readNIRS.name   = 'Read data';
M_readNIRS.tag    = 'M_readNIRS';
M_readNIRS.values = {E_readNIRxscout,E_rawhomer, E_readSNIRF,boxy1, E_GenericDataExportBV, E_readNIRSport, M_readMultimodal}; 
M_readNIRS.help   = {'These modules read NIRS data in different formats.'};

%Module segment or concatenate data 
M_Segment        =   cfg_choice; 
M_Segment.name   =  'Segment/Onset';
M_Segment.tag    = 'M_Segment';
M_Segment.values = {segment ,E_Concatenate_nirsmat,E_Concatenate_file,E_aux2manualtrig,E_manualtriglsl, E_manualtrig,E_createonset_correlationsignal }; 
M_Segment.help   = {'These modules segment or combine data.'};


%Module 3 Preprocessing NIRS data (3a,3b,3c,3d,3e,3f)
M_preprocessing        = cfg_choice; 
M_preprocessing.name   = 'Preprocessing NIRS data';
M_preprocessing.tag    = 'M_preprocessing'; 
M_preprocessing.values = {E_artefactdetection E_chcardiaccontrol E_normalization E_filter E_prewhitening ODtoHbOHbR  E_detrend E_average E_nullifybad}; %preprocess enlever fast preprocessing JT
M_preprocessing.help   = {'These modules apply basic operations on fNIRS data'};


%Module Data Display GUI'
M_GUI        = cfg_choice;
M_GUI.name   = 'Data Display';
M_GUI.tag    = 'M_GUI';
M_GUI.values = {E_GUI E_VIDEO }; 
M_GUI.help   = {'Graphical user interface for data display.'};

%Module Component
M_dataComponent        = cfg_choice;
M_dataComponent.name   = 'Decomposition';
M_dataComponent.tag    = 'M_dataComponent';
M_dataComponent.values = { E_extractcomponent E_substractcomponent E_exportcomponent_list E_exportcomponent_zone E_statcomponent E_MarkCorrectionAsNoise}; 
M_dataComponent.help   = {'Apply operation on selected component identify in display GUI data decomposition'};

%Module Connectivity
M_Connectivity        = cfg_choice; 
M_Connectivity.name   = 'Connectivity';
M_Connectivity.tag    = 'M_Connectivity';
M_Connectivity.values = { E_chcorrMat E_GUI_lookmatrices E_statmatrix}; %
M_Connectivity.help   = {'Connectivity function.'};

%Module  Write external file
M_datawritenirs        =  cfg_choice; 
M_datawritenirs.name   = 'Write external file';
M_datawritenirs.tag    = 'M_datawritenirs';
M_datawritenirs.values = {E_writeNIRSHomer,  E_writeSNIRF, E_NIR_segment,E_writeHMR}; 
M_datawritenirs.help   = {'These modules convert nir file in .nirs, last module of data export, support the export of field such as data (d),coordinate (SD), trigger (s), time (t),do not support noise artifact marking or aux export'};





%Module  Write external file
M_others        =  cfg_choice; 
M_others.name   = 'Additional function';
M_others.tag    = 'M_others';
M_others.values = {E_markCardiac_TargetPCA, E_correctCardiac_TargetPCA, E_correctCardiac_exportBV,E_createonset_correlationsignal, E_eegnirs_MarkMuscular, E_GoNoGotrig}; %,E_createonset_correlationsignal
M_others.help   = {'These modules convert nir file in .nirs, last module of data export, support the export of field such as data (d),coordinate (SD), trigger (s), time (t),do not support noise artifact marking or aux export'};





f_HyperScan_outdir        = cfg_files;
f_HyperScan_outdir.name    = 'Select output folder'; 
f_HyperScan_outdir.tag     = 'f_HyperScan_outdir';       %file names
f_HyperScan_outdir.filter  = {'dir'};
f_HyperScan_outdir.ufilter = '.*';    
f_HyperScan_outdir.num     = [0 1];      % Number of inputs required 
f_HyperScan_outdir.val    = {''};
f_HyperScan_outdir.help    = {'Select the output folder where to save data. By default if you let empty it will save in the first NIRS.mat file selected'}; 

E_HyperScanCombineNIRS    = cfg_exbranch;
E_HyperScanCombineNIRS.name = 'HyperScan Combine';
E_HyperScanCombineNIRS.tag  = 'E_HyperScanCombineNIRS';
E_HyperScanCombineNIRS.val  = {NIRSmat, f_HyperScan_outdir};
E_HyperScanCombineNIRS.prog = @nirs_run_E_HyperScanCombineNIRS;
E_HyperScanCombineNIRS.vout = @nirs_cfg_vout_E_HyperScanCombineNIRS;
E_HyperScanCombineNIRS.help = {['Combine the data of multiple subjects (NIRS.mat) in one file,  detector, and source label will be kept identical for the first-subject and replaced as additional sources and detectors for the second subject; for example, second-subject detectors 1-16 will become detectors 17 to 32,(if the first subject have 16 detectors... ']};

function vout = nirs_cfg_vout_E_HyperScanCombineNIRS(job)
    vout = cfg_dep;                    
    vout.sname      = 'NIRS.mat';       
    vout.src_output = substruct('.','NIRSmat'); 
    vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end



% E_HyperScanDivideNIRS    = cfg_exbranch;
% E_HyperScanDivideNIRS.name = 'HyperScan Separate subject';
% E_HyperScanDivideNIRS.tag  = 'E_HyperScanDivideNIRS';
% E_HyperScanDivideNIRS.val  = {NIRSmat, f_HyperScan_outdir};
% E_HyperScanDivideNIRS.prog = @nirs_run_E_HyperScanDivideNIRS;
% E_HyperScanDivideNIRS.vout = @nirs_cfg_vout_E_HyperScanDivideNIRS;
% E_HyperScanDivideNIRS.help = {['UNCOMPLETE Divide the data of multiple subjects (NIRS.mat) in two file for hypyp compatibility ToDo']};
% 
% function vout = nirs_cfg_vout_E_HyperScanDivideNIRS(job)
%     vout = cfg_dep;                    
%     vout.sname      = 'NIRS.mat';       
%     vout.src_output = substruct('.','NIRSmat'); 
%     vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% end

%Module  Write external file
M_HyperScanNIRS        =  cfg_choice; 
M_HyperScanNIRS.name   = 'HyperScan';
M_HyperScanNIRS.tag    = 'M_HyperScanNIRS';
M_HyperScanNIRS.values = {E_HyperScanCombineNIRS };  %,  E_HyperScanDivideNIRS
M_HyperScanNIRS.help   = {'These modules is design to combine or divide 2 subjects during synchronized hyperscanning session. Combine the sessions to see to participant at once',...
    'After divide as subject for external software as Hypyp export (to complete).'};



%Module Utility
M_Utility        = cfg_choice; 
M_Utility.name   = 'Utility NIRSmat';
M_Utility.tag    = 'M_Utility';
M_Utility.values = {E_NIRSmatdiradjust, E_NIRSmatcreatenewbranch  E_createseedlist E_qualityreport E_zone2channellist E_channellist2zone E_VIDEO E_viewNIRS M_datawritenirs M_HyperScanNIRS  M_others }; %
M_Utility.help   = {'Utility on NIRSmat function.'};

nirsHSJ        = cfg_choice;
nirsHSJ.name   = 'LIONirs';
nirsHSJ.tag    = 'nirsHSJ'; %Careful, this tag nirsHSJ must be the same as
%the name of the toolbox and when called by spm_jobman in nirsHSJ.mM_dataprocessing 
nirsHSJ.values = {M_createMTG M_readNIRS M_Segment  M_preprocessing  M_dataComponent  M_Connectivity M_Utility E_GUI}; 

end
