%-----------------------------------------------------------------------
% Job saved on 26-May-2026 08:55:58 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatdiradjust.NIRSmat = {'C:\Data\ToolboxLIONirs\TutorialHandsOn2026_Visual\Analyse\C01\NIRS.mat'};
matlabbatch{1}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatdiradjust.c_MultimodalPath.b_MultimodalPath_no = struct([]);
matlabbatch{2}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.NIRSmat(1) = cfg_dep('Folder adjustment: NIRS.mat', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{2}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.e_NIRSmatdirnewbranch = 'dCONC';
matlabbatch{2}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.m_newbranchcomponent = 0;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_chcardiaccontrol.NIRSmat(1) = cfg_dep('New branch: NIRS.mat', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_chcardiaccontrol.i_Freq_cardiac = [0.8 2.3];
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_chcardiaccontrol.i_Freq_crossspectrum = 0.1;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_chcardiaccontrol.i_minch_cardiac = 10;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_chcardiaccontrol.i_cardiacwidth = 0;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_chcardiaccontrol.m_cardiacwavelenght = 3;
matlabbatch{4}.spm.tools.nirsHSJ.M_preprocessing.normalization.NIRSmat(1) = cfg_dep('Cardiac Detection: NIRS.mat', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{4}.spm.tools.nirsHSJ.M_preprocessing.normalization.DelPreviousData = 0;
matlabbatch{4}.spm.tools.nirsHSJ.M_preprocessing.normalization.normtype.b_choicenormstim.trigger = [1 2 3 4];
matlabbatch{4}.spm.tools.nirsHSJ.M_preprocessing.normalization.normtype.b_choicenormstim.pretime = '5';
matlabbatch{4}.spm.tools.nirsHSJ.M_preprocessing.normalization.normtype.b_choicenormstim.posttime = '45';
matlabbatch{4}.spm.tools.nirsHSJ.M_preprocessing.normalization.normtype.b_choicenormstim.m_NormType = 0;
matlabbatch{4}.spm.tools.nirsHSJ.M_preprocessing.normalization.normtype.b_choicenormstim.m_choiceNan = 0;
matlabbatch{5}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.NIRSmat(1) = cfg_dep('Normalization dOD = log(I/Io): NIRS.mat', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{5}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.DelPreviousData = 1;
matlabbatch{5}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.lowcutfreq = '0.2';
matlabbatch{5}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.highcutfreq = 'No';
matlabbatch{5}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.filterorder = 4;
matlabbatch{5}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.paddingsymfilter = 1;
matlabbatch{5}.spm.tools.nirsHSJ.M_preprocessing.bpfilt.interpolatebadfilter = 0;
matlabbatch{6}.spm.tools.nirsHSJ.M_preprocessing.ODtoHbOHbR.NIRSmat(1) = cfg_dep('Filter: NIRS.mat', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{6}.spm.tools.nirsHSJ.M_preprocessing.ODtoHbOHbR.DelPreviousData = 1;
matlabbatch{6}.spm.tools.nirsHSJ.M_preprocessing.ODtoHbOHbR.PVF = [1 1];
matlabbatch{6}.spm.tools.nirsHSJ.M_preprocessing.ODtoHbOHbR.C_ODtoHbOHbR_DPF.b_ODtoHbOHbR_DPF1 = struct([]);
matlabbatch{7}.spm.tools.nirsHSJ.E_GUI.NIRSmat(1) = cfg_dep('Modified Beer Lambert Law: NIRS.mat', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{7}.spm.tools.nirsHSJ.E_GUI.DelPreviousData = 0;
