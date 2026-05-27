%-----------------------------------------------------------------------
% Job saved on 26-May-2026 08:52:30 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatdiradjust.NIRSmat = {'C:\Data\ToolboxLIONirs\TutorialHandsOn2026_Visual\Analyse\C01\dCONC\NIRS.mat'};
matlabbatch{1}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatdiradjust.c_MultimodalPath.b_MultimodalPath_no = struct([]);
matlabbatch{2}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.NIRSmat(1) = cfg_dep('Folder adjustment: NIRS.mat', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{2}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.e_NIRSmatdirnewbranch = 'AvgBL';
matlabbatch{2}.spm.tools.nirsHSJ.M_Utility.E_NIRSmatcreatenewbranch.m_newbranchcomponent = 1;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.NIRSmat(1) = cfg_dep('New branch: NIRS.mat', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.DelPreviousData = 0;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.trigger = 1;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.pretime = '5';
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.posttime = '45';
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.badintervalratio = 0.5;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.badchannelratio = 0.5;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.avg_datatype = 2;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.c_baseline_corr.m_defaultbaseline_corr = 1;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.c_rejecttrial.m_noreject_trial = 1;
matlabbatch{3}.spm.tools.nirsHSJ.M_preprocessing.E_average.choiceave.m_Tvalueoption = 2;
matlabbatch{4}.spm.tools.nirsHSJ.E_GUI.NIRSmat(1) = cfg_dep('Epoch averaging: NIRS.mat', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{4}.spm.tools.nirsHSJ.E_GUI.DelPreviousData = 0;
