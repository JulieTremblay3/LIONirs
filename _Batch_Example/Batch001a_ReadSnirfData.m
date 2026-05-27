%-----------------------------------------------------------------------
% Job saved on 11-May-2026 11:28:11 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF.inputSNIRF = '<UNDEFINED>';
matlabbatch{1}.spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF.age1 = 25;
matlabbatch{1}.spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF.c_createImportProjectSnirf.b_createProject = struct([]);
matlabbatch{1}.spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF.output_path = 'c:/Data/Analyzed/Csub01_tapping';
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.NIRSmat(1) = cfg_dep('Read .snirf: NIRS.mat', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.DelPreviousData = 0;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.PrintReport = 1;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_meandiff.m_meandiff = 1;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_meandiff.thresholdstep = 3;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_meandiff.m_thresholdstepzscore = 1;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_meandiff.movingaverage_nbpoint = 1;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_meandiff.too_small_step_dur = 0.4;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_meandiff.printreportthreshold = 0;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_minpourcentagebad.m_minpourcentagebad = 0;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_minpourcentagebad.minpourcentagebad = 3;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_min_subinterval.m_min_subinterval = 1;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_min_subinterval.min_subinterval = 2;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_corr.m_corr = 1;
matlabbatch{2}.spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection.b_corr.corr_thr = 0.8;
matlabbatch{3}.spm.tools.nirsHSJ.E_GUI.NIRSmat(1) = cfg_dep('Artifact Detection: NIRS.mat', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{3}.spm.tools.nirsHSJ.E_GUI.DelPreviousData = 0;
