%-----------------------------------------------------------------------
% Job saved on 26-May-2026 08:32:18 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF.inputSNIRF = {'C:\Data\ToolboxLIONirs\TutorialHandsOn2026_Visual\VISUAL.snirf'};
matlabbatch{1}.spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF.age1 = 25;
matlabbatch{1}.spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF.c_createImportProjectSnirf.b_createProject = struct([]);
matlabbatch{1}.spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF.output_path = 'C:\Data\ToolboxLIONirs\TutorialHandsOn2026_Visual\Analyse\C01';
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
matlabbatch{3}.spm.tools.nirsHSJ.M_dataComponent.E_extractcomponent.c_extractcomponent.b_extractcomponent_PCA.NIRSmat(1) = cfg_dep('Artifact Detection: NIRS.mat', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{3}.spm.tools.nirsHSJ.M_dataComponent.E_extractcomponent.c_extractcomponent.b_extractcomponent_PCA.m_extractcomponent = 'MVTPCA';
matlabbatch{3}.spm.tools.nirsHSJ.M_dataComponent.E_extractcomponent.c_extractcomponent.b_extractcomponent_PCA.m_extractcomponentfigure = 0;
matlabbatch{3}.spm.tools.nirsHSJ.M_dataComponent.E_extractcomponent.c_extractcomponent.b_extractcomponent_PCA.i_extract_pourcentagech = 5;
matlabbatch{3}.spm.tools.nirsHSJ.M_dataComponent.E_extractcomponent.c_extractcomponent.b_extractcomponent_PCA.i_extractnoiseupto_nbPCA = 97;
matlabbatch{3}.spm.tools.nirsHSJ.M_dataComponent.E_extractcomponent.c_extractcomponent.b_extractcomponent_PCA.i_extractnoise_nbPCA = 1;
matlabbatch{4}.spm.tools.nirsHSJ.M_dataComponent.E_substractcomponent.NIRSmat(1) = cfg_dep('Extract: NIRS.mat', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{4}.spm.tools.nirsHSJ.M_dataComponent.E_substractcomponent.i_substractcomponent_label = 'MVTPCA';
matlabbatch{4}.spm.tools.nirsHSJ.M_dataComponent.E_substractcomponent.m_substract_and_offsetcorrection = 1;
matlabbatch{5}.spm.tools.nirsHSJ.E_GUI.NIRSmat(1) = cfg_dep('Subtract: NIRS.mat', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','NIRSmat'));
matlabbatch{5}.spm.tools.nirsHSJ.E_GUI.DelPreviousData = 0;
