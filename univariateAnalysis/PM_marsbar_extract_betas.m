function [marsS] = PM_marsbar_extract_betas()

%roi_dir = '/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/fmri_data2/group_analyses/RetByAD_BalancedTrainingSet_NoLateralPFC_or_Parietal_FINAL18Subs_exMask_copy/CorConfXAD_vs_fix/ROIs/sphereROIs';
roi_dir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/group_analyses/percMnemConj_parModByCorConfAndRT_16Subs/ConfROI';

roi_names = {'PCC_roi.mat' 'rightAnG_roi.mat' 'leftAnG_roi.mat' 'leftTemporalLobe_roi.mat'};
%roi_names = {'RT_Left_AI_-30_26_4_roi.mat'	'RT_RightAI_33_26_-2_roi.mat' 'RT_Left_MFG_-48_20_34_roi.mat'};

analysisDir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/group_analyses/analysis_perc_conf_AcrossPerf_16Subs/';

d = dir([analysisDir '/*conf*']);


for i = 1:length(d)
    for curr_clust=1:length(roi_names) % go through the list of clusters
        
        spm_name = fullfile(analysisDir, d(i).name, 'SPM');
        
        thisSPM = load(spm_name);
        
        D = mardo(spm_name);
        D = autocorr(D,'fmristat',2);
        
        clusters_h =  [roi_dir '/' roi_names{curr_clust}];
        
        % Make marsbar design object
        roi_file = clusters_h;
        
        % Make marsbar ROI object
        R  = maroi(roi_file);
        
        % Fetch data into marsbar data object
        Y  = get_marsy(R, D, 'mean');
        
        % Get contrasts from original design
        xCon = get_contrasts(D);
        
        % Estimate design on ROI data
        E = estimate(D, Y);
        
        % Put contrasts from original design back into design object
        E = set_contrasts(E, xCon);
        
        % get design betas
        b = betas(E);
        
        yS = y_struct(Y);
        
        marsS.clust(curr_clust).cond(i).raw = yS.Y;
        
        marsS.clust(curr_clust).cond(i).betas = b;
        
        marsS.clust(curr_clust).cond(i).names = {d.name};
        
        % get stats and stuff for all contrasts into statistics structure
        marsS.clust(curr_clust).cond(i).stats = compute_contrasts(E, 1:length(xCon));
    end
end


% for i = 1:length(marsS.clust)
%     for j = 1:length(marsS.clust(i).sub);
%         for k = 1:length(marsS.clust(i).sub(j).con)
%             q(j,i,k) = marsS.clust(i).sub(j).con(k);
%         end
%     end
% end

            

