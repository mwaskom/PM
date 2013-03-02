function subjects_tc = PM_marsbar_batch()

[sa, SA ] = PM_SA;

%conds = {'face_cor_highconf' 'house_cor_highconf' 'face_cor_lowconf' 'house_cor_lowconf'};

% conds = {'face_cor_conf1' 'house_cor_conf1' 'face_cor_conf2' 'house_cor_conf2' ...
%      'face_cor_conf3' 'house_cor_conf3' 'face_cor_conf4' 'house_cor_conf4' };

conds = {'face_cor_RT1' 'house_cor_RT1' 'face_cor_RT2'  'house_cor_RT2' ...
    'face_cor_RT3' 'house_cor_RT3' 'face_cor_RT4' 'house_cor_RT4'};
     
%  conds = {'face_RT1_cor' 'house_RT1_cor' 'face_RT2_cor' 'house_RT2_cor' ...
%       'face_RT3' 'house_RT3' 'face_RT4' 'house_RT4' };

% conds = {'face_conf1_corWithZeros' 'house_conf1_corWithZeros' 'face_conf2_corWithZeros' ... 
%     'house_conf2_corWithZeros' 'face_conf3_corWithZeros' ...
%     'house_conf3_corWithZeros' 'face_conf4_corWithZeros' 'house_conf4_corWithZeros'};


%roi_dir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/group_analyses/percMnemConj_parModByConfAndRT_16Subs/ROIs';
%roi_dir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/group_analyses/percMnemConj_parModByCorConfAndRT_16Subs/RT_ROI';

%roi_names = {'rt_p005_leftAI.nii'	'rt_p005_leftMFG.nii'	'rt_p005_rightAI.nii'};
%roi_names = {'conf_p005_left_AnG.nii'	'conf_p005_left_STG.nii' 'conf_p005_PCC.nii'	};
roi_names = {'RT_p005_leftAI.nii'	'RT_p005_leftMFG.nii'	'RT_p005_rightAI.nii'};
%roi_names = {'conf_p005_LeftTemporalLobe.nii' 'conf_p005_leftAnG.nii' 'conf_p005_PCC.nii'	'conf_p005_rightAnG.nii'};
%roi_names = {'AnG_L_-39_-70_40_roi.mat' 'PCC_-6_-52_31_roi.mat' 'STG_L_-63_-46_1_roi.mat'};
%roi_names = {'RT_LeftAI_-30_23_4_roi.mat' 'RT_RightAI_-30_23_4_roi.mat' 'RT_LeftMFG_-42_11_22_roi.mat'};
%roi_names = {'PCC_roi.mat' 'rightAnG_roi.mat' 'leftAnG_roi.mat' 'leftTemporalLobe_roi.mat'};
%roi_names = {'RT_Left_AI_-30_26_4_roi.mat'	'RT_RightAI_33_26_-2_roi.mat' 'RT_Left_MFG_-48_20_34_roi.mat'};

subjects = sa.sa16_Conf;
 
% for r = 1:length(roi_names)
%     clusters_h{r} = [roi_dir '/' roi_names{r}];
% end
% 
% clusters = char(clusters_h{:});

  
for curr_subj=1:length(subjects), % go through the list of subjects
    
    parM = PM_Params(subjects(curr_subj), 'mnem', 0);
    parP = PM_Params(subjects(curr_subj), 'perc', 0);
    roi_dir = fullfile(parM.mnem.analysisdir, 'percMnemConj_parModByConfAndRT_16Subs_1LeftOut');
   

    fprintf(1,'extracting TC from subject %g\n',subjects(curr_subj));
    spm_name = [parM.analysisdir '/SPM.mat'];
    D = mardo(spm_name);
    D = autocorr(D,'fmristat',2);
    for curr_clust=1:length(roi_names) % go through the list of clusters
        
        clusters_h =  [roi_dir '/' roi_names{curr_clust}];
        
        % Make marsbar design object
        roi_file = clusters_h;
        
        % Make marsbar ROI object
        R  = maroi_image(roi_file);
        
        % Fetch data into marsbar data object
        Y  = get_marsy(R, D, 'mean');
        
        % Get contrasts from original design
        xCon = get_contrasts(D);
        
        % Estimate design on ROI data
        E = estimate(D, Y);
        
        % Put contrasts from original design back into design object
        % E = set_contrasts(E, xCon);
        
        % get design betas
        %b = betas(E);
        
        % get stats and stuff for all contrasts into statistics structure
        %marsS = compute_contrasts(E, 1:length(xCon));
        
        % Get definitions of all events in model
        [e_specs, e_names] = event_specs(E);
        n_events = size(e_specs, 2);
        
        % Bin size in seconds for FIR
        % bin_size = tr(E); %
        bin_size = 2;
        
        % Length of FIR in seconds
        fir_length = 20;
        
        % Number of FIR time bins to cover length of FIR
        bin_no = fir_length / bin_size;
        
        % Options - here 'single' FIR model, return estimated % signal change
        opts = struct('single', 1, 'percent', 1);
        
        % Return time courses for all events in fir_tc matrix
        
        %      for e_s = 1:n_events
        %        fir_tc(:, e_s) = event_fitted_fir(E, e_specs(:,e_s), bin_size, ...
        %                                          bin_no, opts);
        %      end
        
        %for cds = 1:length(condType)
            %subjects_tc(curr_subj).rois(curr_clust).conds = condType{cds};
            %subjects_tc(curr_subj).rois(curr_clust).condNames = conds{cds};
            clear fir_tc;
            for c = 1:length(conds)
                e_s = find(strcmp(conds{c}, e_names));
                
                if isempty(e_s)
                    fir_tc(:, c) = nan(fir_length/bin_size,1);
                    subjects_tc(curr_subj).rois(curr_clust).maxWithNeighbors(c) = NaN;
                else
                    
                    fir_tc(:, c) = event_fitted_fir(E, e_specs(:,e_s), bin_size, bin_no, opts);
                    
                    [maxC, maxI] = max(fir_tc(2:((fir_length/bin_size)-1),c));
                    
                    
                    %subjects_tc(curr_subj).rois(curr_clust).maxval(c) = maxC;
                    
                    subjects_tc(curr_subj).rois(curr_clust).maxval(c) = maxC;
                    neighbors = intersect(1:size(fir_tc,2), [maxI-1: maxI+1]);
                    subjects_tc(curr_subj).rois(curr_clust).maxWithNeighbors(c) = mean(fir_tc(neighbors,c));
                    
                end
                subjects_tc(curr_subj).rois(curr_clust).conds = conds;
                
                subjects_tc(curr_subj).rois(curr_clust).roiName = roi_names{curr_clust};
                subjects_tc(curr_subj).rois(curr_clust).raw = fir_tc;
                
            end
            
            save(fullfile('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/group_analyses/percMnemConj_parModByConfAndRT_16Subs/ROIs_RT', 'subjects_tc_mnemConfIndependent.mat'), 'subjects_tc');
    %end
        
        
    end % end clusters loop
    
    subjects_tc(curr_subj).conds = conds;
    subjects_tc(curr_subj).e_specs = e_specs;
    subjects_tc(curr_subj).e_names = e_names;
    subjects_tc(curr_subj).n_events = n_events;
    
    
end % end subjects loop

% end function


%% Scrap post-processing stuff
%[M reh groupPsy idxG] = MnemonicBehAnalysisWholeshebang(sa.sa16_Conf);

%idxEnoughTrials = ~M.ret.vals(:,8:15)>9;

cols = cool(4);

clear w_roi meanAct steAct

for j=1:length(roi_names)
    for i=1:length(subjects_tc)
        w_roi{j}.mat(i,:,:) = subjects_tc(i).rois(j).raw;
        %w_roi{j}.mat(i,:,idxEnoughTrials(i,:)) = NaN; %not large enough to count.
    end
    
    for k = 1:4
    %for k=1:2
        theseCols = (k*2-1):(k*2);
%         if k==1
%             theseCols = 7:8;
%         else
%             theseCols = 1:6;
%         end
        w_roi{j}.matAcrossClass(:,:,k) = nanmean(w_roi{j}.mat(:,:,theseCols),3);
    end
    
    meanAct = squeeze(nanmean(w_roi{j}.matAcrossClass,1));
    steAct = squeeze(nanstd(w_roi{j}.matAcrossClass,1))/sqrt(14);
    
    figure('position',[200 200 200 200]); hold on
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [.5 .5 2 2]); % last 2 are width/height.
    
    for l=1:4
        errorbar(meanAct(:,l), steAct(:,l), 'Color', cols(l,:), 'LineWidth', 1.5)
        a_h = gca;
        set(a_h,'FontSize',12)
        xlim([0 11])
        set(gca,'XTick',[1:2:11'])
        set(gca,'XTickLabel',[1:4:21])
    end
    
    f_h = gcf;
    cd('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/group_analyses/percMnemConj_parModByConfAndRT_16Subs/ROIs_RT');
    print(f_h, '-depsc', [subjects_tc(1).rois(j).roiName(1:(end-4)), '_mnemRTIndependentROIDef.eps'] )
    %errorbar(meanAct, steAct)
end




                                    
                                   
