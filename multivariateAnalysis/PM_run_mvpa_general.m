function [res, results]= PM_run_mvpa_general(subj_array, task, saveName, portion)


for b=(1:length(subj_array))
    
    %% load general parameter information
    tic;
    [S idxTr idxTe par] = PM_mvpa_params(subj_array(b), task);
    if nargin > 3
        S.portion = portion;
    else
        S.portion = [];
    end
    S.idxTr = idxTr;
    S.idxTe = idxTe;
    S.saveName = saveName; %name of classifier results
    S.subj_array = subj_array; %subjects
    
    %% information about which TRs to include in classification
    %which weighted combination of post-stimulus TRs should be used to train the classifier?
    S.TR_weights_train = S.TR_weights_set{1}; % should sum to 1
    S.TRs_to_average_over_train = 1:length(S.TR_weights_train);
    
    %which weighted combination of post-stimulus TRs should be used to test the classifier?
    S.TR_weights_test = S.TR_weights_set{2}; % should sum to 1
    S.TRs_to_average_over_test = 1:length(S.TR_weights_test);
    
    S.TR_weights = S.TR_weights_set;
    S.TRs_to_average_over = 1:length(S.TR_weights);
    
    %% Onsets
    S = PM_mvpa_onsets_and_images(S);
    S.num_conds = size(S.onsets,2);
    
    %% subsample to match RTs
    if S.subsampleToMatch
        [S] = PM_subsample(S);
    end
    
    %% Workspace stuff
    existWorkspace = exist(S.workspace);
    
    % load workspace
    if (S.use_premade_workspace&&existWorkspace) % if we are supposed to use premade workspace, and one with the correct name exists
        load(S.workspace, 'subj');
    else
        [subj] = PM_mvpa_load_and_preprocess_raw_data(S);
    end
    
    %% mask a workspace mask with another mask.
    if ~isempty(S.secondaryMask)
        subj = load_spm_mask(subj, 'secondaryMask', S.secondaryMask);
        subj = intersect_masks(subj,S.roi_name,'secondaryMask');
        subj = create_pattern_from_mask(subj, S.preprocPatName, subj.masks{end}.name , [S.preprocPatName '_masked']);
    end
    
    %% begin classifier loops
    for n = 1: S.num_results_iter
        
        if strcmp(S.patternType, 'raw')
            all_regs = zeros(S.num_conds,S.num_vols); % initialize regs matrix as conditions x timepoints
            
            % convert from seconds to TRs
            for cond = 1:S.num_conds
                for trial = 1: length(S.onsets{cond})
                    time_idx = round(S.onsets{cond}(trial)/S.TR) + 1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
                    all_regs(cond, round(time_idx)) = 1;
                end
            end
            
            % condense regs by removing zeros
            condensed_runs = [];
            condensed_regs_of_interest = [];
            trial_counter = 1;
            for i = 1: size(all_regs,2)
                if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
                    condensed_regs_of_interest(:,trial_counter) = all_regs(:,i);
                    condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
                    trial_counter = trial_counter + 1;
                end
            end
            idx_condense =find(sum(all_regs));
            
            % condense meta_runs using idx_condense
            trial_idx = 1;
            m_runs = 0;
            for r = 1:length(S.meta_runs)
                m_runs(trial_idx:trial_idx+S.meta_runs(r)-1)=r;
                trial_idx = trial_idx+S.meta_runs(r);
            end
            meta_runs_condensed = m_runs(idx_condense);
            
            
            %% select active trials
            S.actives = ones(size(meta_runs_condensed));
            
            % remove artifactual trials in the analysis
            if S.inactivateArtifacts
                [subj S] = PM_inactivateArtifacts(subj, S);
            end
            
            % index active training trials
            allTrainOns = sort([S.onsets_train_in_classifier{:}]);
            allOns = sort([S.onsets{:}]);
            S.trainActives = ismember(allOns, allTrainOns);
            subj = init_object(subj,'selector','trainActives');
            subj = set_mat(subj,'selector','trainActives', S.trainActives);
            
            % select active trials for linear regression
            if S.linReg
                [subj S] = PM_linRegSetup(subj, S);
            end
            
            subj = init_object(subj,'selector','actives');
            subj = set_mat(subj,'selector','actives', S.actives);
            
            
            %% load data
            
            % create meta_runs structure
            all_trials = sum(all_regs,1);
            if strcmp(S.trainTask, S.testTask)
                meta_runs_train = find(all_trials);
                meta_runs_test = [];
                randomNFold = ceil(shuffle(1:length(meta_runs_condensed))/(length(meta_runs_condensed)/S.nFolds));
                subj = init_object(subj,'selector','randomNFold');
                subj = set_mat(subj,'selector','randomNFold', randomNFold);
                subj = set_objfield(subj, 'selector', 'randomNFold', 'group_name', 'randomNFoldGroup');
                subj = PM_create_xvalid_indices_trainActivesOnly(subj,'randomNFold', 'actives_selname', 'trainActives');
            else
                %applies only when train and test data are different 
                TrainTestOneIter = 1*ismember(meta_runs_condensed, 1:length(S.TrainRuns)) + 2*ismember(meta_runs_condensed, length(S.TrainRuns)+1: length(S.runs_vector));
                TrainTestOneIter(TrainTestOneIter==1) = S.trainActives(TrainTestOneIter==1);
                meta_runs_train = idx_condense(find(TrainTestOneIter==1));
                meta_runs_test = idx_condense(find(TrainTestOneIter==2));
                subj = init_object(subj,'selector','TrainTestOneIter');
                subj = set_mat(subj,'selector','TrainTestOneIter', TrainTestOneIter);
                subj = set_objfield(subj, 'selector', 'TrainTestOneIter', 'group_name', 'TrainTestOneIterGroup');
            end
            
            % load training data
            data_by_TR_train = [];
            for dt = 1:length(S.TR_weights_set{1})
                data_by_TR_train(dt,:,:) = S.TR_weights_train(dt)*subj.patterns{end}.mat(:,meta_runs_train+(dt-1));
            end
            temporally_condensed_data_train = squeeze(sum(data_by_TR_train(S.TRs_to_average_over_train,:,:),1));
            clear data_by_TR_train
            
            % load testing data
            data_by_TR_test = [];
            for dt = 1:length(S.TR_weights_set{2})
                data_by_TR_test(dt,:,:) = S.TR_weights_test(dt)*subj.patterns{end}.mat(:,meta_runs_test+(dt-1));
            end
            temporally_condensed_data_test = squeeze(sum(data_by_TR_test(S.TRs_to_average_over_test,:,:),1));
            clear data_by_TR_test
            
            % combine training and testing data
            temporally_condensed_data = horzcat(temporally_condensed_data_train, temporally_condensed_data_test);
            clear temporally_condensed_data_train;
            clear temporally_condensed_data_test;
            
            %% Important note
            % train patterns and onsets are always first, followed
            % by test patterns and onsets.
            
            %% modify subj object
%             if strcmp(S.thisSelector, 'TrainTestOneIterGroup')
%                 %create TrainTestOneIter selector
%                 TrainTestOneIter = 1*ismember(meta_runs_condensed, 1:length(S.TrainRuns)) + 2*ismember(meta_runs_condensed, length(S.TrainRuns)+1: length(S.runs_vector));
%                 TrainTestOneIter(TrainTestOneIter==1) = S.trainActives(TrainTestOneIter==1);
%                 subj = init_object(subj,'selector','TrainTestOneIter');
%                 subj = set_mat(subj,'selector','TrainTestOneIter', TrainTestOneIter);
%                 subj = set_objfield(subj, 'selector', 'TrainTestOneIter', 'group_name', 'TrainTestOneIterGroup');
%             elseif strcmp(S.thisSelector, 'randomNFold_xval')
%                 %create randomNFold selector
%                 randomNFold = ceil(shuffle(1:length(meta_runs_condensed))/(length(meta_runs_condensed)/S.nFolds));
%                 subj = init_object(subj,'selector','randomNFold');
%                 subj = set_mat(subj,'selector','randomNFold', randomNFold);
%                 subj = set_objfield(subj, 'selector', 'randomNFold', 'group_name', 'randomNFoldGroup');
%                 subj = PM_create_xvalid_indices_trainActivesOnly(subj,'randomNFold', 'actives_selname', 'trainActives');
%             end
            
            %create meta_runs_condensed selector
            subj = init_object(subj,'selector','meta_runs_condensed');
            subj = set_mat(subj,'selector','meta_runs_condensed', meta_runs_condensed);
            subj = create_xvalid_indices(subj,'meta_runs_condensed');
            
            % only include 'active' patterns in selector
            grp = find_group(subj, 'selector', S.thisSelector);
            for g = 1:length(grp)
                this_mat = get_mat(subj,'selector',grp{g});
                this_mat(this_mat==1) = this_mat(this_mat==1) .* S.actives(this_mat==1);
                subj = set_mat(subj,'selector',grp{g},this_mat);
            end
            
            % add conditions
            subj = init_object(subj,'regressors','conds');
            subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
            subj = set_objfield(subj,'regressors','conds','condnames',S.condnames);
            
            % add new condensed activation pattern
            subj = duplicate_object(subj,'pattern',S.preprocPatNameFinalMask,S.preprocPatCondensedName);
            subj = set_mat(subj,'pattern',S.preprocPatCondensedName,temporally_condensed_data,'ignore_diff_size',true);
            zhist = sprintf('Pattern ''%s'' created by AG custom code',S.preprocPatCondensedName);
            subj = add_history(subj,'pattern',S.preprocPatCondensedName,zhist,true);
            
            % clean up workspace to save RAM
            subj = remove_mat(subj,'pattern',S.preprocPatNameFinalMask);
            subj = remove_mat(subj,'pattern',S.preprocPatName);
            subj.selectors{1}.mat = condensed_runs;
            subj.selectors{1}.matsize = size(condensed_runs);
            
            S.classSelector = S.thisSelector;
            
        elseif strcmp(S.patternType, 'betas')
            [subj S] = PM_organizeBetasForClassification;
        end
        
        %equate the training set.
        if S.equate_number_of_trials_in_groups
            subj = PM_balanceTrainPats(S, subj);
            S.classSelector = [S.thisSelector 'balanced'];
        end

        S.classifier_pattern = S.preprocPatCondensedName; % data to use for classification.
        S.classifier_mask = subj.masks{end}.name; % mask to use for classification.
        
        %zscore the patterns prior to classification
        if S.perform_second_round_of_zscoring
            display('Performing second round of z-scoring')
            subj = zscore_runs(subj,S.preprocPatCondensedName,'runs'); % Z-score the data
            S.classifier_pattern = [S.preprocPatCondensedName '_z']; % update the classifier data of interest
        end
        
        % run feature selection ANOVA: specify #of voxels (if desired)
        if S.class_args.nVox>0
            display('Performing feature selection')
            statmap_arg = [];
            subj = AG_feature_select_top_N_vox(subj,S.preprocPatCondensedName,'conds',S.classSelector,'nVox_thresh',S.class_args.nVox, 'statmap_funct', S.statmap_funct, 'statmap_arg',statmap_arg);
            S.classifier_mask = subj.masks{end}.name; % use group of masks created by ANOVA
            S.classifier_mask_group = subj.masks{end}.group_name;
        end
        
         %include interactions between voxels
        if S.includeVoxelInteractions
            [subj S] = PM_model_voxel_interactions(subj, S);
        end
        %%
        
        %S.class_args.penalty = S.penaltyParams(pnl);
        
        if S.extractMeanSignal
            % instead of doing classification, just extract
            % the mean signal in a given set of regions and
            % compare them
            [subj results] =  AG_extractMeanSignal(subj,S);
        else
            %run the classification.
            [subj results] = cross_validation(subj,S.classifier_pattern,'conds', ...
            S.classSelector, S.classifier_mask,S.class_args, 'perfmet_functs', S.perfmet_functs);
        end
        
        %set up importance maps.
        if S.generate_importance_maps == 1
            for rif = 1:length(results.iterations);
                thisScratch = results.iterations(rif).scratchpad.w(2:end,:)';
                results_IW{rif}.iterations(1).scratchpad.net.IW{1} = thisScratch;
            end
        end
        
        %store results
        res.subj{b}.penalty(1).nVox(1).weights(1).iter{n} = results;
        res.subj{b}.penalty(1).nVox(1).weights(1).S = S;
        res.subjArray = subj_array;
        
        %save results
        if ~(exist(S.group_mvpa_dir))
            mkdir(S.group_mvpa_dir);
        end
        save (fullfile(S.group_mvpa_dir, S.saveName), 'res');
        
        % display time it took.
        time2finish = toc/60;
        display(['Finished ' S.subj_id ' in ' num2str(time2finish) ' minutes']);
        
        % generate importance maps.
        if S.generate_importance_maps
            PM_generate_importance_maps(subj, results, results_IW, S)
        end
        
        clear subj
    end
end


















