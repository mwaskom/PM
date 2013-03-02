function [res, results]= PM_run_mvpa_general2(subj_array, task, saveName) %#ok<INUSL>



%% subject loop
for b=(1:length(subj_array))
    
    
    %[S idxB par] = PM_mvpa_params(subj_array(b), task);
    S = AG_mvpa_params(subj_array{b});
    
    % pattern weights loop
    for TRW = 1:length(S.TR_weights_set)
        
        tic;
        
        S.TR_weights = S.TR_weights_set{TRW};
        S.TRs_to_average_over = 1:length(S.TR_weights);
        
        S.saveName= saveName; %name of classifier results
        S.subj_array = subj_array; %subjects
        
        %% Onsets
        [S.onsets S.RTs S.enoughTRs S.ixOnsets] = PM_mvpa_onsets_and_images(S);
        S.num_conds = size(S.onsets,2);
        
        %% Workspace stuff
        existWorkspace = exist(S.workspace, 'file');
        
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
        end
        
        
        % if you want to run the classification multiple times, with
        % different initial weights seeding or different random
        % balanced subsets of data
        %% results loop
        for n = 1:S.num_results_iter
            
            % initialize regs matrix as conditions x timepoints
            all_regs = zeros(S.num_conds,S.num_vols);
            
            % convert from seconds to TRs
            for cond = 1: S.num_conds
                for trial = 1: length(S.onsets{cond})
                    % divide by TR length and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
                    time_idx = round(S.onsets{cond}(trial)/S.TR) + 1;
                    all_regs(cond, round(time_idx)) = 1;
                end
                
            end
            
            % condense regs by removing zeros
            condensed_runs =  nan(size(all_regs,1), sum(sum(all_regs)));
            condensed_regs_of_interest = nan(size(all_regs,1), sum(sum(all_regs)));
            trial_counter = 1;
            for i = 1: size(all_regs,2)
                if ~isempty(find(all_regs(:,i),1)) % if not a rest timepoint
                    condensed_regs_of_interest(:,trial_counter) = all_regs(:,i);
                    condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
                    trial_counter = trial_counter + 1;
                end
            end
            idx_condense = logical(sum(all_regs));
            
            % meta runs codes is the label for scan sessions.
            trial_idx = 1;
            m_runs = 0;
            for r = 1:length(S.meta_runs)
                m_runs(trial_idx:trial_idx+S.meta_runs(r)-1)=r;
                trial_idx = trial_idx+S.meta_runs(r);
            end
            meta_runs_condensed = m_runs(idx_condense);
            
            %Train on S.TrainRuns, test on S.TestRuns
            TrainTestOneIter = 1*ismember(meta_runs_condensed, 1:length(S.TrainRuns)) ...
                + 2*ismember(meta_runs_condensed, length(S.TrainRuns)+1: length(S.runs_vector));
            
            
            if S.linReg
                %If we are running a linear regression to predict RT, deactivate the
                %incorrect patterns.
                if strcmp(S.thisSelector, 'TrainTestOneIter')
                    actives = ones(size(TrainTestOneIter));
                    actives(TrainTestOneIter==1) = idxB.corWithZeros';
                    
                    condensed_regs_of_interest = zeros(size(TrainTestOneIter));
                    condensed_regs_of_interest(TrainTestOneIter==1) = idxB.coh_signed';
                else
                    actives = idxB.corWithZeros';
                    condensed_regs_of_interest = idxB.coh_signed';
                end
            else
                %otherwise, all patterns are active
                actives_h = ones(size(meta_runs_condensed));
                actives = actives_h;
            end
            

            
            all_trials = sum(all_regs,1);            
            meta_runs = find(all_trials);
            
            
            %all data, by TR
            %figure out how to initialize this...
            %data_by_TR = nan(length(S.TR_weights), size(S.TR_weights(dt),1), );
            for dt = 1:length(S.TR_weights)
                data_by_TR(dt,:,:) = S.TR_weights(dt)*subj.patterns{end}.mat(:,meta_runs+(dt-1));
            end
            
            %condense the data, maintaining only the TRs we care about
            temporally_condensed_data = squeeze(sum(data_by_TR(S.TRs_to_average_over,:,:),1));
            
            clear data_by_TR
            
            
            %% modify subj object
            
            % store actives variable in subj object
            subj = init_object(subj,'selector','actives');
            subj = set_mat(subj,'selector','actives', actives);
            
            %create TrainTestOneIter selector
            subj = init_object(subj,'selector','TrainTestOneIter');
            subj = set_mat(subj,'selector','TrainTestOneIter', TrainTestOneIter);
            subj = set_objfield(subj, 'selector', 'TrainTestOneIter', 'group_name', 'TrainTestOneIterGroup');
            
            %create meta_runs_condensed selector
            subj = init_object(subj,'selector','meta_runs_condensed');
            subj = set_mat(subj,'selector','meta_runs_condensed', meta_runs_condensed);
            subj = create_xvalid_indices(subj,'meta_runs_condensed');
            
            %create randomNFold selector
            
            if S.nFolds==1
                randomNFold = shuffle(1:length(meta_runs_condensed)); 
            else
                randomNFold_h = shuffle(1:length(meta_runs_condensed));
                maxBinSize = length(randomNFold_h)./S.nFolds;
                randomNFold = ceil(randomNFold_h/maxBinSize)';
            end
            
            subj = init_object(subj,'selector','randomNFold');
            subj = set_mat(subj,'selector','randomNFold', randomNFold);
            subj = create_xvalid_indices(subj,'randomNFold');
            
            grp = find_group(subj, 'selector', S.thisSelector);
            for g = 1:length(grp)
                this_mat = get_mat(subj,'selector',grp{g});
                
                % remove non-train items from train labels
                this_mat(this_mat==1) = this_mat(this_mat==1) .* actives(this_mat==1) .* S.ixOnsets.inTrainSetConcat(this_mat==1)';
               
                % remove non-test items from test labels
                this_mat(this_mat==2) = this_mat(this_mat==2) .* actives(this_mat==2) .* S.ixOnsets.inTestSetConcat(this_mat==2)';
                
                subj = set_mat(subj,'selector',grp{g},this_mat);
            end
            
  
            % initialize regressors object
            subj = init_object(subj,'regressors','conds');
            subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
            subj = set_objfield(subj,'regressors','conds','condnames',S.condnames);
            
            % add new condensed activation pattern
            subj = duplicate_object(subj,'pattern',S.preprocPatName,S.preprocPatCondensedName);
            subj = set_mat(subj,'pattern',S.preprocPatCondensedName,temporally_condensed_data,'ignore_diff_size',true);
            
            zhist = sprintf('Pattern ''%s'' created by AG custom code',S.preprocPatCondensedName);
            subj = add_history(subj,'pattern',S.preprocPatCondensedName,zhist,true);
            
            % clean up workspace to save RAM
            subj = remove_mat(subj,'pattern',S.preprocPatName);
            
            %create condensed_runs selector
            subj.selectors{1}.mat = condensed_runs;
            subj.selectors{1}.matsize = size(condensed_runs);
            
            S.classSelector = S.thisSelector;
            
            %balance the training set.
            if S.equate_number_of_trials_in_cond_1_and_2
                subj = PM_balanceTrainPats(S, subj);
                S.classSelector = [S.thisSelector 'balanced'];
            end
            
            %re-zscore the patterns            
            if S.perform_second_round_of_zscoring
                display('Performing second round of z-scoring')
                subj = zscore_runs(subj,'spiral_d_z_condensed','meta_runs_condensed');
                %subj.patterns{end}.mat = zscore(subj.patterns{end}.mat);
                %subj.patterns{end}.mat(:,active_trials) = zscore(subj.patterns{end}.mat(:,active_trials));
                S.preprocPatCondensedName = 'spiral_d_z_condensed_z';
            end
            
            statmap_arg = [];
            classifier_mask = subj.masks{end}.name; %
            
            % run feature selection ANOVA: specify pvalue (if desired)
            if S.class_args.pThresh ~= 1
                subj = feature_select(subj,S.preprocPatCondensedName,'conds',S.classSelector,'thresh',S.class_args.pThresh,'statmap_arg',statmap_arg);
                classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
            end
            
            % run feature selection ANOVA: specify #of voxels (if desired)
            if S.class_args.nVox>0
                subj = AG_feature_select_top_N_vox(subj,S.preprocPatCondensedName,'conds',S.classSelector,'nVox_thresh',S.class_args.nVox, 'statmap_funct', S.statmap_funct, 'statmap_arg',statmap_arg);
                classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
            end
            
            classifier_pattern = S.preprocPatCondensedName;
            
            %%
%            scrap_save_and_writeout_subjstructs
            
            %penalty parameter loop
            for pnl = 1:length(S.penaltyParams)
                
                S.class_args.penalty = S.penaltyParams(pnl);
                
                if S.extractMeanSignal
                    
                    % instead of doing classification, just extract
                    % the mean signal in a given set of regions and
                    % compare them
                    [subj results] =  AG_extractMeanSignal(subj,S);
                else
                    %run the classification.
                    [subj results] = cross_validation(subj,classifier_pattern,'conds', S.classSelector, classifier_mask,S.class_args, 'perfmet_functs', S.perfmet_functs);               
                end
                
                %set up importance maps structures.
                if ((S.generate_importance_maps == 1) || (S.generateBetaMaps == 1))
                    for rif = 1:length(results.iterations);
                        results_IW{rif}.iterations(1).scratchpad.net.IW{1} = results.iterations(rif).scratchpad.w(2:end,:)';
                    end
                end
                
                if S.generateBetaMaps ==1
                    PM_generateBetaMaps(subj, results_IW, S )
                end
                
                if ~(exist(S.group_mvpa_dir, 'dir'))
                    mkdir(S.group_mvpa_dir);
                end
                
                %store results
                res.subj{b}.penalty(pnl).nVox(1).weights(TRW).iter{n} = results;
                res.subj{b}.penalty(pnl).nVox(1).weights(TRW).S = S;
                res.subjArray = subj_array;
                
                %save results
                save (fullfile(S.group_mvpa_dir, S.saveName), 'res');
                
                time2finish = toc/60;
                display(['Finished ' S.subj_id ' in ' num2str(time2finish) ' minutes']);
                
            end
            
            % generate importance maps.
            if S.generate_importance_maps
                PM_generate_importance_maps(subj, results, results_IW, S)
            end
        end
    end
end



















