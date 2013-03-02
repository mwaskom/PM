function [subj] = PM_mvpa_load_and_preprocess_raw_data_tmp(S)

    % initialize subj structure
    subj = init_subj(S.exp_name,S.subj_id);

    % load mask file
    subj = load_spm_mask(subj,S.roi_name,S.roi_file);
    
    % load functional data
    subj = load_analyze_pattern(subj,'spiral',S.roi_name, S.img_files,'single',true);

    % move pattern to hard drive
    % subj = move_pattern_to_hd(subj, 'spiral');

    % make runs vector
    subj = init_object(subj,'selector','runs');

%     trial_idx = 1;
%     for r = 1:num_runs
%         runs(trial_idx:num_TP_per_run*r) = r;
%         trial_idx = trial_idx + num_TP_per_run;
%     end 

    trial_idx = 1;
    for r = 1:length(S.runs_vector)
        runs(trial_idx:trial_idx+S.runs_vector(r)-1)=r;
        trial_idx = trial_idx+S.runs_vector(r);
    end
    
    
    
    subj = set_mat(subj,'selector','runs',runs);
    
    
 
    
    % detrend the timeseries data
    %subj = detrend_runs(subj,'spiral','runs');  % not in mvpa tutorial, but seems important to do
    
    %subj = remove_mat(subj,'pattern','spiral');  
    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d');
    
    % clean up workspace
    %subj = remove_mat(subj,'pattern','spiral');
    
    % high-pass filter the timeseries data
    subj = hpfilter_runs(subj,'spiral','runs',100,2); % remove frequencies below .01 Hz %change back to 100.  AG 06/22/10.
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral');
    
    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d_hp');
    
    % zscore the data from each run
    subj = zscore_runs(subj,'spiral_hp','runs'); % Z-score the data
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_hp');  
    
    if S.denoise
        %NOTE: use white matter voxels here, too
        
        data_h = single(subj.patterns{end}.mat);
        data_h = reshape(data_h, size(data_h,1), 1, 1, size(data_h,2));
        for i=unique(runs)
            data{i} = data_h(:,:,:,runs==i);
            onsetsTheseRuns = find(runs==i);
            
            thisDesign = zeros(length(onsetsTheseRuns), S.num_conds);
            for j=1:length(S.onsets_train_in_classifier)
                
                onsetsThisCond = sort(round(S.onsets{j}/S.TR));
                theseOnsets = intersect(onsetsThisCond,onsetsTheseRuns);
                theseOnsets_fromBeginningOfRun = theseOnsets - onsetsTheseRuns(1) + 1;
                thisDesign(theseOnsets_fromBeginningOfRun,j) = 1;
                design{i} = thisDesign;
            end
            
        end
        [results,denoiseddata] = GLMdenoisedata(design,data,.5,S.TR,[],[],[],[]);
        
        cat_DND = squeeze(cat(4,denoiseddata{:}));
        subj.patterns{end}.name = 'spiral_dn';
        subj.patterns{end}.mat = cat_DND; %Change this, this is a major hack!
        %subj = duplicate_object(subj,'pattern','spiral_hp_z','spiral_hp_z_dn');
        %subj = set_mat(subj,'pattern','spiral_hp_z_dn',cat_DND,'ignore_diff_size',true);
    end
    
    if ~exist(S.workspace_dir)
        mkdir(S.workspace_dir);
    end
    
    cd (S.workspace_dir);
    
    save (S.workspace);
    %save (fullfile(S.mvpa_dir, [S.subj_id '_' S.roi_name '.mat']));
    
    %save_cmd = ['save ' S.mvpa_dir '/' S.subj_id '_' S.roi_name '.mat'];
    %eval(save_cmd);
    
    % save final pattern in single precision form (8 sig figs) to save RAM and HD space    
    % subj.patterns{end}.mat = single(subj.patterns{end}.mat);