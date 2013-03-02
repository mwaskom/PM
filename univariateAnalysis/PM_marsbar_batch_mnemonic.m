function subjects_tc = PM_marsbar_batch_mnemonic()

conds = { 'face_cor_high' 'face_cor_low' 'house_cor_low' 'house_cor_high'};
%condType = {'face' 'house'};
expt_dir = '/Users/alangordon/mounts/wagner5/alan/perceptMnemonic/fmri_data';
%roi_names = {'faceVsHouse_corOnly_28_-25_-34_roi.mat' 'House_vs_Face_corOnly_25_-18_-26_roi.mat' 'faceVsHouse_corOnly_-41_-25_-30_roi.mat'	'House_vs_Face_corOnly_-31_-15_-26_roi.mat'};
%roi_names = {'left_p001_k5_faceVsHouse_-38_-18_-42_roi.mat' 'left_p001_k5_houseVsFace_-28_-28_-30_roi.mat' 'right_p001_k5_faceVsHouse_31_-21_-46_roi.mat' 'right_p001_k5_houseVsFace_25_-25_-30_roi.mat'};

subjects = { 'pm_120810'    'pm_011611'    'pm_012211'    'pm_012211b'   'pm_012611'};
 
TRs = 3:4;


for curr_subj=1:length(subjects), % go through the list of subjects
    
    roi_dir = fullfile(expt_dir, subjects{curr_subj}, 'analysis_loc');
    
    roi_names_h = dir(fullfile(roi_dir, '*roi.mat'));
    roi_names = {roi_names_h.name}';
    
    fprintf(1,'extracting TC from subject %s\n',subjects{curr_subj});
    spm_name = ['/Users/alangordon/mounts/wagner5/alan/perceptMnemonic/fmri_data/' subjects{curr_subj} '/analysis_ret_2Class_2Perf_2Conf/SPM.mat'];
    D = mardo(spm_name);
    D = autocorr(D,'fmristat',2);
    for curr_clust=1:length(roi_names) % go through the list of clusters
        
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
                
                fir_tc(:, c) = event_fitted_fir(E, e_specs(:,e_s), bin_size, bin_no, opts);
                
                [maxC, maxI] = max(fir_tc(2:((fir_length/bin_size)-1),c));
                
                subjects_tc(curr_subj).rois(curr_clust).conds = conds;
                
                subjects_tc(curr_subj).rois(curr_clust).roiName = roi_names{curr_clust};
                %subjects_tc(curr_subj).rois(curr_clust).maxval(c) = maxC;
                
                subjects_tc(curr_subj).rois(curr_clust).maxval(c) = maxC;
                neighbors = intersect(1:length(fir_tc), [maxI-1: maxI+1]);
                
                subjects_tc(curr_subj).rois(curr_clust).meanTRs(c) = mean(fir_tc(TRs,c));
                
                subjects_tc(curr_subj).rois(curr_clust).raw = fir_tc;
            end
        %end
        
        
    end % end clusters loop
    
    subjects_tc(curr_subj).conds = conds;
    subjects_tc(curr_subj).e_specs = e_specs;
    subjects_tc(curr_subj).e_names = e_names;
    subjects_tc(curr_subj).n_events = n_events;
    
    
end % end subjects loop


%%
for i = 1:length(subjects_tc)
    for j = 1:length(subjects_tc(i).rois)
        for k = 1:length(subjects_tc(i).rois(j).meanTRs)
            M(i,j,k) = subjects_tc(i).rois(j).meanTRs(k);
        end
    end
end

faceTCs = squeeze(mean(M(:,[1 3],:),2));
houseTCs = squeeze(mean(M(:,[2 4],:),2));
            
                                    
                                   
