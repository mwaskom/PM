S.roi_file = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/group_analyses/analysis_perc_parmodByConfAndRT_16subs/corRT/mask.img';
S.roi_name = 'mask.img';
%S.img_files = 

% initialize subj structure
subj = init_subj('pairedT','pairedT15Subs');

% load mask file
subj = load_spm_mask(subj,S.roi_name,S.roi_file);

% load functional data
subj = load_analyze_pattern(subj,'perc',S.roi_name, S.img_files,'single',true);

subj = duplicate_object(subj,'pattern','ret','ret_zscore');

for i=1:17
   thisPat = subj.patterns{2}.mat(:,i);
   thisPatNonZero = thisPat(thisPat~=0);
   zscoredPat = zscore(thisPatNonZero);
   subj.patterns{4}.mat(thisPat~=0,i) = zscoredPat;    
end

[h p ci stats] = ttest(subj.patterns{3}.mat',subj.patterns{4}.mat');

stats.tstat


vol_info = S.vol_info;

thisMap = zeros(vol_info.dim);
voxel_inds = find(subj.masks{end}.mat);
    
thisMap(voxel_inds) = stats.tstat;
    
vol_info.dir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/mvpa_results/ImpMaps_PercVsRet';
vol_info.fname = [ vol_info.dir '/' 'ttestPercVsRet.img'];
vol_info.dt = [16 0];
if isempty(dir([vol_info.dir]))
    mkdir(vol_info.dir);
end

spm_write_vol(vol_info,thisMap);
