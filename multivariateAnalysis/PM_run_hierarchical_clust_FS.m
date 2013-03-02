function subj = PM_run_hierarchical_clust_FS(subj, S)

pats = get_mat(subj, 'pattern', S.preprocPatCondensedName);
pats = get_masked_pattern(subj, S.preprocPatCondensedName, 'spiral_d_z_condensed_thresh1000_1');
regs = subj.regressors{1}.mat(1,:);
thisSel = get_mat(subj, 'selector', S.thisSelector);

idxTrain = (thisSel==1);
X = pats(:,idxTrain)';

T = clusterdata(X',1.1);


CL = 0;
clustPat = [];
for cl = 1:max(T)
    clustIdx = find(T==cl);
    if length(clustIdx)>5
        CL = CL+1;
        clustPat(CL,:) = mean(pats(clustIdx,:));
    end
end

subj = duplicate_object(subj,'pattern',S.preprocPatCondensedName,'hierClust');
subj = set_mat(subj,'pattern','hierClust',clustPat,'ignore_diff_size',true);
subj.patterns{end}.derived_from = '';
subj.patterns{end}.masked_by = 'hierClustMask';

subj = duplicate_object(subj,'mask',S.roi_name,'hierClustMask');
subj.masks{end}.mat = 1:CL;
subj.masks{end}.group_name = 'hierClust_mask_group';
subj.masks{end}.derived_from = '';

