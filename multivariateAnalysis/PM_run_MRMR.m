function subj = PM_run_MRMR(subj, S)
%run mrmr feature selection.

idxTrain = find(subj.selectors{end}.mat==1);

D = subj.patterns{end}.mat(:,idxTrain)';
F = subj.regressors{1}.mat(1,idxTrain);
K = S.nMRMR;


[fea] = mrmr_mid_d(D, F, K);

subj = duplicate_object(subj,'mask','rOT_wfu.m.img','FS_MRMR_mask');
msk_idx = find(subj.masks{1}.mat);
msk_idx_FS = msk_idx(fea);
subj.masks{2}.mat(subj.masks{2}.mat==1) = 0;
subj.masks{2}.mat(msk_idx_FS) = 1;

end

