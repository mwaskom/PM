function [ subj] = PM_run_pca( subj, S )
%Run PCA on the training data, project the testing data onto these
%components.

thisSel = get_mat(subj, 'selector', S.thisSelector);
thisPat = get_mat(subj, 'pattern', S.preprocPatCondensedName);

train_idx = find(thisSel==1);
test_idx = find(thisSel==2);



trainPats = thisPat(:,train_idx);
testPats = thisPat(:,test_idx);

[pca_coeffs pca_scores pca_latent] = lprincomp(trainPats');

pca_pctVarExp = cumsum(pca_latent)/sum(pca_latent);
idxSupThresh = find(pca_pctVarExp>S.pca_proportion_var);

nComponents = idxSupThresh(1);

pc_signal_train = pca_scores(:,1:nComponents)';
pc_signal_test = pca_coeffs(: , 1:nComponents)' * testPats;

subj = duplicate_object(subj,'pattern',S.preprocPatCondensedName,'principleComponents');
subj = set_mat(subj,'pattern','principleComponents',[pc_signal_train pc_signal_test],'ignore_diff_size',true);
subj.patterns{end}.derived_from = '';
subj.patterns{end}.masked_by = 'pca_mask';

subj = duplicate_object(subj,'mask',S.roi_name,'pca_mask');
subj.masks{end}.mat = 1:nComponents;
subj.masks{end}.group_name = 'pca_mask_group';
subj.masks{end}.derived_from = '';



end

