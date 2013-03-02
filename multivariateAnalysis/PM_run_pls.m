function [subj S] = PM_run_pca(subj, S)

pats = get_mat(subj, 'pattern', S.preprocPatCondensedName);
regs = subj.regressors{1}.mat(1,:);
thisSel = get_mat(subj, 'selector', S.thisSelector);

idxTrain = (thisSel==1);
idxTest = (thisSel==2);


X = pats(:,idxTrain)';
Y = regs(1,idxTrain)';
[XL, YL, XS, YS, BETA, PCTVAR] = plsregress(X,Y,S.nPlsComps);

pls_pats = (XL' * pats);


subj = duplicate_object(subj,'pattern',S.preprocPatCondensedName,'pls');
subj = set_mat(subj,'pattern','pls',pls_pats,'ignore_diff_size',true);
subj.patterns{end}.derived_from = '';
subj.patterns{end}.masked_by = 'pls_mask';

subj = duplicate_object(subj,'mask',S.roi_name,'pls_mask');
subj.masks{end}.mat = 1:S.nPlsComps;
subj.masks{end}.group_name = 'pls_mask_group';
subj.masks{end}.derived_from = '';

preds = BETA(1) + BETA(2:end)' * pats(:,idxTest);

pred_bin = (preds>.5);

regsTest = regs(idxTest);
accuracy = (sum((regsTest==1).*(pred_bin==1)) + sum((regsTest==0).*(pred_bin==0)))/length(regsTest);

S.pLSDA_Acc = accuracy;
end