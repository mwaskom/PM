function X = PM_spm_orthWithFirstColumn(X,OPT)
% recursive Gram-Schmidt orthogonalisation of each subsequent column with
% only the first column




x = X(:,1);
j=1;
for i=2:size(X,2)
   D = spm_orth ([X(:,1) X(:,i)]);
   x = [x D(:,2)];
end

X = x;
