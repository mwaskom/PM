function par_coreg(refimg, sourceimg)
% coregisters...  separated from parbatch for modularity...


crX = spm_coreg(refimg, sourceimg);
M  = inv(spm_matrix(crX));
PO = sourceimg;
MM = zeros(4,4,size(PO,1));
for j=1:size(PO,1),
	MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
end;
for j=1:size(PO,1),
	spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
end;
