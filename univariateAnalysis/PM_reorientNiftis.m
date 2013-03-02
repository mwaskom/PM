function [] = PM_reorientNiftis(par)

st.B = par.reorientmat;
P = par.rawNiis;

% par.RO.n = [1 1];
% 
% st.vols{1}.fname = P(1,:);
% st.vols{1}.n = par.RO.n;

mat = spm_matrix(st.B);
if det(mat)<=0
    spm('alert!','This will flip the images',mfilename,0,1);
end;

Mats = zeros(4,4,size(P,1));
spm_progress_bar('Init',size(P,1),'Reading current orientations',...
    'Images Complete');
for i=1:size(P,1),
    Mats(:,:,i) = spm_get_space(P(i,:));
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Init',size(P,1),'Reorienting images',...
    'Images Complete');
for i=1:size(P,1),
    spm_get_space(P(i,:),mat*Mats(:,:,i));
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
% tmp = spm_get_space([st.vols{1}.fname ',' num2str(st.vols{1}.n)]);
% if sum((tmp(:)-st.vols{1}.mat(:)).^2) > 1e-8,
%     spm_image('init',st.vols{1}.fname);
% end;



% Begin the extremely hacky part that nonetheless, seems to ensure that 
% all scans within a .nii file have the same orientation.  If anyone can figure out how to do
% this more elegantly, please let me know!!
for i =1:size(par.rawNiis,1)
    PMatName = strrep(par.rawNiis(i,:), '.nii', '.mat');
    if exist(PMatName, 'file')
        PMat = load(PMatName);
        for i=1:size(PMat.mat,3)
            PMat.mat(:,:,i) =  PMat.mat(:,:,1);
        end
        
        mat = PMat.mat;
        
        save(PMatName,'mat');
    end
end
