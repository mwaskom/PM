file{1} = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/group_analyses/patLinRegAnalysis_Perceptual_standardized_FaceHouseSeparate_cohRidgeRegressionFaceEv_14Subs/percFaces_FaceEv/spmT_0001.img';
file{2} = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/group_analyses/patLinRegAnalysis_Perceptual_standardized_FaceHouseSeparate_cohRidgeRegressionHouseEv_14Subs/percHouses_HouseEv/spmT_0001.img';

outputDir = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/group_analyses/conjAnalysis';

outputName = 'conj_faceByFaceLinRegEv_and_houseByHouseLinRegEv_conjoint_PM_imcalc.img';

vol{1} = spm_vol(file{1});
vol{2} = spm_vol(file{2});

v{1} = spm_read_vols(vol{1});
v{2} = spm_read_vols(vol{2});

vol_info = spm_vol(file{1});
    
outMap=  zeros(vol_info.dim);

voxel_inds = find(v{1}>0);

for i = 1:length(voxel_inds)
    i2 = voxel_inds(i);
    outMap(i2) = (v{1}(i2)>2.031) .* (v{2}(i2)>2.031);    
end


vol_info.dir = outputDir;
vol_info.fname = [ vol_info.dir '/' outputName];

if isempty(dir([vol_info.dir]))
    mkdir(vol_info.dir);
end

spm_write_vol(vol_info,outMap);

fprintf('\n wrote out %s \n', vol_info.fname);

