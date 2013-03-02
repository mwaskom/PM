function  PM_generateBetaMaps(subj, results_IW, S )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if S.generateBetaMaps
    betaDir = [S.importance_maps_dir '/' 'betas'];
    
    if ~exist(betaDir)
       mkdir(betaDir); 
    end
    
    for i = 1:length(results_IW)
        betas_h{i} = 10000*results_IW{i}.iterations(1).scratchpad.net.IW{1}(1,:);
    end
    
    betas = mean(vertcat(betas_h{:}),1);
    
    vol_info = S.vol_info;
    voxel_inds = find(subj.masks{end}.mat);
    vol_info.dir = betaDir;
    vol_info.fname = [ vol_info.dir '/' S.subj_id '_beta.img'];
    
    betaMap = zeros(vol_info.dim);
    betaMap(voxel_inds) = betas;
    spm_write_vol(vol_info,betaMap);
end

end

