function [ ] = PM_copyNiis(par)


%% move  .nii files from raw to functional dir.
if exist(par.rawdir)
    cd (par.rawdir);
    
    if ~exist(par.funcdir)
        mkdir(par.funcdir);
    end
    
    dr = dir('run*.nii');
    dm = dir('mean.nii');
    dv = dir('valid.nii');
    
    for i = 1:length(dr)
        unix (['mv ' dr(i).name ' ' fullfile(par.funcdir, dr(i).name)])
    end
    
    unix (['mv ' dm.name ' ' fullfile(par.funcdir, dm.name)])
    unix (['mv ' dv.name ' ' fullfile(par.funcdir, dv.name)])

else
    warning(['no valid raw dir for subject ' par.substr]);
end

% %% make a copy of the .nii files and prepend them with 'R'.
if exist(par.funcdir)
    cd (par.funcdir);
    
    dr = dir('run*.nii');
    dm = dir('mean.nii');
    dv = dir('valid.nii');
    
    for i = 1:length(dr)
        unix (['cp ' dr(i).name ' R' dr(i).name]);
    end
    
    unix (['cp ' 'valid.nii' ' Rvalid.nii']);
    unix (['cp ' 'mean.nii' ' Rmean.nii']);
    
else
    warning(['no valid func dir for subject ' par.substr]);
end

