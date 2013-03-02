function DK1_meanFuncs(subpar)

if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = par_params(subpar);
end

if ~exist(par.meanfuncdir)
    mkdir(par.meanfuncdir);
end

cd (par.meanfuncdir);

spm_imcalc_ui(par.rascanfiles, fullfile(par.meanfuncdir, 'mean_rascan.img'),'sum(X)', {'dmtx'});