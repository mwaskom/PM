function [res idx ] = fMRIBehAnalysis_Loc(par)

cd (par.behavdir);



loc.cor = {};
loc.cond = [];
loc.RT = [];
loc.coh = [];
loc.sess = [];

d = dir('*Loc*.mat');
d = d(find(cellfun(@(x) ~strcmp(x(1),'.'), {d(:).name})));
locFileIdx = 1:size(d);

numvols_loc = par.(par.task).numvols(par.scansSelect.(par.task).loc);
for r = 1:length(d)
    
    thisLoc = dir(['*Loc' '*' num2str(locFileIdx(r)) '*.mat']);
    thisLoc = thisLoc(find(cellfun(@(x) ~strcmp(x(1),'.'), {thisLoc(:).name})));
    RL(r).dat = load(thisLoc.name);
    
    loc.cor = horzcat(loc.cor, [RL(r).dat.trial_data.result]);
    loc.cond = horzcat(loc.cond, [RL(r).dat.trial_data.stim_group]);
    loc.RT = horzcat(loc.RT, [RL(r).dat.trial_data.rt]);
    loc.coh = horzcat(loc.coh, RL(r).dat.trial_data.stim_coh_seq);
    loc.sess = horzcat(loc.sess, r*ones(size([RL(r).dat.trial_data.result])));
    alltrials_h{r} = [RL(r).dat.trial_data.start_t] - [RL(r).dat.trial_data.time0] - 12 + 2*sum(numvols_loc(1:(locFileIdx(r)-1)));
end

idx.alltrials = round(horzcat(alltrials_h{:}))';

idx.face = (loc.cond==1);
idx.house = (loc.cond==2);


if strcmp(par.substr, 'pm_032811')
    idx.cor = ones(size(idx.alltrials))';
else
    idx.cor = (loc.RT < 2) .* (loc.RT>.2); %strcmp(loc.cor, 'CORRECT');
end

idx.inc = ~idx.cor;

idx.faceCor = idx.face .* idx.cor;
idx.houseCor = idx.house .* idx.cor;

idx.coh = loc.coh;

res.pctCorFace = sum(idx.face .* idx.cor)/sum(idx.face );
res.pctCorHouse = sum(idx.house .* idx.cor)/sum(idx.house );
res.pctCor = mean(idx.cor);

res.rtFaceCor = median(loc.RT(find(idx.cor .* idx.face)));
res.rtHouseCor = median(loc.RT(find(idx.cor .* idx.house)));

idx.sess = loc.sess;

end