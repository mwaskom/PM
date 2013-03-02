function [M reh groupPsy idxG] = MnemonicBehAnalysisWholeshebang(SA)
%analyzes the behavioral data of a group of subjects who have performed the
%mnemonic decision making task

for i = 1:length(SA)
     par = PM_Params(SA(i), 'mnem', 0);
    
    if par.goodSub
    %analyze a given subject's localizer data
    [loc] = fMRIBehAnalysis_Loc(par);
    fnLoc = fieldnames(loc);
    
    for f = 1:length(fnLoc)
        M.loc.vals(i,f) = loc.(fnLoc{f});
        M.loc.names = fnLoc';
    end
    
        
    
    
    %analyze a given subject's retrieval data
    [ret, psyphys, idxG.sub{i}, reh(i)] = Mnemonic_fMRIBehAnalysis_Retrieval(par);
    
    idxG.sub{i}.subNo = par.subNo;
    
    fnRet = fieldnames(ret);
    
    for f = 1:length(fnRet)        
        M.ret.vals(i,f) = ret.(fnRet{f});
        M.ret.names = fnRet';
    end
    
    
    %group psychophysics data
    groupPsy(i).dat = psyphys;
    groupPsy(i).goodSub = par.goodSub;
    end
end

%group plot of psychophysics data
GroupPsychoPhysPlot(groupPsy)