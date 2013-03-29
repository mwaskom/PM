function [M groupPsy gRes idxG] = PerceptualBehAnalysisWholeshebang(SA)
%analyzes the behavioral data of a group of subjects who have performed the
%mnemonic decision making task
i2 = 0;
for i = 1:length(SA)
    par = PM_Params(SA(i), 'perc', 0);
    
    if ~isempty(par.substr)
    i2 = i2+1;    
    %analyze a given subject's localizer data
     [loc] = fMRIBehAnalysis_Loc(par);
     fnLoc = fieldnames(loc);
     
     for f = 1:length(fnLoc)
         M.loc.vals(i,f) = loc.(fnLoc{f});
         M.loc.names = fnLoc';
     end
    
        
    gRes = [];
    %analyze a given subject's perceptual data
    [res, psyphys, idxG.sub{i}] = Perceptual_fMRIBehAnalysis(par);
    %[res, psyphys] = Perceptual_fMRIBehAnalysis(par);
    idxG.sub{i}.subNo = par.subNo;
     
    fnRet = fieldnames(res);
    
    for f = 1:length(fnRet)        
        M.res.vals(i,f) = res.(fnRet{f});
        M.res.names = fnRet';
    end
    
    
    %group psychophysics data
    groupPsy(i).dat = psyphys;
    groupPsy(i).goodSub = par.goodSub;
    end
end

%group plot of psychophysics data
 %Mnemonic_GroupPsychoPhysPlot(groupPsy)