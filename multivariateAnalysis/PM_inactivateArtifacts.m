function [ subj S ] = PM_inactivateArtifacts( subj, S )
% code to inactivate train trials with artifacts.

rawSigInts = S.artStruct.zscoreA_cell([1:length(S.TrainRuns) 1:length(S.TestRuns)]);
rawMotDerivs = S.artStruct.delta_cell([1:length(S.TrainRuns) 1:length(S.TestRuns)]);

catSigInts = vertcat(rawSigInts{:});
catMotDerivs = vertcat(rawMotDerivs{:});

Art = ((abs(catSigInts) > par.art.sigThresh) + (abs(catMotDerivs) > par.art.motThresh)) >0;

ArtWithNeighbors = unique([find(Art); find(Art)-1; find(Art)+1]);

regOnsets = find(all_regs);

artRegs = ismember(regOnsets, ArtWithNeighbors);

actives = actives_h;
actives((artRegs==1)' & (TrainTestOneIter==1))=0; %only remove train trials with artifacts.

S.actives = actives;

end

