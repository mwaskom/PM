function [res] = PM_classificationPostProcessingLoc(qqq)

for s = 1:length(qqq.subj)

% subj_id = qqq.subj{s}.penalty.nVox.weights.iter{1}.header.subj_id;
% par = PM_mvpa_params(subj_id);
% %Design a function that ports the data to R.
% 
% mvpa_ons = load(fullfile(par.onsetsTestDir, 'mvpa_ons'));
% [~, testRegsIdx] = ismember( par.condsTest, mvpa_ons.names);
% 
% [allOns_h aOI] = sort(vertcat(mvpa_ons.onsets{testRegsIdx}));
% 
% idx.enoughTRs_h = vertcat(qqq.subj{s}.S.enoughTRs.train{:});
% idx.enoughTRs = idx.enoughTRs_h(aOI);
% 
% allOns = allOns_h(find(idx.enoughTRs));

for i = 1:length(qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations)
    % for each classification iteration, the vector of whether the
    % classification was correct or not
    correctsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).perfmet.corrects; 
    actsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).acts(1,:);
    desiredsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).perfmet.desireds; 
end

correctsVecCat = [correctsVec{:}];
actsVecCat = [actsVec{:}];
desiredsVecCat = [desiredsVec{:}];


res.cor.class1(s) = mean(correctsVecCat(desiredsVecCat==1));
res.cor.class2(s) = mean(correctsVecCat(desiredsVecCat==2));
end


sprintf('end');







