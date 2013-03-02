function [subj S] = PM_linregsetup(subj S)
if strcmp(S.thisSelector, 'TrainTestOneIterGroup')
    S.actives = ones(size(TrainTestOneIter));
    %actives(TrainTestOneIter==1) = idxB.corWithZeros';
    %idxTrainInClassifier = ismember(idxTr.alltrials, [S.onsets_test_in_classifier{:}]);
    idxTestInClassifier = ismember(idxTe.alltrials, [S.onsets_test_in_classifier{:}]);
    
    condensed_regs_of_interest = zeros(size(TrainTestOneIter));
    condensed_regs_of_interest(TrainTestOneIter==1) = idxTr.coh_signed(S.idxOnsets_train_in_classifier)';
    condensed_regs_of_interest(TrainTestOneIter==2) = idxTe.face(S.idxOnsets_test_in_classifier) - idxTe.house(S.idxOnsets_test_in_classifier);
else
    S.actives = ones(size(idxB.corWithZeros'));
    condensed_regs_of_interest = idxB.coh_signed(idxTr.corWithZeros==1)';
end

end

