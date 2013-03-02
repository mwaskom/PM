function [S] = PM_mvpa_onsets_and_images(S)


if S.xval
    [onsets_train, S.idxOnsets_train_in_classifier] = MakeOnsets(S, 'train');
    [onsets_test, S.idxOnsets_test_in_classifier] = MakeOnsets(S, 'test');
    
    for k=1:length(onsets_train)
        S.onsets{k} = union(onsets_train{k}, onsets_test{k});
    end
else
    [onsets_train, S.idxOnsets_train_in_classifier] = MakeOnsets(S, 'train');
    [onsets_test, S.idxOnsets_test_in_classifier] = MakeOnsets(S, 'test');
    
    for k=1:length(onsets_train)
        S.onsets{k} = union(onsets_train{k}, S.durTrain + onsets_test{k});
    end
end

S.onsets_train_in_classifier = onsets_train;
S.onsets_test_in_classifier = onsets_test;

end

function [onsets idxOnsets_in_classifier] = MakeOnsets(S, part)

if strcmp(part, 'train')
    thisOns = load(fullfile(S.onsetsTrainDir, 'mvpa_ons'));
    thisConds = S.condsTrain;
    theseFiles = S.filenames_train;
    theseTRWeights = S.TR_weights{1};
    thisIdx = S.idxTr;
elseif strcmp(part, 'test')
    thisOns = load(fullfile(S.onsetsTestDir, 'mvpa_ons'));
    thisConds = S.condsTest;
    theseFiles = S.filenames_test;
    theseTRWeights = S.TR_weights{2};
    thisIdx = S.idxTe;
end

onsets = cell(size(thisConds));
for i = 1:length(thisConds)
    for j = 1:length(thisConds{i})
        idxThisCond = find(strcmp (thisOns.names, thisConds{i}{j}));
        if ~isempty(idxThisCond)
            time_idx = floor(thisOns.onsets{idxThisCond}/S.TR) + 1;
            enoughTRs_h = (time_idx + length(theseTRWeights)) <= length(theseFiles);
            %enoughTRs{i} = logical(vertcat(enoughTRs{i}, enoughTRs_h));
            theseOnsets = asRow(thisOns.onsets{idxThisCond});
            if strcmp(S.patternType, 'betas')
                onsets{i} = sort(horzcat(onsets{i}, theseOnsets));
            else
                onsets{i} = sort(horzcat(onsets{i}, theseOnsets(enoughTRs_h)));
            end
        end
    end
end

idxOnsets_in_classifier = asRow(ismember(thisIdx.alltrials, [onsets{:}]));

end
