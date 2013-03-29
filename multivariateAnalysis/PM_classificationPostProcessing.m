function [res, dat, groupPsy, datB] = PM_classificationPostProcessing(qqq, task)

groupPsy = [];
pnl = 1;
weights = 1;
useResiduals = false;
for s = 1:length(qqq.subj)
    
    subj_id = qqq.subjArray(s);
    
    % behavioral analysis
    if strcmp(task, 'mnem')
        par = PM_Params(subj_id, 'mnem', 0);
        [~, ~, idxB] = Mnemonic_fMRIBehAnalysis_Retrieval(par);
    elseif strcmp(task, 'perc')
        par = PM_Params(subj_id, 'perc', 0);
        [~, ~, idxB] = Perceptual_fMRIBehAnalysis(par);
    end
    
    % read in classification data
    resStruct = qqq.subj{s}.penalty(pnl).nVox.weights(weights).iter{1}.iterations;
    resS = qqq.subj{s}.penalty(pnl).nVox.weights(weights).S;
    
    % change path for files created on old filesystem
    if strcmp(resS.onsetsTestDir(1:27), '/Users/alangordon/mounts/w5')
        resS.onsetsTestDir = ['/biac4/wagner/biac3/wagner5' resS.onsetsTestDir(28:end)];
    end
    
    allOns = sort([resS.onsets_test_in_classifier{:}]);
      
    %global signal stuff
    globalSignalDat = load('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/mvpa_files/pre2013/MeanSigIntensityInOTCortex');
    globalSig = globalSignalDat.res.subj{s}.penalty.nVox.weights(weights).iter{1}.iterations.acts(1,:);
    
    %% extract raw data from the resStruct
    for i = 1:length(resStruct)
        
        cv.raw.correctsVec{i} = resStruct(i).perfmet.corrects; %#ok<*AGROW> %correct
        cv.raw.actsVec{i} = resStruct(i).acts(1,:) ; % evidence for class 1      
        cv.raw.actsVec2{i} = resStruct(i).acts(2,:) ; % evidence for class2

        cv.raw.desiredsVec{i} = resStruct(i).perfmet.desireds; % correct class
        
        cv.raw.actsOfTrueClass{i} = cv.raw.actsVec{i}; % evidence for true class
        cv.raw.actsOfTrueClass{i}(cv.raw.desiredsVec{i}~=1) = cv.raw.actsVec2{i}(cv.raw.desiredsVec{i}~=1);
        
        cv.raw.actsOfFalseClass{i} = cv.raw.actsVec2{i}; % evidence for false class
        cv.raw.actsOfFalseClass{i}(cv.raw.desiredsVec{i}~=1) = cv.raw.actsVec{i}(cv.raw.desiredsVec{i}~=1);
        
        cv.raw.actsDifference{i} = cv.raw.actsVec{i} - cv.raw.actsVec2{i}; % difference in evidence for true vs. false class
        cv.raw.actsDifference{i}(cv.raw.desiredsVec{i}~=1) = -1*cv.raw.actsDifference{i}(cv.raw.desiredsVec{i}~=1);
        
        cv.raw.guessesVec{i} = resStruct(i).perfmet.guesses; % binary classifier guesses
        cv.raw.signedActs{i} = abs(cv.raw.actsVec{i} - .5) .* (2*cv.raw.correctsVec{i} - 1); % acts signed in correct direction.
        
        if isfield(resStruct(i), 'test_idx')
            testIdx{i} = resStruct(i).test_idx;
        else
            testIdx{i} = 1:length(cv.raw.guessesVec{i});
        end
        
        cv.raw.signedEv{i} = logit(cv.raw.actsVec{i}); %take the logit of the probabilistic estimate for class 1
        cv.raw.unsignedEv{i} = logit(cv.raw.actsOfTrueClass{i});%take the logit of the probabilistic estimate for the correct class
        
        nTrainingTrials(i) = length(resStruct(i).train_idx);
    end

    testIdxCat = [testIdx{:}];
    [~, testIdxInv] = sort(testIdxCat);
        
    fn = fieldnames(cv.raw);
 
    % consoildate vectors across classifier iterations
    % reorder them to match chronological presentation time.
    for f = 1:length(fn)
        cv.consolidated.(fn{f}) = [cv.raw.(fn{f}){:}];
        cv.reordered.(fn{f}) = cv.consolidated.(fn{f})(testIdxInv);
    end     
    
    %compare onsets from qqq to onsets from idxB
    onsInClassifier = ismember(idxB.alltrials, allOns);
    
    %% get behavioral data
    if strcmp(task, 'mnem')
        idx.cor = idxB.cor(onsInClassifier);
    elseif strcmp(task, 'perc')
        idx.cor = idxB.cor(onsInClassifier);
    end
    
    idx.inc = idxB.inc(onsInClassifier);
    idx.high = idxB.high(onsInClassifier);
    idx.low = idxB.low(onsInClassifier);
    idx.face = idxB.face(onsInClassifier);
    idx.house = idxB.house(onsInClassifier);
    idx.rt = idxB.rt(onsInClassifier);
    idx.unsignedConf = idxB.unsignedConf(onsInClassifier);
    
    dat{s}.idx = idx;
    dat{s}.cv = cv;
    datB{s} = idxB;
    
    %% results summary
    res.cor.class1(s) = nanmean(cv.reordered.correctsVec(cv.reordered.desiredsVec==1));
    res.cor.class2(s) = nanmean(cv.reordered.correctsVec(cv.reordered.desiredsVec==2));
    res.cor.class1Cor(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==1) & (idx.cor==1) ));
    res.cor.class2Cor(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==2) & (idx.cor==1) ));
    res.cor.class1Inc(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==1) & (idx.inc==1) ));
    res.cor.class2Inc(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==2) & (idx.inc==1) ));
    res.cor.class1CorHighConf(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==1) & (idx.cor==1) & (idx.high==1)));
    res.cor.class2CorHighConf(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==2) & (idx.cor==1) & (idx.high==1)));
    res.cor.class1CorLowConf(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==1) & (idx.cor==1) & (idx.high==0)));
    res.cor.class2CorLowConf(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==2) & (idx.cor==1) & (idx.high==0)));
    res.cor.class1IncHighConf(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==1) & (idx.inc==1) & (idx.high==1)));
    res.cor.class2IncHighConf(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==2) & (idx.inc==1) & (idx.high==1)));
    res.cor.class1IncLowConf(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==1) & (idx.inc==1) & (idx.high==0)));
    res.cor.class2IncLowConf(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==2) & (idx.inc==1) & (idx.high==0)));
    
    res.MeanNTrainingTrials(s) = nanmean(nTrainingTrials);
    res.sub(s).S = resS;
   
    
    %% send all variables of interest into a structure that can be ported to R.
    fn = fieldnames(idx);
    for f = 1:length(fn)
       toR.(fn{f}){s} = asColumn(double(idx.(fn{f})));        
    end
    
    fn = fieldnames(cv.reordered);
    for f = 1:length(fn)
       toR.(fn{f}){s} = asColumn(double(cv.reordered.(fn{f})));        
    end
    
    toR.subs{s} = s*ones(size(toR.(fn{f}){s}));

end


%% write out data across subjects
fn = fieldnames(toR);

for f = 1:length(fn)
    out.(fn{f}) = vertcat(toR.(fn{f}){:});
    toPrint(:,f) = [fn{f}; num2cell(out.(fn{f}))];
end
    

cell2csv('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/ROutput/MnemonicData.csv', toPrint, ',', 2000);


end




