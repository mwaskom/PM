function [res, groupPsy, dat, datB] = PM_classificationPostProcessingMnemonicImproved(qqq)

groupPsy = [];
pnl = 1;
weights = 1;
useResiduals = false;
for s = 1:length(qqq.subj)
    
    subj_id = qqq.subjArray(s);
    %[par idxB] = PM_mvpa_params(subj_id, 'mnem');
    
    par = PM_Params(subj_id, 'mnem', 0);
    
    [~, ~, idxB] = Mnemonic_fMRIBehAnalysis_Retrieval(par);
    
    resStruct = qqq.subj{s}.penalty(pnl).nVox.weights(weights).iter{1}.iterations;
    resS = qqq.subj{s}.penalty(pnl).nVox.weights(weights).S;
    
    if strcmp(resS.onsetsTestDir(1:27), '/Users/alangordon/mounts/w5')
        resS.onsetsTestDir = ['/biac4/wagner/biac3/wagner5' resS.onsetsTestDir(28:end)];
    end
    
    mvpa_ons = load(fullfile(resS.onsetsTestDir, 'mvpa_ons'));
    condsTest = [resS.condsTest{:}];
    [~, testRegsIdx] = ismember(condsTest, mvpa_ons.names);
    testRegsIdx(testRegsIdx==0) = [];
    
    %[allOns ~] = sort(vertcat(mvpa_ons.onsets{testRegsIdx}));
    allOns = sort([resS.onsets_test_in_classifier{:}]);
    
    
    globalSignalDat = load('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/mvpa_files/MeanSigIntensityInOTCortex');
    globalSig = globalSignalDat.res.subj{s}.penalty.nVox.weights(weights).iter{1}.iterations.acts(1,:);

    
    for i = 1:length(resStruct)

        correctsVec{i} = resStruct(i).perfmet.corrects; %#ok<*AGROW>
        actsVec{i} = resStruct(i).acts(1,:) ;
       
        if size(resStruct(i).acts,1)==2
            actsVec2{i} = resStruct(i).acts(2,:) ;
        else
            actsVec2{i} = -1*resStruct(i).acts(1,:) ;
        end
        
        desiredsVec{i} = resStruct(i).perfmet.desireds;
        
        actsOfTrueClass{i} = actsVec{i};
        actsOfTrueClass{i}(desiredsVec{i}~=1) = actsVec2{i}(desiredsVec{i}~=1);
        
        actsOfFalseClass{i} = actsVec2{i};
        actsOfFalseClass{i}(desiredsVec{i}~=1) = actsVec{i}(desiredsVec{i}~=1);
        
        actsDifference{i} = actsVec{i} - actsVec2{i};
        actsDifference{i}(desiredsVec{i}~=1) = -1*actsDifference{i}(desiredsVec{i}~=1);
        
        guessesVec{i} = resStruct(i).perfmet.guesses;
        signedActs{i} = abs(actsVec{i} - .5) .* (2*correctsVec{i} - 1);
        
        if isfield(resStruct(i), 'test_idx')
            testIdx{i} = resStruct(i).test_idx;
        else
            testIdx{i} = 1:length(guessesVec{i});
        end
        
        nTrainingTrials(i) = length(resStruct(i).train_idx);
    end
    
    testIdxCat = [testIdx{:}];
    [~, testIdxInv] = sort(testIdxCat);
    
    desiredsVecCat = [desiredsVec{:}];
    signedActsVecCat = [signedActs{:}];
    if length(resS.condsTrain)>2
        actsVecCat_h = [actsVec{:}];
        actsVecCat2_h = [actsVec2{:}];
        
        actsVecCat = actsVecCat_h ./ (actsVecCat_h + actsVecCat2_h );
        
        guessesVecCat = 2 - (actsVecCat>.5);
        correctsVecCat = guessesVecCat==desiredsVecCat;
        actsVecTrueClassCat = [actsOfTrueClass{:}];
        actsVecFalseClassCat = [actsOfFalseClass{:}];
        actsVecDiffCat = [actsDifference{:}];
    else
        actsVecTrueClassCat = [actsOfTrueClass{:}];
        actsVecFalseClassCat = [actsOfFalseClass{:}];
        actsVecDiffCat = [actsDifference{:}];
        actsVecCat = [actsVec{:}];
        correctsVecCat = [correctsVec{:}];
        guessesVecCat = [guessesVec{:}];
    end
    
    % re-order the classifier-derived values to match chronological order
    actsVecTrueClassCat = actsVecTrueClassCat(testIdxInv);
    actsVecFalseClassCat = actsVecFalseClassCat(testIdxInv);
    actsVecDiffCat = actsVecDiffCat(testIdxInv);
    actsVecCat = actsVecCat(testIdxInv);
%     correctsVecCat = correctsVecCat(testIdxInv);
    %guessesVecCat = guessesVecCat(testIdxInv);
    signedActsVecCat = signedActsVecCat(testIdxInv);
    desiredsVecCat = desiredsVecCat(testIdxInv);
    correctsVecCat = correctsVecCat(testIdxInv);
    
    absActsVecCat = abs(actsVecCat - .5);
    actsVecCatLogit = logit(actsVecCat);
    
    actsVecLogitUnsigned = actsVecCatLogit;
    actsVecLogitUnsigned(desiredsVecCat==2) = -1 * actsVecLogitUnsigned(desiredsVecCat==2);
    
    
    %% signal intensity stuff
    if useResiduals        
        rstats = regstats(actsVecTrueClassCat,globalSig);
        actsVecTrueClassCat = rstats.r';        
    end
    
    %compare onsets from qqq to onsets from idxB
    idx.onsInClassifier = ismember(idxB.alltrials, allOns);
    
    classPat = {'face' 'house'};
    class = {'face' 'house'};
    confLevel = {'high' 'low'};
    conf = {'1' '2' '3' '4'};
    
    idx.cor = idxB.cor(idx.onsInClassifier);
    idx.inc = idxB.inc(idx.onsInClassifier);
    idx.high = idxB.high(idx.onsInClassifier);
    idx.low = idxB.low(idx.onsInClassifier);
    idx.face = idxB.face(idx.onsInClassifier);
    idx.house = idxB.house(idx.onsInClassifier);
    idx.rt = idxB.rt(idx.onsInClassifier);
    idx.cresp = idxB.cresp(idx.onsInClassifier);
    idx.unsignedActs = actsVecLogitUnsigned;
    
    dat{s} = idx;
    datB{s} = idxB;
    
    % for organizing by signed binary confidence
    j = 0;
    for i = 1:length(class)
        for c = 1:length(confLevel)
            j = j+1;
            thisClass = class{i};
            thisConf = confLevel{c};
            res.names{j} = [class{i}  '_' confLevel{c} 'conf' ] ;
            if sum(idx.(thisClass) .* idx.(thisConf) .* idx.cor) > 0;
                idx.thisIt{j} = idx.(thisClass) .* idx.(thisConf) .* idx.cor;
                medianAct(j) = median(actsVecCatLogit(idx.thisIt{j}==1));
                classPerf(j) = mean(actsVecLogitUnsigned(idx.thisIt{j}==1)>0);
            else
                medianAct(j) = NaN;
                classPerf(j) = NaN;
            end
        end
    end
    
 
    j=0;
    for i = 1:length(class)
        for c = 1:length(conf)
            j = j+1;
            thisClass = class{i};
            thisConf = conf{c};
            res.names{j} = [class{i}  '_' thisConf 'conf' ] ;
            if sum(idx.(thisClass) .* (idx.cresp==c) .* idx.cor) > 0;
                idx.thisIt{j} = idx.(thisClass) .* (idx.cresp==c) .* idx.cor;
                medianActsFull(j) = median(actsVecCatLogit(idx.thisIt{j}==1));
            else
                medianActsFull(j) = NaN;
            end
        end
    end
    
    
    res.cor.class1(s) = nanmean(correctsVecCat(desiredsVecCat==1));
    res.cor.class2(s) = nanmean(correctsVecCat(desiredsVecCat==2));
    res.cor.class1Cor(s) = nanmean(correctsVecCat((desiredsVecCat==1) & (idx.cor==1) ));
    res.cor.class2Cor(s) = nanmean(correctsVecCat((desiredsVecCat==2) & (idx.cor==1) ));
    res.cor.class1Inc(s) = nanmean(correctsVecCat((desiredsVecCat==1) & (idx.inc==1) ));
    res.cor.class2Inc(s) = nanmean(correctsVecCat((desiredsVecCat==2) & (idx.inc==1) ));
    res.cor.class1CorHighConf(s) = nanmean(correctsVecCat((desiredsVecCat==1) & (idx.cor==1) & (idx.high==1)));
    res.cor.class2CorHighConf(s) = nanmean(correctsVecCat((desiredsVecCat==2) & (idx.cor==1) & (idx.high==1)));
    res.cor.class1CorLowConf(s) = nanmean(correctsVecCat((desiredsVecCat==1) & (idx.cor==1) & (idx.high==0)));
    res.cor.class2CorLowConf(s) = nanmean(correctsVecCat((desiredsVecCat==2) & (idx.cor==1) & (idx.high==0)));
    res.cor.class1IncHighConf(s) = nanmean(correctsVecCat((desiredsVecCat==1) & (idx.inc==1) & (idx.high==1)));
    res.cor.class2IncHighConf(s) = nanmean(correctsVecCat((desiredsVecCat==2) & (idx.inc==1) & (idx.high==1)));
    res.cor.class1IncLowConf(s) = nanmean(correctsVecCat((desiredsVecCat==1) & (idx.inc==1) & (idx.high==0)));
    res.cor.class2IncLowConf(s) = nanmean(correctsVecCat((desiredsVecCat==2) & (idx.inc==1) & (idx.high==0)));
    res.MeanNTrainingTrials(s) = nanmean(nTrainingTrials);
    
    %     
     orderingVec = [3 4 2 1]; %reorder the res variables to go from high confidence house to high confidenc face
     orderingVec8 = [8 7 6 5 1 2 3 4];
     
     res.coh.medianActs(:,s) = medianAct(orderingVec);
     res.coh.medianActsSE(:,s) = .001*ones(size(medianAct(orderingVec)));
%    
     res.coh.classPerf(:,s) = classPerf(orderingVec);
     res.coh.classPerfSE(:,s) = .001*ones(size(classPerf(orderingVec)));
    
     res.coh.medianActsFull(:,s) = medianActsFull(orderingVec8);
     res.coh.medianActsFullSE(:,s) = .001*ones(size(medianActsFull(orderingVec8)));
    %     groupPsy(s).dat.pFaceGuessByConf.mean = res.coh.pGuessFace(orderingVec,s)';
    %     groupPsy(s).dat.pFaceGuessByConf.SE = res.coh.pGuessFaceSE(orderingVec,s);
    %     groupPsy(s).dat.pFaceGuessByConf.ticks.x = 0:5;
    %     groupPsy(s).dat.pFaceGuessByConf.ticks.y = [0 1];
    %     groupPsy(s).dat.pFaceGuessByConf.xlabel = 'Confidence (high house to high face)';
    %     groupPsy(s).dat.pFaceGuessByConf.ylabel = 'pFaceGuess' ;
    %     groupPsy(s).dat.pFaceGuessByConf.xVals = 1:4;
    %     groupPsy(s).dat.pFaceGuessByConf.marker = 'o-';
    
%     
 groupPsy(s).dat.medActByConf.mean = res.coh.medianActs(:,s)';
 groupPsy(s).dat.medActByConf.SE = res.coh.medianActsSE(:,s)';
 groupPsy(s).dat.medActByConf.ticks.x = 0:5;
 groupPsy(s).dat.medActByConf.ticks.y = [-2 2];
 groupPsy(s).dat.medActByConf.xlabel = 'Confidence (high house to high face)';
 groupPsy(s).dat.medActByConf.ylabel = 'Median Classifier Output' ;
 groupPsy(s).dat.medActByConf.xVals = 1:4;
 groupPsy(s).dat.medActByConf.marker = 'o-';
 
 groupPsy(s).dat.medActByConfFull.mean = res.coh.medianActsFull(:,s)';
 groupPsy(s).dat.medActByConfFull.SE = res.coh.medianActsFullSE(:,s)';
 groupPsy(s).dat.medActByConfFull.ticks.x = 0:9;
 groupPsy(s).dat.medActByConfFull.ticks.y = [-2 2];
 groupPsy(s).dat.medActByConfFull.xlabel = 'Confidence (high house to high face)';
 groupPsy(s).dat.medActByConfFull.ylabel = 'Median Classifier Output' ;
 groupPsy(s).dat.medActByConfFull.xVals = 1:8;
 groupPsy(s).dat.medActByConfFull.marker = 'o-';
% 

%groupPsy(s).dat.AccuracyByUnsignedConf.mean = res.coh.classPerf(:,s)';
%groupPsy(s).dat.AccuracyByUnsignedConf.SE = res.coh.classPerfSE(:,s)';
%groupPsy(s).dat.AccuracyByUnsignedConf.ticks.x = 0:5;
%groupPsy(s).dat.AccuracyByUnsignedConf.ticks.y = [0 1];
%groupPsy(s).dat.AccuracyByUnsignedConf.xlabel = 'Confidence (high house to high face)';
%groupPsy(s).dat.AccuracyByUnsignedConf.ylabel = 'pCorrect' ;
%groupPsy(s).dat.AccuracyByUnsignedConf.xVals = 1:4;
%groupPsy(s).dat.AccuracyByUnsignedConf.marker = 'o-';
    
    groupPsy(s).goodSub = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    toR.allBinaryConf_h{s} = double(idx.high');
    toR.allSubs_h{s} = s*ones(size(idx.high))';
    toR.allClass_h{s} = double(idx.face');
    toR.allUnsignedAct_h{s} = actsVecLogitUnsigned';
    toR.allSignedAct_h{s} = actsVecCatLogit';
    toR.allRT_h{s} = idx.rt';
    toR.allCresp{s} = idx.cresp';
    toR.allCorResp{s} = idx.cor';
    toR.allActsTrueClass{s} = actsVecTrueClassCat';
    toR.allActsFalseClass{s} = actsVecFalseClassCat';
    toR.allActsDifference{s} = actsVecDiffCat';
    toR.allActsVec1{s} = [actsVec{:}]';
    toR.allActsVec2{s} = [actsVec2{:}]';
    toR.allGlobalSig{s} = globalSig';
end


%%
out.allSubs =  vertcat(toR.allSubs_h{:});
out.allBinaryConf =  vertcat(toR.allBinaryConf_h{:});
out.allClass =  vertcat(toR.allClass_h{:});
out.allUnsignedAct =  vertcat(toR.allUnsignedAct_h{:});
out.allSignedAct =  vertcat(toR.allSignedAct_h{:});
out.allRT =  vertcat(toR.allRT_h{:});
out.allCresp = vertcat(toR.allCresp{:});
out.allCorResp = vertcat(toR.allCorResp{:});
out.allActsTrueClass = vertcat(toR.allActsTrueClass{:});
out.allActsFalseClass = vertcat(toR.allActsFalseClass{:});
out.allActDiffs = vertcat(toR.allActsDifference{:});
out.allActs1 = vertcat(toR.allActsVec1{:});
out.allActs2 = vertcat(toR.allActsVec2{:});
out.allGlobalSig = vertcat(toR.allGlobalSig{:});

toPrint(:,1) = ['subs'; num2cell(out.allSubs)];
toPrint(:,2) = ['binaryConf'; num2cell(out.allBinaryConf)];
toPrint(:,3)= ['class'; num2cell(out.allClass)];
toPrint(:,4)= ['unsignedAct'; num2cell(out.allUnsignedAct)];
toPrint(:,5)= ['signedAct'; num2cell(out.allSignedAct)];
toPrint(:,6)= ['RT'; num2cell(out.allRT)];
toPrint(:,7)= ['unsignedConf'; num2cell(out.allCresp)];
toPrint(:,8)= ['corResp'; num2cell(out.allCorResp)];
toPrint(:,10)= ['actsTrueClass'; num2cell(out.allActsTrueClass)];
toPrint(:,11)= ['actsFalseClass'; num2cell(out.allActsFalseClass)];
toPrint(:,12)= ['actsDiff'; num2cell(out.allActDiffs)];
toPrint(:,13)= ['acts1'; num2cell(out.allActs1)];
toPrint(:,14)= ['acts2'; num2cell(out.allActs2)];
%toPrint(:,15) = ['globalSig'; num2cell(out.allGlobalSig)];

cell2csv('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/ROutput/MnemonicData.csv', toPrint, ',', 2000);



end



% testing percentage of high/low conf in the classifier 
% for i=1:length(y.res.subj)
%     
%     onsInClassifier = y.res.subj{i}.penalty.nVox.weights.S.idxOnsets_train_in_classifier;
%     for j=1:4
%         thisConf = y.res.subj{i}.penalty.nVox.weights.S.idxTr.cresp==j;
%     
%     pctInClassifier(i,j) = sum(thisConf .* onsInClassifier)/sum(onsInClassifier);
%     end
% end


