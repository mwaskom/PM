function [res, groupPsy, dat, datB] = PM_classificationPostProcessingPerceptualImproved(qqq)

groupPsy = [];
pnl = 1;
weights = 1;
useResiduals = false;
for s = 1:length(qqq.subj)
    
    subj_id = qqq.subjArray(s);
    %[par idxB] = PM_mvpa_params(subj_id, 'mnem');
    
    par = PM_Params(subj_id, 'perc', 0);

    [~, ~, idxB] = Perceptual_fMRIBehAnalysis(par);
    
    resStruct = qqq.subj{s}.penalty(pnl).nVox.weights(weights).iter{1}.iterations;
    resS = qqq.subj{s}.penalty(pnl).nVox.weights(weights).S;
    
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
    correctsVecCat = correctsVecCat(testIdxInv);
    guessesVecCat = guessesVecCat(testIdxInv);
    signedActsVecCat = signedActsVecCat(testIdxInv);
    desiredsVecCat = desiredsVecCat(testIdxInv);
    
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
    idx.onsInClassifier = ismember(floor(idxB.alltrials), floor(allOns));
    
    classPat = {'face' 'house'};
    class = {'face' 'house'};
    conf = {'1' '2' '3' '4'};
    
    idx.cor = idxB.cor(idx.onsInClassifier)';
    idx.inc = idxB.inc(idx.onsInClassifier)';
    idx.face = idxB.face(idx.onsInClassifier)';
    idx.house = idxB.house(idx.onsInClassifier)';
    idx.rt = idxB.rt(idx.onsInClassifier)';
    idx.conf_unsigned = idxB.conf_unsigned(idx.onsInClassifier)';
    idx.unsignedActs = actsVecLogitUnsigned;
    idx.corWithZeros = idxB.corWithZeros(idx.onsInClassifier)';
    
    dat{s} = idx;
    datB{s} = idxB; 
    
    % for organizing by signed binary confidence
    j = 0;
    for i = 1:length(class)
            j = j+1;
            thisClass = class{i};
            res.names{j} = [class{i}] ;
            if sum(idx.(thisClass) .* idx.cor) > 0;
                idx.thisIt{j} = idx.(thisClass) .* idx.cor;
                medianAct(j) = median(actsVecLogitUnsigned(idx.thisIt{j}==1));
            else
                medianAct(j) = NaN;
            end
    end
    
    res.cor.class1(s) = nanmean(correctsVecCat(desiredsVecCat==1));
    res.cor.class2(s) = nanmean(correctsVecCat(desiredsVecCat==2));
    res.cor.class1Cor(s) = nanmean(correctsVecCat((desiredsVecCat==1) & (idx.cor==1) ));
    res.cor.class2Cor(s) = nanmean(correctsVecCat((desiredsVecCat==2) & (idx.cor==1) ));
    res.cor.class1CorHighConf(s) = nanmean(correctsVecCat((desiredsVecCat==1) & (idx.cor==1) & (idx.conf_unsigned>2)));
    res.cor.class2CorHighConf(s) = nanmean(correctsVecCat((desiredsVecCat==2) & (idx.cor==1) & (idx.conf_unsigned>2)));
    res.MeanNTrainingTrials(s) = nanmean(nTrainingTrials);
    
    %     
%     orderingVec = [3 4 2 1]; %reorder the res variables to go from high confidence house to high confidenc face
%     res.coh.medianActs(:,s) = medianAct(orderingVec);
%     res.coh.medianActsSE(:,s) = .001*ones(size(medianAct(orderingVec)));
%     
    
    %     groupPsy(s).dat.pFaceGuessByConf.mean = res.coh.pGuessFace(orderingVec,s)';
    %     groupPsy(s).dat.pFaceGuessByConf.SE = res.coh.pGuessFaceSE(orderingVec,s);
    %     groupPsy(s).dat.pFaceGuessByConf.ticks.x = 0:5;
    %     groupPsy(s).dat.pFaceGuessByConf.ticks.y = [0 1];
    %     groupPsy(s).dat.pFaceGuessByConf.xlabel = 'Confidence (high house to high face)';
    %     groupPsy(s).dat.pFaceGuessByConf.ylabel = 'pFaceGuess' ;
    %     groupPsy(s).dat.pFaceGuessByConf.xVals = 1:4;
    %     groupPsy(s).dat.pFaceGuessByConf.marker = 'o-';
    
%     
%     groupPsy(s).dat.medActByConf.mean = res.coh.medianActs';
%     groupPsy(s).dat.medActByConf.SE = res.coh.medianActsSE';
%     groupPsy(s).dat.medActByConf.ticks.x = 0:5;
%     groupPsy(s).dat.medActByConf.ticks.y = [-1 5];
%     groupPsy(s).dat.medActByConf.xlabel = 'Confidence (high house to high face)';
%     groupPsy(s).dat.medActByConf.ylabel = 'Median Classifier Output' ;
%     groupPsy(s).dat.medActByConf.xVals = 1:4;
%     groupPsy(s).dat.medActByConf.marker = 'o-';
    
    
    %     groupPsy(s).dat.AccuracyByUnsignedConf.mean = res.coh.classperf(orderingVec,s)';
    %     groupPsy(s).dat.AccuracyByUnsignedConf.SE = res.coh.classperfSE(orderingVec,s);
    %     groupPsy(s).dat.AccuracyByUnsignedConf.ticks.x = 0:5;
    %     groupPsy(s).dat.AccuracyByUnsignedConf.ticks.y = [0 1];
    %     groupPsy(s).dat.AccuracyByUnsignedConf.xlabel = 'Confidence (high house to high face)';
    %     groupPsy(s).dat.AccuracyByUnsignedConf.ylabel = 'pCorrect' ;
    %     groupPsy(s).dat.AccuracyByUnsignedConf.xVals = 1:4;
    %     groupPsy(s).dat.AccuracyByUnsignedConf.marker = 'o-';
    
    groupPsy(s).goodSub = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    toR.allSubs_h{s} = s*ones(size(idx.cor))';
    toR.allClass_h{s} = double(idx.face');
    toR.allUnsignedAct_h{s} = actsVecLogitUnsigned';
    toR.allSignedAct_h{s} = actsVecCatLogit';
    toR.allRT_h{s} = idx.rt';
    toR.allConf_unsigned{s} = idx.conf_unsigned';
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
out.allClass =  vertcat(toR.allClass_h{:});
out.allUnsignedAct =  vertcat(toR.allUnsignedAct_h{:});
out.allSignedAct =  vertcat(toR.allSignedAct_h{:});
out.allRT =  vertcat(toR.allRT_h{:});
out.allConf_unsigned = vertcat(toR.allConf_unsigned{:});
out.allCorResp = vertcat(toR.allCorResp{:});
out.allActsTrueClass = vertcat(toR.allActsTrueClass{:});
out.allActsFalseClass = vertcat(toR.allActsFalseClass{:});
out.allActDiffs = vertcat(toR.allActsDifference{:});
out.allActs1 = vertcat(toR.allActsVec1{:});
out.allActs2 = vertcat(toR.allActsVec2{:});
out.allGlobalSig = vertcat(toR.allGlobalSig{:});

toPrint(:,1) = ['subs'; num2cell(out.allSubs)];
toPrint(:,2)= ['class'; num2cell(out.allClass)];
toPrint(:,3)= ['unsignedAct'; num2cell(out.allUnsignedAct)];
toPrint(:,4)= ['signedAct'; num2cell(out.allSignedAct)];
toPrint(:,5)= ['RT'; num2cell(out.allRT)];
toPrint(:,6)= ['unsignedConf'; num2cell(out.allConf_unsigned)];
toPrint(:,7)= ['corResp'; num2cell(out.allCorResp)];
toPrint(:,8)= ['actsTrueClass'; num2cell(out.allActsTrueClass)];
toPrint(:,9)= ['actsFalseClass'; num2cell(out.allActsFalseClass)];
toPrint(:,10)= ['actsDiff'; num2cell(out.allActDiffs)];
toPrint(:,11)= ['acts1'; num2cell(out.allActs1)];
toPrint(:,12)= ['acts2'; num2cell(out.allActs2)];
%toPrint(:,15) = ['globalSig'; num2cell(out.allGlobalSig)];

cell2csv('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/ROutput/PercData.csv', toPrint, ',', 2000);



end







