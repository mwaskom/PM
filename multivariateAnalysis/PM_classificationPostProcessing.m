function [res groupPsy idx idxB] = PM_classificationPostProcessing(qqq)


res.groupCatDat.acts = [];
res.groupCatDat.cor = [];
res.groupCatDat.desireds = [];
res.groupCatDat.guesses = [];
res.groupCatDat.coh_signed = [];
res.groupCatDat.conf =[];
res.groupCatDat.sub = [];
res.groupCatDat.RT = [];


for s = 1:length(qqq.subjArray)
    subj_id = qqq.subjArray(s);
    par = PM_Params(subj_id, 'perc');
    
    %par = PM_Params(qqq.subj{s}.S.subj_array{s});
    [~, ~, idxB, gRes] = Perceptual_fMRIBehAnalysis(par);
    
    res.gRes(s) = gRes;
    
    mvpa_ons = load(fullfile(qqq.subj{s}.penalty.nVox.weights.S.onsetsTestDir, 'mvpa_ons'));
    
    %Design a function that ports the data to R.
    
    allOns = sort(unique(vertcat(mvpa_ons.onsets{:})));
    
    for i = 1:length(qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations)
        % for each classification iteration, the vector of whether the
        % classification was correct or not
        correctsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).perfmet.corrects;
        actsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).acts(1,:);
        actsVec2{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).acts(2,:);
        desiredsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).perfmet.desireds;
        guessesVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).perfmet.guesses;
    end
    
    
    desiredsVecCat = [desiredsVec{:}];
    
    
    
    if length(qqq.subj{s}.penalty.nVox.weights.S.condnames)>2
        actsVecCat_h = [actsVec{:}];
        actsVecCat2_h = [actsVec2{:}];
        actsVecCat = actsVecCat_h ./ (actsVecCat_h + actsVecCat2_h );
        
        guessesVecCat = 2 - (actsVecCat>.5);
        correctsVecCat = guessesVecCat==desiredsVecCat;
    else
        actsVecCat = [actsVec{:}];
        correctsVecCat = [correctsVec{:}];
        guessesVecCat = [guessesVec{:}];
    end
    
    absActsVecCat = abs(actsVecCat - .5);
    actsVecCatLogit = log(actsVecCat ./ (1- actsVecCat));

    actsVecLogitUnsigned = actsVecCatLogit;
    actsVecLogitUnsigned(desiredsVecCat==2) = -1 * actsVecLogitUnsigned(desiredsVecCat==2);
    
    cohs = [20 35 45 60 100];
    %cohs = [0 20 35 60 100];
    
    classPat = {'face' 'house'};
    class = {'face' 'house'};
    conf = 1:4;
%     confSigned = {'house_conf4_cor' 'house_conf3_cor' 'house_conf2_cor' 'house_conf1_cor' ...
%         'face_conf1_cor' 'face_conf3_cor' 'face_conf2_cor' 'face_conf1_cor'};
    
    idx.CorReg = ~cellfun('isempty', strfind(mvpa_ons.names,'cor'));
    idx.CohReg = ~cellfun('isempty', strfind(mvpa_ons.names,'coh'));
    idx.coh0 = ~cellfun('isempty', strfind(mvpa_ons.names,'coh0'));
    idx.inc = ~cellfun('isempty', strfind(mvpa_ons.names,'inc'));
    idx.unsignedActs = actsVecLogitUnsigned;
    
    idx.conf = ~cellfun('isempty', strfind(mvpa_ons.names,'conf'));
    
    corByCohRegVals = find(idx.CohReg .* ~idx.inc) ;
    confRegVals = find(idx.conf);
    
    j = 0;
    for i = 1:length(class)
        for c = 1:length(cohs)   
            
            idx.thisClass = ~cellfun('isempty', strfind(mvpa_ons.names, classPat{i}));
            idx.thisCoh = ~cellfun('isempty', strfind(mvpa_ons.names, ['coh' num2str(cohs(c))] ));
            
            if sum(idx.thisClass .* idx.thisCoh .* ~idx.inc) > 0;
                
                j = j+1;
                idx.thisIt{j} = idx.thisClass .* idx.thisCoh .* ~idx.inc;
                
                mo_idx{j} = find(idx.thisIt{j});
                res.sub(s).names{j} = [class{i} '_' num2str(cohs(c))] ;
                %res.sub(s).coh_signed(j) = ((cohs(c)/1.00001 +
                %.00000000001) * (2*strcmp(classPat{i}, 'face') - 1))/100;
                % to separate 0% chosen face and house
                res.sub(s).coh_signed(j) = ((cohs(c)) * (2*strcmp(classPat{i}, 'face') - 1))/100;
            end
            
        end
    end
    idx.thisCoh = ~cellfun('isempty', strfind(mvpa_ons.names, ['coh0']));
    j = j+1;
    idx.thisIt{j} = idx.thisCoh;
    mo_idx{j} = find(idx.thisIt{j});
    res.sub(s).names{j} = ['coh0'] ;
    res.sub(s).coh_signed(j) = 0;
    
    
    orderedLabs = zeros(size(allOns));
    for i = 1:length(mo_idx)
        
        [~, ix1] = ismember(vertcat(mvpa_ons.onsets{mo_idx{i}}), allOns);
        orderedLabs(ix1) = i;
    end
    
    if qqq.subj{1}.penalty.nVox.weights.S.linReg
        [res.sub(s).coh.corrPredictedAndActualCoh.r  res.sub(s).corrPredictedAndActualCoh.p] = ...
            corr(desiredsVecCat', actsVecCatLogit');
        
%         [res.sub(s).coh.corrPredictedAndActualCoh_corOnly.r  res.sub(s).corrPredictedAndActualCoh_corOnly.p] = ...
%             corr(desiredsVecCat(find(idxB.cor))', actsVecCat(find(idxB.cor))');
%         
%         [res.sub(s).coh.corrPredictedAndActualCoh_faceCor.r  res.sub(s).corrPredictedAndActualCoh_faceCor.p] = ...
%             corr(desiredsVecCat(find(idxB.cor .* idxB.face))', actsVecCat(find(idxB.cor .* idxB.face))');
%         
%         [res.sub(s).coh.corrPredictedAndActualCoh_houseCor.r  res.sub(s).corrPredictedAndActualCoh_houseCor.p] = ...
%             corr(desiredsVecCat(find(idxB.cor .* idxB.house))', actsVecCat(find(idxB.cor .* idxB.house))');
        
        res.sub(s).coh.pctClassCorFace = sum((desiredsVecCat>0) .* (actsVecCatLogit>0))/sum(desiredsVecCat>0);
        res.sub(s).coh.pctClassCorHouse = sum((desiredsVecCat<0) .* (actsVecCatLogit<0))/sum(desiredsVecCat<0);
        
    else
    
    for i = 1:length(mo_idx)
        regIX = orderedLabs==i;
        regCor = correctsVecCat(find(regIX));
        res.sub(s).coh.classperf(i) = mean(regCor);
        res.sub(s).coh.N(i) = length(regCor);
        
        res.sub(s).coh.pFaceGuess(i) = mean(actsVecCatLogit(regIX) > .5);
        
        regActCor = actsVecCatLogit(find(regIX .* correctsVecCat'));
        regActAll = actsVecCatLogit(find(regIX ));
        
        res.sub(s).coh.meanActsCor(i) = mean(regActCor);
        res.sub(s).coh.meanActsAll(i) = mean(regActAll);
        res.sub(s).coh.seActsAll(i) = ste(regActAll);
        
    end
    
    
    
    [cs_sorted, ix_cs] = sort(res.sub(s).coh_signed);
    
    
    j = 0;
    
    for i = ix_cs
        
        j = j + 1;
        groupPsy(s).dat.pCorByCoh.mean(j) = res.sub(s).coh.classperf(i);
        groupPsy(s).dat.pCorByCoh.SE(j) = sqrt((res.sub(s).coh.classperf(i) * (1-res.sub(s).coh.classperf(i))) / res.sub(s).coh.N(i) );
        
        groupPsy(s).dat.actsByCoh.mean(j) = res.sub(s).coh.meanActsAll(i);
        groupPsy(s).dat.actsByCoh.SE(j) = res.sub(s).coh.seActsAll(i);
        
        groupPsy(s).dat.pFaceGuessByCoh.mean(j) = res.sub(s).coh.pFaceGuess(i);
        groupPsy(s).dat.pFaceGuessByCoh.SE(j) = sqrt((res.sub(s).coh.pFaceGuess(i) * (1-res.sub(s).coh.pFaceGuess(i))) / res.sub(s).coh.N(i) );
    end
    
    groupPsy(s).dat.pCorByCoh.ticks.x =  [-1:.1:1];
    groupPsy(s).dat.pCorByCoh.ticks.y = [0 1];
    groupPsy(s).dat.pCorByCoh.xlabel = 'Coherence (high house to high face)';
    groupPsy(s).dat.pCorByCoh.ylabel = 'Percent Correct' ;
    groupPsy(s).dat.pCorByCoh.xVals = cs_sorted;
    groupPsy(s).dat.pCorByCoh.marker = 'o-';
    
    groupPsy(s).dat.actsByCoh.ticks.x = [-1:.1:1];
    groupPsy(s).dat.actsByCoh.ticks.y = [-8 8];
    groupPsy(s).dat.actsByCoh.xlabel = 'Coherence (high house to high face)';
    groupPsy(s).dat.actsByCoh.ylabel = 'Classifier Output' ;
    groupPsy(s).dat.actsByCoh.xVals = cs_sorted;
    groupPsy(s).dat.actsByCoh.marker = 'o-';
    
    groupPsy(s).dat.pFaceGuessByCoh.ticks.x = [-1:.1:1];
    groupPsy(s).dat.pFaceGuessByCoh.ticks.y = [0 1];
    groupPsy(s).dat.pFaceGuessByCoh.xlabel = 'Coherence (high house to high face)';
    groupPsy(s).dat.pFaceGuessByCoh.ylabel = 'Probability of a Face Guess' ;
    groupPsy(s).dat.pFaceGuessByCoh.xVals = cs_sorted;
    groupPsy(s).dat.pFaceGuessByCoh.marker = 'o';

    corCohSigned = idxB.coh_signed(idxB.cor==1);
    guessFaceCor = actsVecCatLogit(idxB.cor==1) > .5; %of correctly retrieved trials, which were given 'face' guesses?
    
    g_coh_signed = -1:0.002:1;
    groupPsy(s).dat.pFaceGuessByCoh.g_coh_signed = g_coh_signed;
    beta = logistfit([ones(length(corCohSigned),1) corCohSigned guessFaceCor']);
    groupPsy(s).dat.pFaceGuessByCoh.g_fc = 1./(1+exp(-(beta(1)+beta(2)*g_coh_signed)));
    
    groupPsy(s).goodSub = true;
    res.sub(s).pFaceGuessByCoh.beta = beta;
    %% confidence
    conf_ons = load(fullfile(par.subdir, 'analysis_percDMByConf_3d', 'ons'));
    
    for f=conf;
        idx.thisConf = ~cellfun('isempty', strfind(conf_ons.names, sprintf('conf%g_cor', f)));
        theseOns{f} = vertcat(conf_ons.onsets{find(idx.thisConf)});
        conf_h{f}= f*ones(1,length(theseOns{f}));
    end
    
    confCat = [conf_h{:}];
    
    z2 = vertcat(theseOns{:});
    [~, ix2] = sort(z2);
    
    conf_sorted = confCat(ix2);
    
    res.sub(s).conf.conf = conf_sorted;
    
    for i = 1:length(conf_h)
        regIX = conf_sorted==i;
        regCor = correctsVecCat(find(regIX));
        res.sub(s).conf.classperf(i) = mean(regCor);
        
        regAct = actsVecCatLogit(find(regIX));
        res.sub(s).conf.meanActs(i) = mean(regAct);
    end
    
    
    
   
    
    for i = 1:8
        regIX = ((idxB.conf'==i) .* ~idxB.inc);
        regCor = correctsVecCat(find(regIX));
        res.sub(s).signedconf.classperf(i) = mean(regCor);
        
        regAct = actsVecCatLogit(find(regIX));
        res.sub(s).signedconf.meanActs(i) = mean(regAct);
        
        res.sub(s).signedconf.pFaceGuess(i) = mean(actsVecCatLogit(find(regIX)) > .5);
        res.sub(s).signedconf.N(i) = sum(regIX);
        
        groupPsy(s).dat.pCorByConf.mean(i) = res.sub(s).signedconf.classperf(i);
        groupPsy(s).dat.pCorByConf.SE(i) = sqrt((mean(regCor) * (1-mean(regCor))) / length(regCor));
        
        groupPsy(s).dat.actsByConf.mean(i) = res.sub(s).signedconf.meanActs(i);
        groupPsy(s).dat.actsByConf.SE(i) = ste(regAct);
        
        
        groupPsy(s).dat.pFaceGuessByConf.mean(i) = res.sub(s).signedconf.pFaceGuess(i);
        groupPsy(s).dat.pFaceGuessByConf.SE(i) = sqrt((res.sub(s).signedconf.pFaceGuess(i) * (1-res.sub(s).signedconf.pFaceGuess(i))) / res.sub(s).signedconf.N(i) );
    end
    
    
    groupPsy(s).dat.pCorByConf.ticks.x = 0:9;
    groupPsy(s).dat.pCorByConf.ticks.y = [0 1];
    groupPsy(s).dat.pCorByConf.xlabel = 'Confidence (high house to high face)';
    groupPsy(s).dat.pCorByConf.ylabel = 'Percent Correct' ;
    groupPsy(s).dat.pCorByConf.xVals = 1:8;
    groupPsy(s).dat.pCorByConf.marker = 'o-';
    
    
    groupPsy(s).dat.actsByConf.ticks.x = 0:9;
    groupPsy(s).dat.actsByConf.ticks.y = [-8 8];
    groupPsy(s).dat.actsByConf.xlabel = 'Confidence (high house to high face)';
    groupPsy(s).dat.actsByConf.ylabel = 'Classifier Output' ;
    groupPsy(s).dat.actsByConf.xVals = 1:8;
    groupPsy(s).dat.actsByConf.marker = 'o-';
    
    
    groupPsy(s).dat.pFaceGuessByConf.ticks.x = 0:9;
    groupPsy(s).dat.pFaceGuessByConf.ticks.y = [0 1];
    groupPsy(s).dat.pFaceGuessByConf.xlabel = 'Confidence (high house to high face)';
    groupPsy(s).dat.pFaceGuessByConf.ylabel = 'Probability of a Face Guess' ;
    groupPsy(s).dat.pFaceGuessByConf.xVals = 1:8;
    groupPsy(s).dat.pFaceGuessByConf.marker = 'o';
    
    corConfSigned = idxB.conf(idxB.cor==1)';
    guessFaceCor = actsVecCatLogit(idxB.cor==1)' > .5; %of correctly retrieved trials, which were given 'face' guesses?
    
    g_coh_signed = 1:0.002:8;
    groupPsy(s).pFaceGuessByConf.g_coh_signed = g_coh_signed;
    beta = logistfit([ones(length(corConfSigned),1) corConfSigned guessFaceCor]);
    groupPsy(s).pFaceGuessByConf.g_fc = 1./(1+exp(-(beta(1)+beta(2)*g_coh_signed)));
    res.sub(s).pFaceGuessByConf.beta = beta;
    

    %% evidence and rt
    
    %figure; 
   % scatter(actsVecCatLogit(idxB.cor==1), idxB.rt(idxB.cor==1));
    
    

    %% correct vs. incorrect retrieval
    res.corClass(s).corRet = mean(correctsVecCat(idxB.cor==1));
    res.corClass(s).incRet = mean(correctsVecCat(idxB.inc==1));
    
    res.corClass(s).corSE = sqrt(((res.corClass(s).corRet)*(1-res.corClass(s).corRet) )/ sum(idxB.cor));
    res.corClass(s).incSE = sqrt(((res.corClass(s).incRet)*(1-res.corClass(s).incRet) )/ sum(idxB.inc));
    
    res.actsClass(s).corRet = median(absActsVecCat(idxB.cor==1));
    res.actsClass(s).incRet = median(absActsVecCat(idxB.inc==1));
    
    res.actsClass(s).corSE = ste(absActsVecCat(idxB.cor==1));
    res.actsClass(s).incSE = ste(absActsVecCat(idxB.inc==1));
    
    
    res.groupCatDat.acts = vertcat(res.groupCatDat.acts, actsVecCatLogit');
    res.groupCatDat.cor = vertcat(res.groupCatDat.cor, idxB.cor);
    res.groupCatDat.desireds = vertcat(res.groupCatDat.desireds, desiredsVecCat');
    res.groupCatDat.guesses = vertcat(res.groupCatDat.guesses, guessesVecCat');
    res.groupCatDat.coh_signed = vertcat(res.groupCatDat.coh_signed, idxB.coh_signed);
    res.groupCatDat.conf = vertcat(res.groupCatDat.conf, idxB.conf');
    res.groupCatDat.sub = vertcat(res.groupCatDat.sub, s*ones(size(actsVecCatLogit')));
    res.groupCatDat.RT = vertcat(res.groupCatDat.RT, idxB.rt);
    
    %% res histograms
    
    res.hist.actsCorFaceProb{s} =  actsVecCatLogit((idxB.cor .*  idxB.face)==1);
    res.hist.actsCorHouseProb{s} = actsVecCatLogit((idxB.cor .*  idxB.house)==1);    
    
    res.hist.actsCorFaceLogit{s} =  actsVecCatLogit ((idxB.cor .*  idxB.face)==1);
    res.hist.actsCorHouseLogit{s} = actsVecCatLogit((idxB.cor .*  idxB.house)==1); 
    %% stats
    res.linreg(s) = regstats(actsVecCatLogit', [idxB.conf' idxB.coh_signed  idxB.cor desiredsVecCat']);
    
    [B dev res.logreg(s).stats]  = mnrfit([idxB.conf' idxB.coh_signed  idxB.cor desiredsVecCat'], (guessesVecCat'));
    
    
    [B dev res.sub(s).pFaceGuessByConf.stats] = mnrfit(corConfSigned', (guessFaceCor+1)');
    [B dev res.sub(s).pFaceGuessByCoh.stats]  = mnrfit(corCohSigned', (guessFaceCor+1)');
    
    res.sub(s).signedEvByCoh.stats = regstats(actsVecCatLogit(idxB.cor==1)'  , corCohSigned);
    res.sub(s).signedEvByConf.stats = regstats(actsVecCatLogit(idxB.cor==1)' , corConfSigned);
    
    %[~,~,~,res.sub(s).actsByPercPerf.stats] = ttest2(absActsVecCat(idxB.cor==1), absActsVecCat(idxB.inc==1));
    
    [table,chi2,res.sub(s).ClassPerfByPercPerf.p] = crosstab(correctsVecCat, idxB.cor);

    [~,~,~,res.sub(s).actsByPercPerf.stats] = ttest2(absActsVecCat(idxB.cor==1), absActsVecCat(idxB.inc==1));
    
    idx.legitRT = ~isnan(idxB.rt);
    idx.corLegitRT = (idx.legitRT .* idxB.cor)==1;
    res.sub(s).evByConfCorAndRT.stats = regstats(actsVecLogitUnsigned(idx.corLegitRT)' , [idxB.conf_unsigned(idx.corLegitRT) idxB.coh_(idx.corLegitRT) idxB.rt(idx.corLegitRT) ]);
    res.sub(s).ConfByEvAndCoh.stats = regstats(corConfSigned  , [actsVecCatLogit(idxB.cor==1)' corCohSigned ]);
    res.sub(s).evByConfAndCor_allPercTrials.stats = regstats(actsVecCatLogit' , [idxB.conf_signed idxB.coh_signed idxB.cor]);
    
    [B dev res.sub(s).pFaceClass_ByConfAndCor_allPercTrials.stats] = mnrfit([corCohSigned corConfSigned], (guessFaceCor+1)');
    
   
    [res.sub(s).rtCor.stats] = regstats(idxB.rt(idx.legitRT) , actsVecCatLogit(idx.legitRT)'  );
    [res.sub(s).rtCor.r,res.sub(s).rtCor.p] =  corr(idxB.rt((idx.legitRT .* idxB.cor)==1) , actsVecLogitUnsigned((idx.legitRT .* idxB.cor)==1)');
    [res.sub(s).rtCorFace.r,res.sub(s).rtCorFace.p] =  corr(idxB.rt((idx.legitRT .* idxB.cor .* idxB.face)==1) , actsVecLogitUnsigned((idx.legitRT .* idxB.cor .* idxB.face)==1)');
    [res.sub(s).rtCorHouse.r,res.sub(s).rtCorHouse.p] =  corr(idxB.rt((idx.legitRT .* idxB.cor .* idxB.house)==1) , actsVecLogitUnsigned((idx.legitRT .* idxB.cor .* idxB.house)==1)');
   
    end
end



    

