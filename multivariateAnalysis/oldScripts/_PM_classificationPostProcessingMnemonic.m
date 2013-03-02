function [res, groupPsy] = PM_classificationPostProcessingMnemonic(qqq)

for s = 1:length(qqq.subj)

subj_id = qqq.subj{s}.penalty.nVox.weights.iter{1}.header.subj_id;
par = PM_mvpa_params(subj_id, 'mnem');
%Design a function that ports the data to R.
%par2 = PM_params(subj_id);

%ret = Mnemonic_fMRIBehAnalysis_Retrieval(par2);

mvpa_ons = load(fullfile(par.onsetsTestDir, 'mvpa_ons'));
[~, testRegsIdx] = ismember([par.condsTest{:}], mvpa_ons.names);

[allOns_h aOI] = sort(vertcat(mvpa_ons.onsets{testRegsIdx}));

%idx.enoughTRs_h = vertcat(qqq.subj{s}.S.enoughTRs.test{:});
%idx.enoughTRs = idx.enoughTRs_h(aOI);

%allOns = allOns_h(find(idx.enoughTRs));
allOns = allOns_h;


for i = 1:length(qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations)
    % for each classification iteration, the vector of whether the
    % classification was correct or not
    correctsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).perfmet.corrects; 
    actsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).acts(1,:);
    desiredsVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).perfmet.desireds; 
    guessesVec{i} = qqq.subj{s}.penalty.nVox.weights.iter{1}.iterations(i).perfmet.guesses; 
    signedActs{i} = abs(actsVec{i} - .5) .* (2*correctsVec{i} - 1);
end

correctsVecCat = [correctsVec{:}];
actsVecCat = [actsVec{:}];
desiredsVecCat = [desiredsVec{:}];
absActsVecCat = abs(actsVecCat-.5);
guessesVecCat = [guessesVec{:}];
signedActsVecCat = [signedActs{:}];

classPat = {'face' 'house'};
class = {'face' 'house'};
confLevel = {'high' 'low'};

idx.cor = ~cellfun('isempty', strfind(mvpa_ons.names,'cor'));
idx.inc = ~cellfun('isempty', strfind(mvpa_ons.names,'inc'));
idx.high = ~cellfun('isempty', strfind(mvpa_ons.names,'high'));
idx.low = ~cellfun('isempty', strfind(mvpa_ons.names,'low'));
idx.face = ~cellfun('isempty', strfind(mvpa_ons.names,'face'));
idx.house = ~cellfun('isempty', strfind(mvpa_ons.names,'house'));

% corByCohRegVals = find(idx.CohReg .* ~idx.inc) ;
% confRegVals = find(idx.conf);

j = 0;
for i = 1:length(class)
    for c = 1:length(confLevel)
        
        idx.thisClass = ~cellfun('isempty', strfind(mvpa_ons.names, classPat{i}));
        idx.thisConfLevel = ~cellfun('isempty', strfind(mvpa_ons.names, confLevel{c} ));
        
         if sum(idx.thisClass .* idx.thisConfLevel .* ~idx.inc) > 0;
             
             j = j+1;
             idx.thisIt{j} = idx.thisClass .* idx.thisConfLevel .* ~idx.inc;
             
             mo_idx(j) = find(idx.thisIt{j});
             res.names{j} = [class{i} '_' confLevel{c} ] ;
         end
    end
end



orderedLabs = zeros(size(allOns));
for i = 1:length(mo_idx)
    %classCoh_h{i}= i*ones(1,length(mvpa_ons.onsets{mo_idx(i)})); 
    
    [~, ix1] = ismember(mvpa_ons.onsets{mo_idx(i)}, allOns);
    
    ixNoZeros = ix1(ix1>0);
    orderedLabs(ixNoZeros) = i;
end


for i = 1:length(mo_idx)
regIX = orderedLabs==i;
regCor = correctsVecCat(find(regIX));
res.coh.classperf(i,s) = mean(regCor);

res.N.classperf(i,s) = length(regCor);

regAct = actsVecCat(find(regIX));
res.coh.medianActs(i,s) = nanmedian(regAct);

res.coh.pGuessFace(i,s) = nanmean(guessesVecCat(find(regIX))==1);

end

groupPsy(1).pFaceGuessByConf.mean = nanmean(res.coh.pGuessFace([3 4 2 1],:)');
groupPsy(1).pFaceGuessByConf.SE = ste(res.coh.pGuessFace([3 4 2 1], :)');

groupPsy(1).pFaceGuessByConf.ticks.x = 0:5;
groupPsy(1).pFaceGuessByConf.ticks.y = [0 1];
groupPsy(1).pFaceGuessByConf.xlabel = 'Confidence (high house to high face)';
groupPsy(1).pFaceGuessByConf.ylabel = 'pFaceGuess' ;
groupPsy(1).pFaceGuessByConf.xVals = 1:4;
groupPsy(1).pFaceGuessByConf.marker = 'o-';


groupPsy(1).medActByConf.mean = nanmean(res.coh.medianActs([3 4 2 1],:)');
groupPsy(1).medActByConf.SE = ste(res.coh.medianActs([3 4 2 1], :)');

groupPsy(1).medActByConf.ticks.x = 0:5;
groupPsy(1).medActByConf.ticks.y = [.4 .6];
groupPsy(1).medActByConf.xlabel = 'Confidence (high house to high face)';
groupPsy(1).medActByConf.ylabel = 'Classifier Output' ;
groupPsy(1).medActByConf.xVals = 1:4;
groupPsy(1).medActByConf.marker = 'o-';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
orderedBinaryConf = 3 - ceil(abs(orderedLabs-2.5)); % 1 for hi conf, 2 for lo

for c = 1:length(confLevel)
    idx.thisConfLevel = orderedBinaryConf == c;
    
    res.cor.conf(c,s) = nanmean(correctsVecCat(idx.thisConfLevel==1));
    res.medSignedActs.conf(c,s) = nanmedian(signedActsVecCat(idx.thisConfLevel==1));
    res.medActs.conf(c,s) = nanmedian(absActsVecCat(idx.thisConfLevel==1));
end


res.cor.class1(s) = nanmean(correctsVecCat(desiredsVecCat==1));
res.cor.class2(s) = nanmean(correctsVecCat(desiredsVecCat==2));

res.medianAbsActs(s) = nanmean(signedActsVecCat);

end

%%%%



sprintf('end');







