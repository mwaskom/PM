function PM_MakeRegs(par)


cd (par.behavdir);

dFN = dir('*Perc*');
fileNames = {dFN.name}';

countdown = 12; %length of countdown period (in secs) at the beginning of each run 

if strcmp(par.task, 'mnem')
    [~, ~, idx] = Mnemonic_fMRIBehAnalysis_Retrieval(par);
elseif strcmp(par.task, 'perc')
    [~, ~, idx] = Perceptual_fMRIBehAnalysis(par);
end

qqq = load(par.classmat);
%[res dat] = PM_classificationPostProcessing(qqq.res, par.task);
thisSub = find(par.subNo == qqq.res.subjArray);

perf_set = {'cor' 'inc'};
perf_set_perc = {'corWithZeros', 'inc'};
class_set = {'face' 'house'};
resp_set = {'resp_face' 'resp_house'};
hand_set = {'respLeft' 'respRight'};

i = 0;

%idx.inClassifier = res.sub(thisSub).S.idxOnsets_test_in_classifier;
%idx.ev = dat{thisSub}.cv.reordered.unsignedEv;

%%
if strcmp(par.thisAnalysis, 'AnalysisRetByRREv')
    for l=1:length(class_set)
        for p=1:length(perf_set)
            
            idx.thisOns = logical(idx.(perf_set{p}) .* idx.(class_set{l}) .* idx.inClassifier);
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', class_set{l}, perf_set{p});
                
                stimOnsets{i}= idx.alltrials(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                idxThisOnsInClassifierOnsSpace = idx.thisOns(idx.inClassifier);
                thisEv = idx.ev(idxThisOnsInClassifierOnsSpace);
                
                pmod(i).name{1} = 'Ev';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = thisEv;
            end
        end
    end
%%
elseif strcmp(par.thisAnalysis, 'AnalysisRetByRREv_signedByHand')
    for l=1:length(class_set)
        for p=1:length(perf_set)
            
            idx.thisOns = logical(idx.(perf_set{p}) .* idx.(hand_set{l}) .* idx.inClassifier);
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', hand_set{l}, perf_set{p});
                
                stimOnsets{i}= idx.alltrials(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                idxThisOnsInClassifierOnsSpace = idx.thisOns(idx.inClassifier);
                thisEv = idx.ev(idxThisOnsInClassifierOnsSpace);
                
                pmod(i).name{1} = 'Ev';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = thisEv;
            end
        end
    end
    %%
elseif strcmp(par.thisAnalysis, 'analysisRetByConf')
    for l=1:length(class_set)
        for p=1:length(perf_set)
            
            idx.thisOns = logical(idx.(perf_set{p}) .* idx.(class_set{l}));
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', class_set{l}, perf_set{p});
                
                stimOnsets{i}= idx.alltrials(idx.thisOns);
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                pmod(i).name{1} = 'Conf';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = idx.unsignedConf(idx.thisOns);
            end
        end
    end
    %%
elseif strcmp(par.thisAnalysis, 'analysisRetByConf_signedByHand')
    for l=1:length(class_set)
        for p=1:length(perf_set)
            
            idx.thisOns = logical(idx.(perf_set{p}) .* idx.(hand_set{l}));
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', hand_set{l}, perf_set{p});
                
                stimOnsets{i}= idx.alltrials(idx.thisOns);
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                pmod(i).name{1} = 'Conf';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = idx.unsignedConf(idx.thisOns);
            end
        end
    end
    %%
elseif strcmp(par.thisAnalysis, 'analysisRetByRT')
    for l=1:length(class_set)
        for p=1:length(perf_set)
            
            idx.thisOns = logical(idx.(perf_set{p}) .* idx.(class_set{l}) .* ~isnan(idx.rt));
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', class_set{l}, perf_set{p});
                
                stimOnsets{i}= idx.alltrials(idx.thisOns);
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                pmod(i).name{1} = 'RT';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = idx.rt(idx.thisOns);
            end
        end
    end
    %%
elseif strcmp(par.thisAnalysis, 'analysisRetByRT_signedByHand')
    for l=1:length(class_set)
        for p=1:length(perf_set)
            
            idx.thisOns = logical(idx.(perf_set{p}) .* idx.(hand_set{l}) .* ~isnan(idx.rt));
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', hand_set{l}, perf_set{p});
                
                stimOnsets{i}= idx.alltrials(idx.thisOns);
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                pmod(i).name{1} = 'RT';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = idx.rt(idx.thisOns);
            end
        end
    end
elseif strcmp(par.thisAnalysis, 'analysisPercByConf')
    for l=1:length(class_set)
        for p=1:length(perf_set_perc)
            
            idx.thisOns = logical(idx.(perf_set_perc{p}) .* idx.(class_set{l}) .* idx.validTrials);
            
            thisPMod = idx.unsignedConf(find(idx.thisOns));
            isConstantPMod = length(unique(thisPMod))==1;
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', class_set{l}, perf_set_perc{p});
                
                stimOnsets{i}= idx.alltrials(idx.thisOns);
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                if ~isConstantPMod
                    pmod(i).name{1} = 'Conf';
                    pmod(i).poly{1} = 1;
                    pmod(i).param{1} = idx.unsignedConf(idx.thisOns);
                end
            end
        end
    end
    %%
elseif strcmp(par.thisAnalysis, 'analysisPercByConf_signedByHand')
    for l=1:length(class_set)
        for p=1:length(perf_set_perc)
            
            idx.thisOns = logical(idx.(perf_set_perc{p}) .* idx.(hand_set{l}) .* idx.validTrials);
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', hand_set{l}, perf_set_perc{p});
                
                stimOnsets{i}= idx.alltrials(idx.thisOns);
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                pmod(i).name{1} = 'Conf';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = idx.unsignedConf(idx.thisOns);
            end
        end
    end
elseif strcmp(par.thisAnalysis, 'analysisPercByRT')
    for l=1:length(class_set)
        for p=1:length(perf_set_perc)
            
            idx.thisOns = logical(idx.(perf_set_perc{p}) .* idx.(class_set{l}) .* ~isnan(idx.rt) .* idx.validTrials);
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', class_set{l}, perf_set_perc{p});
                
                stimOnsets{i}= idx.alltrials(idx.thisOns);
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                pmod(i).name{1} = 'RT';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = idx.rt(idx.thisOns);
            end
        end
    end
    %%
elseif strcmp(par.thisAnalysis, 'analysisPercByRT_signedByHand')
    for l=1:length(class_set)
        for p=1:length(perf_set_perc)
            
            idx.thisOns = logical(idx.(perf_set_perc{p}) .* idx.(hand_set{l}) .* ~isnan(idx.rt) .* idx.validTrials);
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', hand_set{l}, perf_set_perc{p});
                
                stimOnsets{i}= idx.alltrials(idx.thisOns);
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                pmod(i).name{1} = 'RT';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = idx.rt(idx.thisOns);
            end
        end
    end
elseif strcmp(par.thisAnalysis, 'AnalysisPercByPercEv')
    for l=1:length(class_set)
        for p=1:length(perf_set)
            
            idx.thisOns = logical(idx.(perf_set{p}) .* idx.(class_set{l}) .* idx.inClassifier);
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', class_set{l}, perf_set{p});
                
                stimOnsets{i}= idx.alltrials(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                idxThisOnsInClassifierOnsSpace = idx.thisOns(idx.inClassifier);
                thisEv = idx.ev(idxThisOnsInClassifierOnsSpace);
                
                pmod(i).name{1} = 'Ev';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = thisEv;
            end
        end
    end
    %%
elseif strcmp(par.thisAnalysis, 'AnalysisPercByPercEv_signedByHand')
    for l=1:length(class_set)
        for p=1:length(perf_set)
            
            idx.thisOns = logical(idx.(perf_set{p}) .* idx.(hand_set{l}) .* idx.inClassifier);
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', hand_set{l}, perf_set{p});
                
                stimOnsets{i}= idx.alltrials(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                idxThisOnsInClassifierOnsSpace = idx.thisOns(idx.inClassifier);
                thisEv = idx.ev(idxThisOnsInClassifierOnsSpace);
                
                pmod(i).name{1} = 'Ev';
                pmod(i).poly{1} = 1;
                pmod(i).param{1} = thisEv;
            end
        end
    end
elseif strcmp(par.thisAnalysis, 'analysisPercByAllConf')
    for l=1:length(class_set)
        
        
        idx.thisOns = logical(idx.(class_set{l}) .* idx.validTrials);
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_%s', class_set{l});
            
            stimOnsets{i}= idx.alltrials(idx.thisOns);
            stimNames{i} = fName;
            stimDurations{i} = 0;
            
            pmod(i).name{1} = 'Conf';
            pmod(i).poly{1} = 1;
            pmod(i).param{1} = idx.unsignedConf(idx.thisOns);
            
            pmod(i).name{2} = 'Acc';
            pmod(i).poly{2} = 1;
            pmod(i).param{2} = idx.corWithZeros(idx.thisOns);
        end
        
    end
elseif strcmp(par.thisAnalysis, 'analysisPercByAllRT')
    for l=1:length(class_set)
        
        
        idx.thisOns = logical(idx.(class_set{l}) .* idx.validTrials .* ~isnan(idx.rt));
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_', class_set{l});
            
            stimOnsets{i}= idx.alltrials(idx.thisOns);
            stimNames{i} = fName;
            stimDurations{i} = 0;
            
            pmod(i).name{1} = 'RT';
            pmod(i).poly{1} = 1;
            pmod(i).param{1} = idx.rt(idx.thisOns);
            
            pmod(i).name{2} = 'Acc';
            pmod(i).poly{2} = 1;
            pmod(i).param{2} = idx.corWithZeros(idx.thisOns);
        end
        
    end
elseif strcmp(par.thisAnalysis, 'analysisRetByAllConf')
    for l=1:length(class_set)
          
        idx.thisOns = logical(idx.(class_set{l})  );
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_', class_set{l});
            
            stimOnsets{i}= idx.alltrials(idx.thisOns);
            stimNames{i} = fName;
            stimDurations{i} = 0;
            
            pmod(i).name{1} = 'Conf';
            pmod(i).poly{1} = 1;
            pmod(i).param{1} = idx.unsignedConf(idx.thisOns);
            
            pmod(i).name{2} = 'Acc';
            pmod(i).poly{2} = 1;
            pmod(i).param{2} = idx.cor(idx.thisOns);
        end
    end
    
elseif strcmp(par.thisAnalysis, 'analysisRetByAllRT')
    for l=1:length(class_set)
        
        
        idx.thisOns = logical(idx.(class_set{l}) .*  ~isnan(idx.rt));
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_', class_set{l});
            
            stimOnsets{i}= idx.alltrials(idx.thisOns);
            stimNames{i} = fName;
            stimDurations{i} = 0;
            
            pmod(i).name{1} = 'RT';
            pmod(i).poly{1} = 1;
            pmod(i).param{1} = idx.rt(idx.thisOns);
            
            pmod(i).name{2} = 'Acc';
            pmod(i).poly{2} = 1;
            pmod(i).param{2} = idx.cor(idx.thisOns);
        end
        
    end
    
    %%
else
    error('no regressor scheme establibshed for this analysis')
end
%%
otherTrials = setdiff(idx.alltrials, [stimOnsets{:}]);

if ~isempty(otherTrials)
    i = i+1;
    stimOnsets{i} = otherTrials;
    stimNames{i} = 'junk';
    stimDurations{i} = 0;
end

onsets = stimOnsets;
names = stimNames;
durations = stimDurations;


if ~exist(par.analysisdir)
    mkdir (par.analysisdir);
end

sessReg = zeros(sum(par.numvols),length(par.numvols)-1);
for i = 1:(length(par.numvols)-1)
    if (i==1)
        sessReg(1:par.numvols(i),i) = ones(par.numvols(i),1);
    else
        sessIdx = (1+sum(par.numvols(1:(i-1)))):sum(par.numvols(1:i));
        sessReg(sessIdx,i) = ones(par.numvols(i),1);
    end
end

R = horzcat(sessReg);

cd (par.analysisdir);

save ons.mat onsets durations names pmod;
save regs.mat R