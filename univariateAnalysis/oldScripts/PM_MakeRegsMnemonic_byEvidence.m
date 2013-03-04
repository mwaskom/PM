function PM_MakeRegsMnemonic_byEvidence(par)


qqq_h = load(par.classmat);
qqq = qqq_h.res;
[~,~, dat, datB] = PM_classificationPostProcessingMnemonicImproved(qqq);

[~,subNo] = ismember(par.subNo, qqq.subjArray);

idx = dat{subNo};
idxB = datB{subNo};

perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
conf_set = {'low' 'high'};


i = 0;

for l=1:length(class_set)
    
    for p=1:length(perf_set)
        
        %for c = 1:length(conf_set)
            %idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l}) .* idx.(conf_set{c});
            idx.thisOns = idxB.(perf_set{p}) .* idxB.(class_set{l}) .* idx.onsInClassifier;
            
            if ~isempty(idxB.alltrials(find(idx.thisOns)))
                
                i = i+1;
                %fName = sprintf('%s_%s_%s', class_set{l},  perf_set{p},  conf_set{c});
                fName = sprintf('%s_%s', class_set{l},  perf_set{p});
                
                stimOnsets{i}= idxB.alltrials(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                idx.thisClassOns = idx.(perf_set{p}) .* idx.(class_set{l});
                pmod(i).name = {'unsignedClassOut'}; 
                pmod(i).param = {idx.unsignedActs(find(idx.thisClassOns))}; 
                pmod(i).poly = {1};
            end
            
        %end
    end
end

if any(~idx.onsInClassifier)
    stimOnsets{i+1} = idxB.alltrials(find(~idx.onsInClassifier));
    stimNames{i+1} = 'junk';
    stimDurations{i+1} = 0;
end

onsets = stimOnsets;
durations = stimDurations;
names = stimNames;

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

R = sessReg;

cd (par.analysisdir);
save ons.mat onsets durations names pmod;
save regs.mat R