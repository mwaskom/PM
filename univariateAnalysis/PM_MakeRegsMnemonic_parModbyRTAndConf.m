function PM_MakeRegsMnemonic_parModbyRTAndConf(par)

meanAcrossIters = 0;

[~, ~, idx] = Mnemonic_fMRIBehAnalysis_Retrieval(par);


perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};


i = 0;

for l=1:length(class_set)
    
    for p=1:length(perf_set)
        
        idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l}) .* ~isnan(idx.rt);
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_%s', class_set{l},  perf_set{p} );
            
            stimOnsets{i}= idx.alltrials(find(idx.thisOns .* idx.cleanResp));
            stimNames{i} = fName;
            stimDurations{i} = 0;
            
            if length(stimOnsets{i})>1
                pmod(i).param{1} = idx.cresp(find(idx.thisOns .* idx.cleanResp));
                pmod(i).param{2} = idx.rt(find(idx.thisOns .* idx.cleanResp));
                
                pmod(i).name{1} = 'conf';
                pmod(i).name{2} = 'RT';
                
                pmod(i).poly{1} = 1;
                pmod(i).poly{2} = 1;
            end
        end
        
    end
end

if sum(idx.junk + isnan(idx.rt))>0
    stimOnsets{i+1} = idx.alltrials(find(idx.junk + isnan(idx.rt)));
    stimNames{i+1} = 'junk';
    stimDurations{i+1} = 0;
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


R = sessReg;

cd (par.analysisdir);
save ons.mat onsets durations names pmod;
save regs.mat R