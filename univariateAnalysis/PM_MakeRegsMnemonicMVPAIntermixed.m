function PM_MakeRegsMnemonicMVPAIntermixed(par)

[~, ~, idx] = Mnemonic_fMRIBehAnalysis_RetrievalIntermixed(par);


perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
conf_set = {'low' 'high'};


i = 0;

for l=1:length(class_set)
    
    for p=1:length(perf_set)
        
        for c = 1:length(conf_set)
            idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l}) .* idx.(conf_set{c});
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s_%s', class_set{l},  perf_set{p},  conf_set{c});
                
                stimOnsets{i}= idx.alltrials(find(idx.thisOns))';
                stimNames{i} = fName;
                stimDurations{i} = 0;
            end
            
        end
        
        idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l});
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_%s', class_set{l},  perf_set{p});
            
            stimOnsets{i}= idx.alltrials(find(idx.thisOns))';
            stimNames{i} = fName;
            stimDurations{i} = 0;
        end
    end
end

stimOnsets{i+1} = idx.alltrials(find(idx.junk))';
stimNames{i+1} = 'junk';
stimDurations{i+1} = 0;

onsets = stimOnsets;
names = stimNames;
durations = stimDurations;

if ~exist(par.analysisdir)
    mkdir (par.analysisdir);
end



sessReg = zeros(sum(par.numvols),par.numscans-1);
for i = 1:(par.numscans - 1)
    sessReg((i-1)*par.numvols(i)+1:i*par.numvols(i),i) = ones(par.numvols(i),1);
end


for i = 1:par.numscans
    cd(fullfile(par.subdir, 'functional', ['scan' (prepend(num2str(par.scans_to_include(i)),2))]   ));
    motTxt = dir('rp*');
    motRegs_h{i} = textread(motTxt.name);
end

motRegs = vertcat(motRegs_h{:});

R = horzcat(sessReg, motRegs);




cd (par.analysisdir);
save mvpa_ons.mat onsets names;
save regs.mat R