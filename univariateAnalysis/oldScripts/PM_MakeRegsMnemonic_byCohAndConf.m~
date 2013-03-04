function PM_MakeRegsMnemonic_byEvidence(par)

meanAcrossIters = 0;

[~, ~, idx] = Mnemonic_fMRIBehAnalysis_Retrieval(par);


perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
conf_set = {'low' 'high'};


i = 0;

for l=1:length(class_set)
    
    for p=1:length(perf_set)
        
        %for c = 1:length(conf_set)
            idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l});
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_%s', class_set{l},  perf_set{p} );
                
                stimOnsets{i}= idx.alltrials(find(idx.thisOns .* idx.cleanResp));
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
                pmod(i).param{1} = idx.cresp(idx.thisOns .* idx.cleanResp);
                pmod(i).param{1} = idx.cresp(idx.thisOns .* idx.cleanResp);
            end
            
        %end
    end
end
stimOnsets{i+1} = idx.alltrials(find(idx.junk));
stimNames{i+1} = 'junk';
stimDurations{i+1} = 0;





pmod(1).param = {actDiffsPersonCorGuess};
pmod(1).name = {'GreaterEvPersonCorGuess'};
pmod(1).poly = {1};

pmod(2).param = {actDiffsHouseCorGuess};
pmod(2).name = {'GreaterEvSceneCorGuess'};
pmod(2).poly = {1};

fn = fieldnames(ons.AD);

for f = 1:length(fn);
    onsets{f} = ons.AD.(fn{f});
    durations{f} = 0;
    names{f} = fn{f};
end


if ~exist(par.analysisdir)
    mkdir (par.analysisdir);
end



sessReg = zeros(sum(par.numvols),par.numscans-1);
for i = 1:(par.numscans - 1)
    sessReg((i-1)*par.numvols(i)+1:i*par.numvols(i),i) = ones(par.numvols(i),1);
end

R = sessReg;


cd (par.analysisdir);
save ons.mat onsets durations names pmod;
save regs.mat R