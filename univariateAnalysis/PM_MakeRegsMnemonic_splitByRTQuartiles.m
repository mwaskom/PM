function PM_MakeRegsMnemonic_splitByRTQuartiles(par)

[~, ~, idx] = Mnemonic_fMRIBehAnalysis_Retrieval(par);


%perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
%conf_set = {'low' 'high'};
RT_set = 1:4;

qBins = quantile(idx.rt, 3);
[~, idx.RTBins] = histc(idx.rt, [0 qBins Inf]);


i = 0;

for l=1:length(class_set)
    
    %for p=1:length(perf_set)
        
        for c = 1:length(RT_set)
            %idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l}) .* idx.RTBins==c;
            idx.thisOns =  idx.(class_set{l}) .* idx.RTBins==c;
            %idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l}) .* idx.(conf_set{c});
            %idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l});
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                fName = sprintf('%s_RT%g', class_set{l},  c);
                %fName = sprintf('%s_%s_RT%g', class_set{l},  perf_set{p},  c);
                %fName = sprintf('%s_%s_%sconf', class_set{l},  perf_set{p},  conf_set{c});
                %fName = sprintf('%s_%s', class_set{l},  perf_set{p});
                
                stimOnsets{i}= idx.alltrials(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = 0;
                
%                 pmod(i).param = {idx.rt(find(idx.thisOns))};
%                 pmod(i).name = {['RT_' fName]};
%                 pmod(i).poly = {1};
            end
        end
    %end
end

if sum(idx.junk>0)
    stimOnsets{i+1} = idx.alltrials(find(idx.junk));
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
save ons.mat onsets durations names; %pmod;
save regs.mat R