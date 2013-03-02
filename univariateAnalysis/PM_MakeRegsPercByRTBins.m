function PM_MakeRegsPercByRTBins(par)


cd (par.behavdir);

dFN = dir('*Perc*');
fileNames = {dFN.name}';

countdown = 12; %length of countdown period (in secs) at the beginning of each run 
         
         


[res psy idx] = Perceptual_fMRIBehAnalysis(par);

perf_set = {'corWithZeros' 'inc'};
class_set = {'face' 'house'};
resp_set = {'resp_face' 'resp_house'};

idxLegitRTs = (idx.rt>.3)& idx.corWithZeros;

RT_set = 1:4;


i = 0;

for l=1:length(class_set)
    %for p=1:length(perf_set)
    qBins = quantile(idx.rt(idxLegitRTs & idx.(resp_set{l})), 3);
    [~, idx.RTBins] = histc(idx.rt, [0 qBins Inf]);
    
    for c=RT_set
        idx.thisOns = (idx.RTBins==c) .* idx.(resp_set{l}) .* idxLegitRTs;
        %idx.thisOns = (idx.RTBins==c) .* idx.(resp_set{l});
        %idx.thisOns = (idx.RTBins==c) .* idx.(perf_set{p}) .* idx.(resp_set{l});
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            %fName = sprintf('%s_RT%s_%s', class_set{l}, num2str(c), perf_set{p});
            fName = sprintf('%s_cor_RT%s', class_set{l}, num2str(c));
            
            stimOnsets{i}= idx.alltrials(find(idx.thisOns));
            stimNames{i} = fName;
            stimDurations{i} = 0;
        end
    end
    %end
end

for l=1:length(class_set)
    idx.thisOns = idx.inc .* idx.(resp_set{l});
    if ~isempty(idx.alltrials(find(idx.thisOns)))
           i = i+1;
            fName = sprintf('%s_inc', class_set{l});
            
            stimOnsets{i}= idx.alltrials(find(idx.thisOns));
            stimNames{i} = fName;
            stimDurations{i} = 0;
    end
end

idxJunk = setdiff(idx.alltrials, vertcat(stimOnsets{:}));
if ~isempty(idxJunk)
    i = i+1;
    stimOnsets{i}= idxJunk;
    stimNames{i} = 'junk';
    stimDurations{i} = 0;
end

onsets = stimOnsets;
names = stimNames;
durations = stimDurations;

if ~exist(par.analysisdir)
    mkdir (par.analysisdir);
end



sessReg = zeros(sum(par.numvols),length(par.numvols) - 1);
for i = 1:(length(par.numvols) - 1)
    sessReg((i-1)*par.numvols(i)+1:i*par.numvols(i),i) = ones(par.numvols(i),1);
end

% for i = 1:par.numscans
%     cd(fullfile(par.subdir, 'functional', ['scan' (prepend(num2str(par.scans_to_include(i)),2))]   ));
%     motTxt = dir('rp*');
%     motRegs_h{i} = textread(motTxt.name);
% end
% 
% motRegs = vertcat(motRegs_h{:});

R = horzcat(sessReg);




cd (par.analysisdir);

save ons.mat onsets durations names;
save regs.mat R