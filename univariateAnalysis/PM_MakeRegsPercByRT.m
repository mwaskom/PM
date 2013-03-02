function PM_MakeRegsPercByConf(par)


cd (par.behavdir);

dFN = dir('*Perc*');
fileNames = {dFN.name}';

countdown = 12; %length of countdown period (in secs) at the beginning of each run 
         
         


[res psy idx] = Perceptual_fMRIBehAnalysis(par);

perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
resp_set = {'respFace' 'respHouse'};


i = 0;

for l=1:length(class_set)
    for p=1:length(perf_set)
        
        idx.thisOns = idx.(perf_set{p}) .* idx.(class_set{l}) .* (idx.rt > .25);
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_%s', class_set{l}, perf_set{p});
            
            
            stimOnsets{i}= idx.alltrials(find(idx.thisOns));
            stimNames{i} = fName;
            stimDurations{i} = 0;
            
            pmod(i).param = {idx.rt(find(idx.thisOns))};
            pmod(i).name = {['RT_' fName]};
            pmod(i).poly = {1};
        end
    end
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
save ons.mat onsets durations names pmod;
save regs.mat R