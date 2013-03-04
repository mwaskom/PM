function PM_MakeRegsPerc_parModByEvConfAndRT(par)

qqq_h = load(par.classmat);
qqq = qqq_h.res;
[~,~,dat,datB] = PM_classificationPostProcessingPerceptualImproved(qqq);

[~,subNo] = ismember(par.subNo, qqq.subjArray);

idx = dat{subNo};
idxB = datB{subNo};

cd (par.behavdir);

dFN = dir('*Perc*');
fileNames = {dFN.name}';

countdown = 12; %length of countdown period (in secs) at the beginning of each run 
         

%[res psy idx] = Perceptual_fMRIBehAnalysis(par);

perf_set = {'corWithZeros' 'inc'};
class_set = {'face' 'house'};
%resp_set = {'resp_face' 'resp_house'};


i = 0;

for l=1:length(class_set)
    for p=1:length(perf_set)
        
        idx.thisOns = idxB.(perf_set{p}) .* idxB.(class_set{l}) .* ~isnan(idxB.rt) .* idx.onsInClassifier;
        
        if ~isempty(idxB.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_%s', class_set{l}, perf_set{p});
            
            stimOnsets{i}= idxB.alltrials(find(idx.thisOns));
            stimNames{i} = fName;
            stimDurations{i} = 0;
            
            thisConf = idxB.conf_unsigned(find(idx.thisOns));
            thisRT = idxB.rt(find(idx.thisOns));
            
            idx.thisClassOns = idx.(perf_set{p}) .* idx.(class_set{l}) .* ~isnan(idx.rt);
            
            %pmod(i).param{3} = thisConf;
            %pmod(i).param{2} = thisRT;
            pmod(i).param{1} = idx.unsignedActs(find(idx.thisClassOns)); 
            
            %pmod(i).name{3} = 'Conf';
            %pmod(i).name{2} = 'RT';
            pmod(i).name{1} = 'Ev';
            
            %pmod(i).poly{3} = 1;
            %pmod(i).poly{2} = 1;
            pmod(i).poly{1} = 1;
        end
    end
end

idx.extraOns = (idxB.junk | ~idx.onsInClassifier | isnan(idxB.rt));
if sum(idx.extraOns)>0
    stimOnsets{i+1} = idxB.alltrials(idx.extraOns);
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


R = horzcat(sessReg);

cd (par.analysisdir);

save ons.mat onsets durations names pmod;
save regs.mat R