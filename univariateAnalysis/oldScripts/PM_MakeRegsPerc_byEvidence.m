function PM_MakeRegsEvidence(par)


cd (par.behavdir);

dFN = dir('*Perc*.mat');
fileNames = {dFN.name}';

countdown = 12; %length of countdown period (in secs) at the beginning of each run 
         
[S] = PM_PatternParams(par.subNo,'perc');
[res psy idx] = Perceptual_fMRIBehAnalysis(par);

%%

% load classmat
classMat = load(S.classMatFile);
subjArrayIdx = classMat.qqq.subjArray==par.subNo ;

% load relevant performance vectors
for k = 1:length(classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations)
corTrials_h{k} = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.corrects;
allActs_h{k} = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).acts(1,:);
allDesireds_h{k} = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.desireds;
end

allActs = horzcat(allActs_h{:});
corTrials = horzcat(corTrials_h{:});
allDesidreds = horzcat(allDesireds_h{:});
allClass = idx.face + 2*idx.house;
allRespClass = idx.resp_face + 2*idx.resp_house;

allRT = idx.rt;
allCoh = idx.coh_;
allConf = idx.conf_unsigned; 


allActsLogit = log(allActs ./ (1- allActs)); %use logit instead of probability

% turn this from face evidence into general evidence in the correct
% direction, by flipping evidence for scenes.

%allActs(allClass==2) = -allActs(allClass==2);
allActsLogit(allClass==2) = -allActsLogit(allClass==2);

%%

perf_set = {'corWithZeros' 'inc'};
class_set = {'face' 'house'};
resp_set = {'respFace' 'respHouse'};


i = 0;

for p=1:length(perf_set)
    for l=1:length(class_set)
        idx.thisOns =  idx.(perf_set{p}) .* idx.(class_set{l});
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            fName = sprintf('%s_%s', class_set{l} , perf_set{p});
            
            stimOnsets{i}= idx.alltrials(find(idx.thisOns));
            stimNames{i} = fName;
            stimDurations{i} = 0;
        end
    end
end


evFace = allActsLogit((allClass==1 .* idx.corWithZeros)==1);
evHouse = allActsLogit((allClass==2 .* idx.corWithZeros)==1);

pmod(1).param = {evFace};
pmod(1).name = {'evidenceFace'};
pmod(1).poly = {1};

pmod(2).param = {evHouse};
pmod(2).name = {'evidenceHouse'};
pmod(2).poly = {1};

if sum(idx.junk>0)
    i= i+1;
    stimOnsets{i} = idx.alltrials(find(idx.junk));
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




R = horzcat(sessReg);


cd (par.analysisdir);
save ons.mat onsets durations names pmod;
save regs.mat R