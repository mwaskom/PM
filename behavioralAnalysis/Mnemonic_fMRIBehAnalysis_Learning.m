function [lrn] = Mnemonic_fMRIBehAnalysis_Learning(par)

lrn.item = {}; %presented words
lrn.cond = []; % for responses while the item is still displayed
lrn.judgeResp = {}; % normal RT
lrn.judgeRT = {}; % early RT (collected during .5 s of stimulus presentation

cd (par.behavdir);

if strcmp(par.substr, 'pm_052011')
    rIdx = 10:17;
else
    rIdx = 10:18;
end

for r = rIdx   
    % load the relevant data file
    thisLrn = dir(['*rehearse' '*(' num2str(r) ')*']);
    thisLrn = thisLrn(find(cellfun(@(x) ~strcmp(x(1),'.'), {thisLrn(:).name})));
    
    RL(r).dat = load(thisLrn.name);
    
    lrn.item = vertcat(lrn.item, RL(r).dat.theData.item);
    lrn.cond = horzcat(lrn.cond, RL(r).dat.theData.cond);
    lrn.judgeResp = horzcat(lrn.judgeResp, RL(r).dat.theData.judgeResp);
    lrn.judgeRT = horzcat(lrn.judgeRT,  RL(r).dat.theData.judgeRT); %becuase the stimulus is presented for .5 s before judgeResp collection 
end

%takes the first recorded response; also, changes '1!' into '1'
judgeRespProcessed =  arrayfun(@(x) x{1}(1), lrn.judgeResp, 'UniformOutput',false);

%takes the first recorded RT
judgeRTProcessed = cell2mat(arrayfun(@(x) x{1}(1), lrn.judgeRT, 'UniformOutput',false));

respVec = judgeRespProcessed;
RTVec = judgeRTProcessed;

lrn.resp = respVec;
lrn.RT = RTVec;


%% relevant indices

hands = {'q' 'p'};
% hand number is counterbalanced across subjects
HN = RL(10).dat.S.encHandNum;

% when subjects responded 'face' or 'house'
lrn.respOld = ismember(lrn.resp, hands{HN});
lrn.respNew =  ismember(lrn.resp, hands{3-HN});

