function [res psyphys idx reh] = Mnemonic_fMRIBehAnalysis_Retrieval(par)

cwd = pwd;
rtThresh = .25; %RT less than this thresh (in sec) are excluded from analysis


%% read in appropriate data
hands = {{'1' '2' '3' '4'} {'6' '7' '8' '9'} };
certRatTrans = {[0 8 7 6 5 1 2 3 4 ] [0 1 2 3 4 8 7 6 5] [0 4 3 2 1 5 6 7 8]};


if ismember(par.substr, {'pm_120810' 'pm_060111'})
    retFileIdx = [1 3 4];
    %for this subject, runs 1, 3, and 4 correspond to retrieval
elseif ismember(par.substr, { 'pm_042811' 'pm_050411_2' 'pm_052011' 'pm_060111'})
    retFileIdx = 2:4;
elseif ismember(par.substr, {'pm_040312'})
    retFileIdx = [2 3 4 6 ];
else
    retFileIdx = 1:3;
    %for all other subjects, it's runs 1, 2, and 3.
end


reh.wordList = {}; %presented words
reh.stimResp = {}; % for responses while the item is still displayed
reh.judgeResp = {}; % for responses after the item has left the screen
reh.judgeRT = {}; % normal RT
reh.stimRT = {}; % early RT (collected during .5 s of stimulus presentation
reh.cond = []; % 1 = face, 2 = house
reh.allTrials = [];
reh.sess = [];

cd (par.behavdir);

if strcmp(par.substr, 'pm_011611')
    %this subject was collected under unusal circumstances (a slightly
    %older version of the task, and 180 words, including 30 lures, were
    %tested, instead of the usual 150 non-lures).  This subject's results
    %were manually tabulated, and then shuffled back into a matlab script
    %that can be read by the present script.  
    
    RL.dat = load('retrievalDataManual');
    
    idx.lure = RL.dat.retData.lure==1;
    
    reh.wordList = RL.dat.retData.word;
    reh.cond = RL.dat.retData.cond;
    reh.resp = RL.dat.retData.resp;
    reh.RT = RL.dat.retData.RT;
    
    idx.alltrials = round(RL.dat.retData.onset)';
    
    reh.cond(find(idx.lure)) = 0;  %remove responses to lures
    reh.resp(find(idx.lure)) = {'0'}; %remove responses to lures

else
    
    for r = 1:length(par.scans_to_include)

        % load the relevant data file
        thisReh = dir(['*retrieve' '*(' num2str(retFileIdx(r)) ')*']);
        thisReh = thisReh(find(cellfun(@(x) ~strcmp(x(1),'.'), {thisReh(:).name})));
        
        RL(r).dat = load(thisReh.name);
        
        
        %fields of interest from the data file
        
        if (strcmp(par.substr, 'pm_032811') && (r==1))
            trialIdx = 1:30;
        elseif (strcmp(par.substr, 'pm_040312') && (r==2))
            trialIdx = 1:45;
        else
            trialIdx = find(RL(r).dat.theData.onset>0);
        end
        
        reh.wordList = vertcat(reh.wordList, RL(r).dat.list.words(trialIdx));
        reh.stimResp = horzcat(reh.stimResp, RL(r).dat.theData.stimresp(trialIdx));
        reh.judgeResp = horzcat(reh.judgeResp, RL(r).dat.theData.judgeResp(trialIdx));
        reh.judgeRT = horzcat(reh.judgeRT,  RL(r).dat.theData.judgeRT(trialIdx)); %becuase the stimulus is presented for .5 s before judgeResp collection
        reh.stimRT = horzcat(reh.stimRT, RL(r).dat.theData.stimRT(trialIdx));
        reh.cond = horzcat(reh.cond, RL(r).dat.theData.cond(trialIdx));
        reh.sess = horzcat(reh.sess, r*ones(size(RL(r).dat.theData.cond(trialIdx))));
        
        alltrials_h{r} = RL(r).dat.theData.onset(trialIdx)-12 + 2*sum(par.numvols(1:r-1));        
    end
 
    
    idx.alltrials = round(horzcat(alltrials_h{:}));
    
    %takes the first recorded response; also, changes '1!' into '1'
    judgeRespProcessed =  arrayfun(@(x) x{1}(1), reh.judgeResp, 'UniformOutput',false);
    
    %takes the first recorded RT
    stimRespProcessed =  arrayfun(@(x) x{1}(1), reh.stimResp, 'UniformOutput',false);
    
    %add .5 to recorded RTs, because the stims were presented for .5 s
    %before the judge responses were collected
    judgeRTProcessed_h = cell2mat(arrayfun(@(x) x{1}(1), reh.judgeRT, 'UniformOutput',false));
    judgeRTProcessed = judgeRTProcessed_h +.5;
    judgeRTProcessed(judgeRTProcessed==.5)=0;
    
    stimRTProcessed = cell2mat(arrayfun(@(x) x{1}(1), reh.stimRT, 'UniformOutput',false));
    
    idx.earlyResp = find(~strcmp(stimRespProcessed, 'n'));
    respVec = judgeRespProcessed;
    RTVec = judgeRTProcessed;
    
    %when subjects responded early (i.e. during stimulus presentation,
    %count their response and RT)
    respVec(idx.earlyResp) = stimRespProcessed(idx.earlyResp);
    RTVec(idx.earlyResp) = stimRTProcessed(idx.earlyResp);
    
    %vectors of cleaned responses and RT
    reh.resp = respVec;
    reh.RT = RTVec;
end




%% relevant indices

%when faces and houses were presented
idx.face = reh.cond==1;
idx.house = reh.cond==2;

% hand number is counterbalanced across subjects
HN = RL(1).dat.S.retHandNum;
CRTFunc = HN;

 if strcmp(par.substr, 'pm_040611')
     %for this subject during this session, I'm presuming that 
     %the left box was in the right hand, and vice versa.  
     %the data are completely consistent with this assumption...
     HN = 2;
     CRTFunc = 3;
 end


% takes a raw reh.resp, a cell array of responses (1-4 and 5-9) and turn
% it into reh.certRat, a mat of numbers 1-8, where 1 signifies highest confidence house
% and 8 highest confidence face.
respNum_h = reh.resp;
respNum_h(strcmp(reh.resp, 'n'))= {'0'};
respNum_h =  arrayfun(@(x) str2num(x{1}), respNum_h, 'UniformOutput',false);
respNum_h = cell2mat(respNum_h);
respNum_h(respNum_h>4) = respNum_h(respNum_h>4)-1; %change '6-9' to '5-8';

for i = [0:8];
    reh.certRat(respNum_h==i)=certRatTrans{CRTFunc}(i+1);
end

% when subjects responded 'face' or 'house'
idx.respFace = ismember(reh.resp, hands{HN});
idx.respHouse =  ismember(reh.resp, hands{3-HN});

h.idx.cor.Face = idx.respFace.*idx.face;
h.idx.cor.House = idx.respHouse.*idx.house;

% index of all correct and all incorrect trials
idx.cor = idx.respFace.*idx.face + idx.respHouse.*idx.house;
idx.inc = idx.respFace.*idx.house + idx.respHouse.*idx.face;

% index of trials for which an interpretable response was recorded
idx.exceedsRTThresh = (reh.RT>rtThresh); %is the RT above a minimum threshold?
idx.cleanResp = (idx.respFace + idx.respHouse) .* idx.exceedsRTThresh;

idx.certRat = reh.certRat;

idx.junk = ~idx.cleanResp;

%from low to high confidence
idx.C{1} = ismember(reh.certRat, [4 5]);
idx.C{2} = ismember(reh.certRat, [3 6]);
idx.C{3} = ismember(reh.certRat, [2 7]);
idx.C{4} = ismember(reh.certRat, [1 8]);

idx.cresp = 1*idx.C{1} + 2*idx.C{2} + 3*idx.C{3} + 4*idx.C{4};

% binary hi vs. low
idx.low = ismember(reh.certRat, 2:7);
idx.high = ismember(reh.certRat, [1 8]);

% unsigned confidence from 1-4 regardless of face/house status
idx.unsignedConf = abs(round(idx.certRat - 4.5));

idx.sess = reh.sess;
%% data from learning phase 
[idxL] = Mnemonic_fMRIBehAnalysis_Learning(par);

[~, encToRetTransform] = ismember(reh.wordList, idxL.item);
[~, retToEncTransform] = ismember(idxL.item, reh.wordList);

encToRetTransform(encToRetTransform==0) = [];
retToEncTransform(retToEncTransform==0) = [];

lureIdx = setdiff(1:length(idxL.item), encToRetTransform);
sum(idxL.respOld(lureIdx)) ./ sum(idxL.respOld(lureIdx) + idxL.respNew(lureIdx));

idx.enc.respOld = zeros(size(idx.cor));
idx.enc.respNew = zeros(size(idx.cor));

idx.enc.respOld(retToEncTransform) = idxL.respOld(encToRetTransform);
idx.enc.respNew(retToEncTransform) = idxL.respNew(encToRetTransform);

res.encHits = sum(idx.enc.respOld) ./ (sum(idx.enc.respOld + idx.enc.respNew));
res.encFAs = sum(idxL.respOld(lureIdx)) ./ sum(idxL.respOld(lureIdx) + idxL.respNew(lureIdx));

%-- Calculate dprime
NOld = sum(idx.enc.respOld + idx.enc.respNew);
NLure = sum(idxL.respOld(lureIdx) + idxL.respNew(lureIdx));

if res.encHits==1 % if 100% Hits
    pHit = 1-(1/(2*NOld)); %
else
    pHit = res.encHits;
end

if res.encFAs == 0 % if 0% FA
    pFA = 1/(2*NLure);
elseif res.encFAs == 1
    pFA = 1 - 1/(2*NLure);
else
    pFA = res.encFAs;
end

%-- Convert to Z scores
zHit = norminv(pHit) ;
zFA = norminv(pFA) ;

%-- Calculate d-prime
res.encDPrime = zHit - zFA ;

if length(idx.cor) ~= length(idx.enc.respOld)
   error('length of idx.enc variables is different from regular idx variables') 
end

%% results section

%accuracy results
res.pctCorFace = sum(idx.face .* idx.cor .*idx.cleanResp)/sum(idx.face .*idx.cleanResp);
res.pctCorHouse = sum(idx.house .* idx.cor .*idx.cleanResp)/sum(idx.house .*idx.cleanResp);
res.pctHouseResps = sum(idx.respHouse .* idx.cleanResp)/sum(idx.cleanResp);

%how many trials have a recorded response?
res.usableN = sum(idx.cleanResp); 

for CI = [1:8]
   thisCI = ['confBin' num2str(CI)];
   res.(thisCI) =  sum(reh.certRat == CI);
end

res.meanConfCorFace = mean(idx.cresp(logical(idx.face .* idx.cor .*idx.cleanResp)));
res.meanConfCorHouse = mean(idx.cresp(logical(idx.house .* idx.cor .*idx.cleanResp)));

%rt results
res.rtCor = median(reh.RT(find(idx.cor .*idx.cleanResp)));
res.rtInc = median(reh.RT(find(idx.inc .*idx.cleanResp)));
res.rtFace = median(reh.RT(find(idx.cleanResp .* idx.face)));
res.rtHouse = median(reh.RT(find(idx.cleanResp .* idx.house)));
res.rtFaceCor = median(reh.RT(find(idx.cor .*idx.cleanResp .* idx.face)));
res.rtHouseCor = median(reh.RT(find(idx.cor .*idx.cleanResp .* idx.house)));
res.rtFaceInc = median(reh.RT(find(idx.inc .*idx.cleanResp .* idx.face)));
res.rtHouseInc = median(reh.RT(find(idx.inc .*idx.cleanResp .* idx.house)));

idx.rtCor = reh.RT(find(idx.cor ));
idx.rt = reh.RT;

%% psychophysics plots - signed
for i = 1:8
    
    %iHand converts the index 'i' to the subject's response.  For instance, for some
    %subjects '1' is highest confidence face, for others, '1' is highest
    %confidence 'house'
    %Here, 1 is highest confidence house, to 8 which is highest confidence
    %Face
    %iHand = handIdxHouseToFace{HN}(i);
    
    
    %percent correct for responses ranging from high conf fact to high conf
    %house
    psyphys.AccuracyBySignedConf.mean(i) = sum(ismember(reh.certRat, i) .* idx.cor .* idx.cleanResp) / sum(ismember(reh.certRat, i) .* idx.cleanResp);
    psyphys.AccuracyBySignedConf.SE(i) = sqrt((psyphys.AccuracyBySignedConf.mean(i) * (1-(psyphys.AccuracyBySignedConf.mean(i)))) / sum(ismember(reh.certRat, i).* idx.cleanResp) );
    psyphys.AccuracyBySignedConf.ticks.x = 0:9;
    psyphys.AccuracyBySignedConf.ticks.y = [0 1];
    psyphys.AccuracyBySignedConf.xlabel = 'Conf (high House to high Face)';
    psyphys.AccuracyBySignedConf.ylabel = 'percent cor' ;
    psyphys.AccuracyBySignedConf.xVals = 1:8;
    psyphys.AccuracyBySignedConf.marker = 'o-';
    
    %choice for responses ranging from high conf fact to high conf
    %house
     psyphys.pFace.mean(i) = sum(ismember(reh.certRat, i) .* idx.respFace .* idx.cleanResp) / sum(ismember(reh.certRat, i) .* idx.cleanResp);
     psyphys.pFace.SE(i) = sqrt((psyphys.pFace.mean(i) * (1-(psyphys.pFace.mean(i)))) / sum(ismember(reh.certRat, i).* idx.cleanResp) );
     psyphys.pFace.ticks.x = 0:9;
     psyphys.pFace.ticks.y = [0 1];
     psyphys.pFace.xlabel = 'Conf (high House to high Face)';
     psyphys.pFace.ylabel = 'percent face choice' ;
     psyphys.pFace.xVals = 1:8;
     psyphys.pFace.marker = 'o-';
%     
    %RT for responses ranging from high conf fact to high conf
    %house
    psyphys.CorrectRTBySignedConf.mean(i) = median(reh.RT(find((ismember(reh.certRat, i) .* idx.cor .* idx.cleanResp))));
    psyphys.CorrectRTBySignedConf.SE(i) = ste(reh.RT(find((ismember(reh.certRat, i) .* idx.cor .* idx.cleanResp))));
    psyphys.CorrectRTBySignedConf.ticks.x = 0:9;
    psyphys.CorrectRTBySignedConf.ticks.y = [0 4];
    psyphys.CorrectRTBySignedConf.xlabel = 'Conf (highH to highF)';
    psyphys.CorrectRTBySignedConf.ylabel = 'RT (correct trials only)' ;
    psyphys.CorrectRTBySignedConf.xVals = 1:8;
    psyphys.CorrectRTBySignedConf.marker = 'o-';
    
    %distribution of responses ranging from high conf fact to high conf
    %house
    psyphys.RespDistBySignedConf.mean(i) = sum(ismember(reh.certRat, i));
    psyphys.RespDistBySignedConf.SE(i) = 0;
    psyphys.RespDistBySignedConf.ticks.x = 0:9;
    psyphys.RespDistBySignedConf.ticks.y = [0 60];
    psyphys.RespDistBySignedConf.xlabel = 'Conf (highH to highF)';
    psyphys.RespDistBySignedConf.ylabel = 'Number of responses' ;
    psyphys.RespDistBySignedConf.xVals = 1:8;
    psyphys.RespDistBySignedConf.marker = 'o-';
end




%% psychophysics plots - unsigned

for i = 1:4
    
    %which button presses correspond to this level of confidence (1 =
    %lowest confidence, 4 = highest)
    iHand = [(5-i) (4+i)];
    
    %percent correct for lowest to highest confidence
    psyphys.AccuracyByUnsignedConf.mean(i) = sum(ismember(reh.certRat, iHand) .* idx.cor) / sum(ismember(reh.certRat, iHand));
    psyphys.AccuracyByUnsignedConf.SE(i) = sqrt((psyphys.AccuracyByUnsignedConf.mean(i) * (1-(psyphys.AccuracyByUnsignedConf.mean(i)))) / sum(ismember(reh.certRat, iHand)));
    psyphys.AccuracyByUnsignedConf.ticks.x = 0:5;
    psyphys.AccuracyByUnsignedConf.ticks.y = [0 1];
    psyphys.AccuracyByUnsignedConf.xlabel = 'Conf (low to high)';
    psyphys.AccuracyByUnsignedConf.ylabel = 'percent cor' ;
    psyphys.AccuracyByUnsignedConf.xVals = 1:4;
    psyphys.AccuracyByUnsignedConf.marker = 'o-';
    
    
    %RT for lowest to highest confidence
    psyphys.CorRTByUnsignedConf.mean(i) = median(reh.RT(find((ismember(reh.certRat, iHand) .* idx.cor))));
    psyphys.CorRTByUnsignedConf.SE(i) = ste(reh.RT(find((ismember(reh.certRat, iHand) .* idx.cor))));
    psyphys.CorRTByUnsignedConf.ticks.x = 0:5;
    psyphys.CorRTByUnsignedConf.ticks.y = [0 4];
    psyphys.CorRTByUnsignedConf.xlabel = 'Conf (low to high)';
    psyphys.CorRTByUnsignedConf.ylabel = 'RT (correct trials only)' ;
    psyphys.CorRTByUnsignedConf.xVals = 1:4;
    psyphys.CorRTByUnsignedConf.marker = 'o-';
    
        %distribution of responses
    psyphys.RespDistByUnsignedConf.mean(i) = sum(ismember(reh.certRat, iHand) );
    psyphys.RespDistByUnsignedConf.SE(i) = 0;
    psyphys.RespDistByUnsignedConf.ticks.x = 0:5;
    psyphys.RespDistByUnsignedConf.ticks.y = [0 100];
    psyphys.RespDistByUnsignedConf.xlabel = 'Conf (low to high)';
    psyphys.RespDistByUnsignedConf.ylabel = 'Number of responses' ;
    psyphys.RespDistByUnsignedConf.xVals = 1:4;
    psyphys.RespDistByUnsignedConf.marker = 'o-';
end



cd (cwd)
%% stats

[~, reh.stats.pctCorByConf.p] = corr(idx.cor, idx.cresp);

 [~, ~, reh.stats.pctCorByConf.stats] = mnrfit(idx.cresp', (idx.cor+1)');
 [~, reh.stats.RTByConf.p] = corr(reh.RT', idx.cresp');
