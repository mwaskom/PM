function [res psyphys idx] = Mnemonic_fMRIBehAnalysis_Retrieval(par)


%%This script is not functional yet; it's basically a copy of
%%Mnemonic_fMRIBehAnalysis_Retrieval so far.  The reason is I just
%%remembered that we did not record responses during the rehearsal phase:
%%we must fix this.  


%% read in appropriate data
hands = {{'1' '2' '3' '4'} {'5' '6' '7' '8'}};
handIdx = {1:8, 8:-1:1};

rehFileIdx = 1:18;

reh.wordList = {}; %presented words
reh.stimResp = {}; % for responses while the item is still displayed
reh.judgeResp = {}; % for responses after the item has left the screen
reh.judgeRT = {}; % normal RT
reh.stimRT = {}; % early RT (collected during .5 s of stimulus presentation
reh.cond = []; % 1 = face, 2 = house

cd (par.behavdir);



for r = 1:18
    
    % load the relevant data file
    thisReh = dir(['*rehearse' '*(' num2str(rehFileIdx(r)) ')*']);
    RL(r).dat = load(thisReh.name);
    
    
    %fields of interest from the data file
    reh.wordList = vertcat(reh.wordList, RL(r).theData.item);
    reh.stimResp = horzcat(reh.stimResp, RL(r).dat.theData.stimresp);
    reh.judgeResp = horzcat(reh.judgeResp, RL(r).dat.theData.judgeResp);
    reh.judgeRT = horzcat(reh.judgeRT,  RL(r).dat.theData.judgeRT); %becuase the stimulus is presented for .5 s before judgeResp collection
    reh.stimRT = horzcat(reh.stimRT, RL(r).dat.theData.stimRT);
    reh.cond = horzcat(reh.cond, RL(r).dat.theData.cond);
    
    alltrials_h{r} = RL(r).dat.theData.onset-12 + 2*sum(par.numvols(1:(retFileIdx(r)-1)));
    
    
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


%% relevant indices

%when faces and houses were presented
idx.face = reh.cond==1;
idx.house = reh.cond==2;

% hand number is counterbalanced across subjects
HN = RL(1).dat.S.retHandNum;

% when subjects responded 'face' or 'house'
idx.respFace = ismember(reh.resp, hands{HN});
idx.respHouse =  ismember(reh.resp, hands{3-HN});

h.idx.cor.Face = idx.respFace.*idx.face;
h.idx.cor.House = idx.respHouse.*idx.house;

% index of all correct and all incorrect trials
idx.cor = idx.respFace.*idx.face + idx.respHouse.*idx.house;
idx.inc = idx.respFace.*idx.house + idx.respHouse.*idx.face;

% index of trials for which an interpretable response was recorded
idx.cleanResp = (idx.respFace + idx.respHouse) ;

idx.junk = ~idx.cleanResp;

%from low to high confidence
idx.C{1} = ismember(reh.resp, {'4' '5'});
idx.C{2} = ismember(reh.resp, {'3' '6'});
idx.C{3} = ismember(reh.resp, {'2' '7'});
idx.C{4} = ismember(reh.resp, {'1' '8'});

% binary hi vs. low
idx.low = ismember(reh.resp, {'3' '4' '5' '6'});
idx.high = ismember(reh.resp, {'1' '2' '7' '8'});

%% results section

%accuracy results
res.pctCorFace = sum(idx.face .* idx.cor .*idx.cleanResp)/sum(idx.face .*idx.cleanResp);
res.pctCorHouse = sum(idx.house .* idx.cor .*idx.cleanResp)/sum(idx.house .*idx.cleanResp);
res.pctHouseResps = sum(idx.respHouse .* idx.cleanResp)/sum(idx.cleanResp);

%how many trials have a recorded response?
res.usableN = sum(idx.cleanResp); 

%rt results
res.rtCor = median(reh.RT(find(idx.cor .*idx.cleanResp)));
res.rtInc = median(reh.RT(find(idx.inc .*idx.cleanResp)));
res.rtFace = median(reh.RT(find(idx.cleanResp .* idx.face)));
res.rtHouse = median(reh.RT(find(idx.cleanResp .* idx.house)));
res.rtFaceCor = median(reh.RT(find(idx.cor .*idx.cleanResp .* idx.face)));
res.rtHouseCor = median(reh.RT(find(idx.cor .*idx.cleanResp .* idx.house)));
res.rtFaceInc = median(reh.RT(find(idx.inc .*idx.cleanResp .* idx.face)));
res.rtHouseInc = median(reh.RT(find(idx.inc .*idx.cleanResp .* idx.house)));



%% psychophysics plots - signed
for i = 1:8
    
    %iHand converts the index 'i' to the subject's response.  For instance, for some
    %subjects '1' is highest confidence face, for others, '1' is highest
    %confidence 'house'
    iHand = handIdx{HN}(i);
    
    
    %percent correct for responses ranging from high conf fact to high conf
    %house
    psyphys.pCor.mean(i) = sum(ismember(reh.resp, num2str(iHand)) .* idx.cor) / sum(ismember(reh.resp, num2str(iHand)));
    psyphys.pCor.SE(i) = sqrt((psyphys.pCor.mean(i) * (1-(psyphys.pCor.mean(i)))) / sum(ismember(reh.resp, num2str(iHand))));
    psyphys.pCor.ticks.x = 0:9;
    psyphys.pCor.ticks.y = [0 1];
    psyphys.pCor.xlabel = 'Conf (high Face to high House)';
    psyphys.pCor.ylabel = 'percent cor' ;
    
    %RT for responses ranging from high conf fact to high conf
    %house
    psyphys.RTCor.mean(i) = median(reh.RT(find((ismember(reh.resp, num2str(iHand)) .* idx.cor))));
    psyphys.RTCor.SE(i) = ste(reh.RT(find((ismember(reh.resp, num2str(iHand)) .* idx.cor))));
    psyphys.RTCor.ticks.x = 0:9;
    psyphys.RTCor.ticks.y = [0 4];
    psyphys.RTCor.xlabel = 'Conf (high Face to high House)';
    psyphys.RTCor.ylabel = 'RT (correct trials only)' ;
end


% % percent correct across confidence levels
% figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
% errorbar(1:8, psyphys.pCor.mean, psyphys.pCor.SE, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
% set(gca, 'XLim', [0 9], 'XTick', 0:1:9, 'XTickLabel', makeTickLabel(0:1:9,1), ...
%     'YLim', [0 1], 'TickDir', 'out');
% grid on;
% 
% xlabel('Confidence (highest conf face to highest conf house)');
% ylabel('percent cor');
% 
% 
% % RT across confidence levels (correct trials)
% figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
% errorbar(1:8, psyphys.RTCor.mean, psyphys.RTCor.SE, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
% set(gca, 'XLim', [0 9], 'XTick', 0:1:9, 'XTickLabel', makeTickLabel(0:1:9,1), ...
%     'YLim', [0 4], 'TickDir', 'out');
% grid on;
% 
% xlabel('Confidence (highest conf face to highest conf house)');
% ylabel('RT (correct trials only)');

%% psychophysics plots - unsigned

for i = 1:4
    
    %which button presses correspond to this level of confidence (1 =
    %lowest confidence, 4 = highest)
    iHand = {num2str(5-i) num2str(4+i)};
    
    %percent correct for lowest to highest confidence
    psyphys.unsignedpCor.mean(i) = sum(ismember(reh.resp, iHand) .* idx.cor) / sum(ismember(reh.resp, iHand));
    psyphys.unsignedpCor.SE(i) = sqrt((psyphys.pCor.mean(i) * (1-(psyphys.pCor.mean(i)))) / sum(ismember(reh.resp, iHand)));
    psyphys.unsignedpCor.ticks.x = 0:5;
    psyphys.unsignedpCor.ticks.y = [0 1];
    psyphys.unsignedpCor.xlabel = 'Confidence (low to high)';
    psyphys.unsignedpCor.ylabel = 'percent cor' ;
    
    
    %RT for lowest to highest confidence
    psyphys.unsignedRTCor.mean(i) = median(reh.RT(find((ismember(reh.resp, iHand) .* idx.cor))));
    psyphys.unsignedRTCor.SE(i) = ste(reh.RT(find((ismember(reh.resp, iHand) .* idx.cor))));
    psyphys.unsignedRTCor.ticks.x = 0:5;
    psyphys.unsignedRTCor.ticks.y = [0 4];
    psyphys.unsignedRTCor.xlabel = 'Confidence (low to high)';
    psyphys.unsignedRTCor.ylabel = 'RT (correct trials only)' ;
end


% % plot percent correct across confidence levels
% figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
% errorbar(1:4, psyphys.unsignedpCor.mean, psyphys.unsignedpCor.SE, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
% set(gca, 'XLim', [0 5], 'XTick', 0:1:5, 'XTickLabel', makeTickLabel(0:1:5,1), ...
%     'YLim', [0 1], 'TickDir', 'out');
% grid on;
% xlabel('Confidence (low to high)');
% ylabel('percent cor');
% 
% 
% % plot RT across confidence levels (correct trials)
% figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
% errorbar(1:4, psyphys.unsignedRTCor.mean, psyphys.unsignedRTCor.SE, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
% set(gca, 'XLim', [0 5], 'XTick', 0:1:5, 'XTickLabel', makeTickLabel(0:1:5,1), ...
%     'YLim', [0 4], 'TickDir', 'out');
% grid on;
% 
% xlabel('Confidence (low to high)');
% ylabel('RT (correct trials only)');
