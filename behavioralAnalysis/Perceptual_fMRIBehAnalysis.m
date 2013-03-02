function [res psyphys idx gRes] = Perceptual_fMRIBehAnalysis(par)

gRes=[];
rtThresh = .25; %RT less than this thresh (in sec) are excluded from analysis
excludeHighRTs = 0;
%% read in appropriate data
% hands = {[1 2 3 4] [5 6 7 8]};
% handIdx = {1:8, 8:-1:1};


cd (par.behavdir)
dFN = dir('*Perc*.mat');
fileNames = {dFN.name};

dotFiles = cell2mat(cellfun(@(x) x(1)=='.', fileNames, 'UniformOutput', false));
fileNames(find(dotFiles)) = []; %remove hidden files that are prepended with dots.

localizer = 0;
countdown = 12;

trial_data = combineDataFile(fileNames, par.behavdir);
trial_num = length(trial_data);


%extract stim presentation and behavioral variables
if isfield(trial_data,'response')
    resp = cat(1, trial_data.response);
    
    if strcmp(par.substr, 'pm_031711')
        resp(isnan(resp))=9; %'9' button presses were not recorded for this subject.  Assume that Nans reflect '9' button presses.
        resp(resp>4) = resp(resp>4)-1; %change '6 7 8 9' to '5 6 7 8';
    end
    
    if ismember(par.substr, {'pm_050811_2'  'pm_052611'})
        % for these subjects, the hand that indicates face/house status is
        % flipped.
        resp_ = ceil(resp/4);
        cresp_ = 4- (mod(resp-1, 4));
    elseif ismember(par.substr, { 'pm_031711'})
        % for this subject confidence ratings are inverted.
        resp_ = 3-ceil(resp/4);
        cresp_ = mod(resp-1, 4) + 1;
    elseif ismember(par.substr, { 'pm_042811_2'}) 
        % for this subject, the mapping is unique...
        resp_ = 3-ceil(resp/4);
        cresp_ = 5-ceil(abs(resp-4.5));
    else
        resp_ = 3-ceil(resp/4);
        cresp_ = 4- (mod(resp-1, 4) );
    end
    
end


%Transform resp to coordinates Alan can understand.  
%From 1 (most house like) to 4(least house like) to 5( least face like) to 8
%(most face like).

%certRatTrans = [8 7 6 5 1 2 3 4];
%certRatTrans = [4 3 2 1 5 6 7 8];
%[4 3 2 1 8 7 6 5]

certRatTrans = [1 2 3 4 8 7 6 5];
%certRatTrans = [4 3 2 1 8 7 6 5];

if strcmp(par.substr, 'pm_031711')
    certRatTrans = [4 3 2 1 5 6 7 8];
end

if ismember(par.substr, {'pm_050811_2' 'pm_052611'})
    certRatTrans = [8 7 6 5 1 2 3 4];
end
 
certRat = zeros(size(resp))';
for i = 1:8
    certRat(resp==i)=certRatTrans(i);
end

     


if all(isnan(resp_))
        %accept all the trials if we did not collect response
    valid_trials = 1 : trial_num;
elseif ~all(isnan(resp_)) && all(isnan(cresp_))
        %if we collected direction choices but not certainty responses
    valid_trials = find(~isnan(resp_));
else
        %if we collected both direction choices and certainty responses 
    valid_trials = find(~isnan(resp_) & ~isnan(cresp_));
end

resp = resp(valid_trials);
resp_ = resp_(valid_trials);
cresp_ = cresp_(valid_trials);
trial_data = trial_data(valid_trials);
time0 = cat(1, trial_data.time0);
start_t = cat(1, trial_data.start_t);
event_name = cat(1, trial_data.event_name);
if size(trial_data(1).event_time,1)==1 
    event_time = cat(1, trial_data.event_time);
else
    event_time = cat(2, trial_data.event_time)';
end
event_time = event_time + repmat( start_t-time0-countdown ,[1 size(event_time,2)]);
stim_on = event_time(strmatch('stim_on',event_name));
stim_off = event_time(strmatch('stim_off',event_name));
dur_ = stim_off - stim_on;
stim_ = cat(1, trial_data.stim_group);
coh_ = cat(1, trial_data.stim_coh_seq);
coh_signed = coh_;
coh_signed(stim_==2) = -coh_signed(stim_==2);
scan_ = cat(1, trial_data.ownership);       %which scan each trial belongs to
coh_set = unique(coh_);
coh_set_signed = unique(coh_signed);

rt = cat(1,trial_data.rt);

certRat = certRat(valid_trials)';
%********************** XXXXXXXXXXXXXXX    
%%
idx.face = (stim_ == 1) .* (coh_~=0);
idx.house = (stim_ == 2) .* (coh_~=0);

idx.resp_face = resp_ == 1;
idx.resp_house = resp_ == 2;


idx.cor = idx.face .* idx.resp_face + idx.house .* idx.resp_house;
idx.inc = idx.face .* idx.resp_house + idx.house .* idx.resp_face;

idx.corWithZeros = idx.cor + (coh_==0);

idx.coh_signed = coh_signed;

idx.conf_signed = certRat;
idx.conf_unsigned = ceil(abs(idx.conf_signed-4.5));
idx.conf_binary = 1+(abs(idx.conf_signed-4.5)>3);

idx.rt = rt;

idx.exceedsRTThresh = rt > rtThresh;

if excludeHighRTs
    idx.cleanResp = (idx.resp_face + idx.resp_house) .* idx.exceedsRTThresh;
else
    idx.cleanResp = (idx.resp_face + idx.resp_house);
end

idx.coh_set = coh_set;
idx.coh_set_pct = coh_set*100;
idx.coh_ = coh_;
idx.scan_ = scan_;

idx.alltrials =  (scan_-1).*par.numvols(scan_)' * par.TR + stim_on; 
idx.junk = ~(idx.cor + idx.inc);
%% results section

%accuracy results
res.pctCorFace = sum(idx.face .* idx.cor .*idx.cleanResp)/sum(idx.face .*idx.cleanResp);
res.pctCorHouse = sum(idx.house .* idx.cor .*idx.cleanResp)/sum(idx.house .*idx.cleanResp);
res.pctHouseResps = sum(idx.resp_house .* idx.cleanResp)/sum(idx.cleanResp);

res.usableN = sum(idx.cleanResp); 

for CI = 1:8 
   thisCI = ['confBin' num2str(CI)];
   res.(thisCI) =  sum(ismember(certRat, CI));
   idx.conf(ismember(certRat, CI)) = CI;
end


%rt results
res.rtCor = nanmedian(rt(find(idx.cor .*idx.cleanResp)));
res.rtInc = nanmedian(rt(find(idx.inc .*idx.cleanResp)));
res.rtFace = nanmedian(rt(find(idx.cleanResp .* idx.face)));
res.rtHouse = nanmedian(rt(find(idx.cleanResp .* idx.house)));
res.rtFaceCor = nanmedian(rt(find(idx.cor .*idx.cleanResp .* idx.face)));
res.rtHouseCor = nanmedian(rt(find(idx.cor .*idx.cleanResp .* idx.house)));
res.rtFaceInc = nanmedian(rt(find(idx.inc .*idx.cleanResp .* idx.face)));
res.rtHouseInc = nanmedian(rt(find(idx.inc .*idx.cleanResp .* idx.house)));




% %% confidence psychophysics plots - signed
for i = 1:8
    
    %iHand converts the index 'i' to the subject's response.  For instance, for some
    %subjects '1' is highest confidence face, for others, '1' is highest
    %confidence 'house'
    %iHand = handIdx{HN}(i);
    
    
    %percent correct for responses ranging from high conf house to high conf
    %face
    psyphys.AccuracyBySignedConf.mean(i) = sum(ismember(certRat, i) .* idx.cor .* idx.cleanResp) / sum(ismember(certRat, i) .* idx.cleanResp);
    psyphys.AccuracyBySignedConf.SE(i) = sqrt((psyphys.AccuracyBySignedConf.mean(i) * (1-(psyphys.AccuracyBySignedConf.mean(i)))) / sum(ismember(certRat, i).* idx.cleanResp));
    psyphys.AccuracyBySignedConf.ticks.x = 0:9;
    psyphys.AccuracyBySignedConf.ticks.y = [0 1];
    psyphys.AccuracyBySignedConf.xlabel = 'Conf (high House to high Face)';
    psyphys.AccuracyBySignedConf.ylabel = 'Percent Correct' ;
    psyphys.AccuracyBySignedConf.xVals = 1:8;
    psyphys.AccuracyBySignedConf.marker = 'o-';
    
    %percent correct for responses ranging from high conf house to high conf
    %face
    psyphys.pFaceRespBySignedConf.mean(i) = sum(ismember(certRat, i) .* idx.resp_face .* idx.cleanResp) / sum(ismember(certRat, i) .* idx.cleanResp);
    psyphys.pFaceRespBySignedConf.SE(i) = sqrt((psyphys.AccuracyBySignedConf.mean(i) * (1-(psyphys.AccuracyBySignedConf.mean(i)))) / sum(ismember(certRat, i).* idx.cleanResp));
    psyphys.pFaceRespBySignedConf.ticks.x = 0:9;
    psyphys.pFaceRespBySignedConf.ticks.y = [0 1];
    psyphys.pFaceRespBySignedConf.xlabel = 'Conf (high House to high Face)';
    psyphys.pFaceRespBySignedConf.ylabel = 'Percent FaceResp' ;
    psyphys.pFaceRespBySignedConf.xVals = 1:8;
    psyphys.pFaceRespBySignedConf.marker = 'o-';

    %RT for responses ranging from high conf face to high conf
    %house
    psyphys.RTBySignedConf.mean(i) = nanmedian(rt(find((ismember(certRat, i) .* idx.cor .* idx.cleanResp))));
    psyphys.RTBySignedConf.SE(i) = ste(rt(find((ismember(certRat, i) .* idx.cor .* idx.cleanResp))));
    psyphys.RTBySignedConf.ticks.x = 0:9;
    psyphys.RTBySignedConf.ticks.y = [0 4];
    psyphys.RTBySignedConf.xlabel = 'Conf (high House to high Face)';
    psyphys.RTBySignedConf.ylabel = 'RT (correct trials only)' ;
    psyphys.RTBySignedConf.xVals = 1:8;
    psyphys.RTBySignedConf.marker = 'o-';
    
        %distribution of responses ranging from high conf fact to high conf
    %house
    psyphys.RespDistBySignedConf.mean(i) = sum(ismember(certRat, i));
    psyphys.RespDistBySignedConf.SE(i) = 0;
    psyphys.RespDistBySignedConf.ticks.x = 0:9;
    psyphys.RespDistBySignedConf.ticks.y = [0 60];
    psyphys.RespDistBySignedConf.xlabel = 'Conf (highH to highF)';
    psyphys.RespDistBySignedConf.ylabel = 'Number of responses' ;
    psyphys.RespDistBySignedConf.xVals = 1:8;
    psyphys.RespDistBySignedConf.marker = 'o-';
end


%% confidence psychophysics plots - unsigned

for i = 1:4
    
    %percent correct for lowest to highest confidence
    psyphys.AccuracyByUnsignedConf.mean(i) = sum(ismember(cresp_, i) .* idx.cor) / sum(ismember(cresp_, i));
    psyphys.AccuracyByUnsignedConf.SE(i) = sqrt((psyphys.AccuracyByUnsignedConf.mean(i) * (1-(psyphys.AccuracyByUnsignedConf.mean(i)))) / sum(ismember(cresp_, i)));
    psyphys.AccuracyByUnsignedConf.ticks.x = 0:5;
    psyphys.AccuracyByUnsignedConf.ticks.y = [0 1];
    psyphys.AccuracyByUnsignedConf.xlabel = 'Confidence (low to high)';
    psyphys.AccuracyByUnsignedConf.ylabel = 'Percent Correct' ;
    psyphys.AccuracyByUnsignedConf.xVals = 1:4;
    psyphys.AccuracyByUnsignedConf.marker = 'o-';

    %RT for lowest to highest confidence
    psyphys.RTByUnsignedConf.mean(i) = nanmedian(rt(find((ismember(cresp_, i) .* idx.cor))));
    psyphys.RTByUnsignedConf.SE(i) = nanstd(rt(find((ismember(cresp_, i) .* idx.cor)))) / sqrt(sum(~isnan(rt(find((ismember(cresp_, i) .* idx.cor))))));
    psyphys.RTByUnsignedConf.ticks.x = 0:5;
    psyphys.RTByUnsignedConf.ticks.y = [0 4];
    psyphys.RTByUnsignedConf.xlabel = 'Confidence (low to high)';
    psyphys.RTByUnsignedConf.ylabel = 'RT (correct trials only)' ;
    psyphys.RTByUnsignedConf.xVals = 1:4;
    psyphys.RTByUnsignedConf.marker = 'o-';
    
                    %distribution of responses
    psyphys.RespDistByUnsignedConf.mean(i) = sum(ismember(cresp_, i) );
    psyphys.RespDistByUnsignedConf.SE(i) = 0;
    psyphys.RespDistByUnsignedConf.ticks.x = 0:5;
    psyphys.RespDistByUnsignedConf.ticks.y = [0 100];
    psyphys.RespDistByUnsignedConf.xlabel = 'Conf (low to high)';
    psyphys.RespDistByUnsignedConf.ylabel = 'Number of responses' ;
    psyphys.RespDistByUnsignedConf.xVals = 1:4;
    psyphys.RespDistByUnsignedConf.marker = 'o-';
end

%% coherence by pct correct psychophysics plots - unsigned

for i = 1:length(coh_set)
    iCoh = coh_set(i);
    
    psyphys.AccuracyByCoherence.mean(i) = sum(ismember(coh_, iCoh) .* idx.cor .* idx.cleanResp) / sum(ismember(coh_, iCoh) .* idx.cleanResp);
    psyphys.AccuracyByCoherence.SE(i) = sqrt((psyphys.AccuracyByCoherence.mean(i) * (1-(psyphys.AccuracyByCoherence.mean(i)))) / sum(ismember(coh_, iCoh).* idx.cleanResp) );
    psyphys.AccuracyByCoherence.ticks.x = 0:.1:1;
    psyphys.AccuracyByCoherence.ticks.y = [0 1];
    psyphys.AccuracyByCoherence.xlabel = 'Coherence';
    psyphys.AccuracyByCoherence.ylabel = 'Percent Correct' ;
    psyphys.AccuracyByCoherence.xVals = coh_set;
    psyphys.AccuracyByCoherence.marker = 'o-';

    psyphys.CorRTByCoherence.mean(i) = nanmedian(rt(find((ismember(coh_, iCoh) .* idx.cor .* idx.cleanResp))));
    psyphys.CorRTByCoherence.SE(i) = ste(rt(find((ismember(coh_, iCoh) .* idx.cor .* idx.cleanResp))));
    psyphys.CorRTByCoherence.ticks.x = 0:.1:1;
    psyphys.CorRTByCoherence.ticks.y = [0 4];
    psyphys.CorRTByCoherence.xlabel = 'Stimulus Coherence';
    psyphys.CorRTByCoherence.ylabel = 'RT (correct trials only)' ;
    psyphys.CorRTByCoherence.xVals = coh_set;
    psyphys.CorRTByCoherence.marker = 'o-';
    

end

%% coherence by pct correct psychophysics plots - signed

for i = 1:length(coh_set_signed)
    iCoh = coh_set_signed(i);
    
    psyphys.AccuracyBySignedCoherence.mean(i) = sum(ismember(coh_signed, iCoh) .* idx.cor .* idx.cleanResp) / sum(ismember(coh_signed, iCoh) .* idx.cleanResp);
    psyphys.AccuracyBySignedCoherence.SE(i) = sqrt((psyphys.AccuracyBySignedCoherence.mean(i) * (1-(psyphys.AccuracyBySignedCoherence.mean(i)))) / sum(ismember(coh_signed, iCoh).* idx.cleanResp) );
    psyphys.AccuracyBySignedCoherence.ticks.x = -1:.1:1;
    psyphys.AccuracyBySignedCoherence.ticks.y = [0 1];
    psyphys.AccuracyBySignedCoherence.xlabel = 'Coherence';
    psyphys.AccuracyBySignedCoherence.ylabel = 'Percent Correct' ;
    psyphys.AccuracyBySignedCoherence.xVals = coh_set_signed;
    psyphys.AccuracyBySignedCoherence.marker = 'o-' ;    

    psyphys.CorRTBySignedCoherence.mean(i) = nanmedian(rt(find((ismember(coh_signed, iCoh) .* idx.cor .* idx.cleanResp))));
    psyphys.CorRTBySignedCoherence.SE(i) = ste(rt(find((ismember(coh_signed, iCoh) .* idx.cor .* idx.cleanResp))));
    psyphys.CorRTBySignedCoherence.ticks.x = -1:.1:1;
    psyphys.CorRTBySignedCoherence.ticks.y = [0 4];
    psyphys.CorRTBySignedCoherence.xlabel = 'Stimulus Coherence';
    psyphys.CorRTBySignedCoherence.ylabel = 'RT (correct trials only)' ;
    psyphys.CorRTBySignedCoherence.xVals = coh_set_signed;
    psyphys.CorRTBySignedCoherence.marker = 'o-';
  
end

%% choice by coherence - signed

for i = 1:length(coh_set_signed)
    iCoh = coh_set_signed(i);
    
    psyphys.ChoiceByCoherence.mean(i) = sum(ismember(coh_signed, iCoh) .* idx.resp_face .* idx.cleanResp) / sum(ismember(coh_signed, iCoh) .* idx.cleanResp);
    psyphys.ChoiceByCoherence.SE(i) = sqrt((psyphys.ChoiceByCoherence.mean(i) * (1-(psyphys.ChoiceByCoherence.mean(i)))) / sum(ismember(coh_signed, iCoh).* idx.cleanResp) );
    psyphys.ChoiceByCoherence.ticks.x = -1:.1:1;
    psyphys.ChoiceByCoherence.ticks.y = [0 1];
    psyphys.ChoiceByCoherence.xlabel = 'Coherence';
    psyphys.ChoiceByCoherence.ylabel = 'Percent Face Responses' ;
    psyphys.ChoiceByCoherence.xVals = coh_set_signed;
    psyphys.ChoiceByCoherence.marker = 'o';
    
    g_coh_signed = -1:0.002:1;
    psyphys.ChoiceByCoherence.g_coh_signed = g_coh_signed;
    beta = logistfit([ones(length(coh_signed ),1) (coh_signed .* idx.cleanResp) (idx.resp_face .* idx.cleanResp)]);
    psyphys.ChoiceByCoherence.g_fc = 1./(1+exp(-(beta(1)+beta(2)*g_coh_signed)));
    
end

%% confidence by coherence - signed

for i = 1:length(coh_set_signed)
    iCoh = coh_set_signed(i);
    
    psyphys.SignedConfByCoherence.mean(i) = nanmean(cresp_(find((ismember(coh_signed, iCoh) .* idx.cor .* idx.cleanResp))));
    psyphys.SignedConfByCoherence.SE(i) = ste(cresp_(find((ismember(coh_signed, iCoh) .* idx.cor .* idx.cleanResp))));
    psyphys.SignedConfByCoherence.ticks.x = -1:.1:1;
    psyphys.SignedConfByCoherence.ticks.y = [0 4];
    psyphys.SignedConfByCoherence.xlabel = 'Coherence';
    psyphys.SignedConfByCoherence.ylabel = 'Confidence (Correct Trials)' ;
    psyphys.SignedConfByCoherence.xVals = coh_set_signed;
    psyphys.SignedConfByCoherence.marker = 'o-';
    

    
end

%% rt by coherence - signed

for i = 1:length(coh_set_signed)
    iCoh = coh_set_signed(i);
    
    psyphys.CorRTByCoherence.mean(i) = nanmean(rt(find((ismember(coh_signed, iCoh) .* idx.cor .* idx.cleanResp))));
    psyphys.CorRTByCoherence.SE(i) = ste(rt(find((ismember(coh_signed, iCoh) .* idx.cor .* idx.cleanResp))));
    psyphys.CorRTByCoherence.ticks.x = -1:.1:1;
    psyphys.CorRTByCoherence.ticks.y = [0 4];
    psyphys.CorRTByCoherence.xlabel = 'Coherence';
    psyphys.CorRTByCoherence.ylabel = 'RT (Correct Trials)' ;
    psyphys.CorRTByCoherence.xVals = coh_set_signed;
    psyphys.CorRTByCoherence.marker = 'o-';
        
end

%% stats
[B dev gRes.stats.pFaceRespByCoh]  = mnrfit(coh_signed, idx.resp_face+1);
gRes.stats.RTByConf.stats = regstats(rt, cresp_);
gRes.stats.RTByCoh.stats = regstats(rt, coh_);
gRes.stats.cohByConf.stats = regstats(coh_, cresp_);

[B dev gRes.stats.AccuracyByConf]  = mnrfit(cresp_, idx.cor+1);
[B dev gRes.stats.AccuracyByCoh]  = mnrfit(coh_, idx.cor+1);