function res = PM_classificationExportToCSV(par, qqq)
  
             
cd /Users/alangordon/mounts/wagner5/alan/perceptMnemonic/fmri_data/pm_100910/behav/

fileNames = {'10-09-2010_Yula_cert01.mat'
    '10-09-2010_Yula_cert02.mat'
    '10-09-2010_Yula_cert03.mat'
    '10-09-2010_Yula_cert04.mat'
    '10-09-2010_Yula_cert05.mat'
    '10-09-2010_Yula_cert06.mat'
    '10-09-2010_Yula_cert07.mat'};

countdown = 12; %length of countdown period (in secs) at the beginning of each run 
         
         
% responses = {[4 4 4 3 3 4 4 4 4 4 3 3 4 3 4 4 3 3 4]
%              [4 3 4 4 3 3 3 4 3 3 3 4 3 4 3 3 4 3 4]
%              [4 4 4 3 3 3 3 3 3 4 3 3 4 3 3 4 3 3 4 4]
%              [4 3 4 3 4 4 3 4 4 4 4 4 4 3 4 4 4 3 3]
%              [4 3 4 3 3 3 4 3 3 3 3 4 4 4 4 4 4 4 4 4]
%              [4 4 4 4 4 4 3 4 4 3 4 4 4 3 4 4 3 4 4 4]
%              [3 3 3 3 4 3 3 3 3 3 4 4 4 3 3 4 3 4 3 3]
%              [3 4 4 4 3 3 3 3 4 3 4 3 4 3 3 3 4 3 3]
%              [3 3 4 4 4 3 3 4 3 3 4 4 3 4 4 4 3 4 4 3]
%              [4 4 4 3 3 4 4 4 4 3 4 4 3 4 3 4 3 4 3]};
% for i = 1 : length(fileNames)
%     load(fileNames{i});
%     for j = 1 : length(trial_data)
%         if (trial_data(j).stim_group==1 && responses{i}(j)==4) || ...
%            (trial_data(j).stim_group==2 && responses{i}(j)==3)
%             trial_data(j).result = {'CORRECT'}; 
%             trial_data(j).response = trial_data(j).stim_group;
%         else
%             trial_data(j).result = {'WRONG'};
%             trial_data(j).response = 3-trial_data(j).stim_group;
%         end
%     end
%     save(fileNames{i}, 'eyelink_struct', 'screen_struct', 'task_struct', 'trial_data');
% end



%********************** XXXXXXXXXXXXXXX    
        %load behavioral data files and generate 
trial_data = combineDataFile(fileNames, '../RawData');
trial_num = length(trial_data);


    %extract stim presentation and behavioral variables
if isfield(trial_data,'response')
    resp_ = cat(1, trial_data.response);
else
    resp_ = nan(trial_num, 1);
end
if isfield(trial_data, 'cresponse')
    cresp_ = cat(1, trial_data.cresponse);
else
    cresp_ = nan(trial_num, 1);
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
resp_ = resp_(valid_trials);
cresp_ = cresp_(valid_trials);
trial_data = trial_data(valid_trials);
time0 = cat(1, trial_data.time0);
start_t = cat(1, trial_data.start_t);
event_name = cat(1, trial_data.event_name);
event_time = cat(1, trial_data.event_time);
event_time = event_time + repmat( start_t-time0-countdown ,[1 size(event_time,2)]);
stim_on = event_time(strmatch('stim_on',event_name));
stim_off = event_time(strmatch('stim_off',event_name));
dur_ = stim_off - stim_on;
stim_ = cat(1, trial_data.stim_group);
coh_ = cat(1, trial_data.stim_coh_seq);
scan_ = cat(1, trial_data.ownership);       %which scan each trial belongs to
coh_set = unique(coh_);

%%
%alltrials =  (scan_-1)*318 + stim_on; 
alltrials =  (scan_-1).*par.numvols(scan_)' * par.TR + stim_on;
coh_set_pct = coh_set*100;



idx.face = stim_==1;
idx.house = stim_==2;

idx.resp_face = ismember(resp_ , 5:8);
idx.resp_house = ismember(resp_ , 1:4);

idx.cor = idx.face .* idx.resp_face + idx.house .* idx.resp_house;
idx.inc = idx.face .* idx.resp_house + idx.house .* idx.resp_face;

idx.all  = idx.cor + idx.inc;

perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
resp_set = {'resp_face' 'resp_house'};
conf_set = 1:4;


res(:,1) = ['decidCor'; num2cell(1*idx.cor + 0* idx.inc -1*(~idx.cor.*~idx.inc))];
res(:,2) = ['coh'; num2cell(coh_)];
res(:,3) = ['stimType'; num2cell(stim_)];
res(:,4) = ['conf'; num2cell(1+mod(resp_-1,4))];
res(:,5) = ['RT'; num2cell([trial_data.rt]')];
res(:,6) = ['act'; num2cell(qqq.subj{1}.iter{1}.iterations.acts(1,:)')];
res(:,7) = ['classCor'; num2cell(double(qqq.subj{1}.iter{1}.iterations(1).perfmet.corrects'))]; 

cell2csv('/Users/alangordon/mounts/wagner5/alan/perceptMnemonic/fmri_data/pm_100910/classifierAnalysis/pm_100910.csv', res, ',', 2000);


