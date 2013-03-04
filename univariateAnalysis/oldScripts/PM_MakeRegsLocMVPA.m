function PM_MakeRegsLocMVPA(par)


cd (par.behavdir)

dFN = dir('*Loc*');
fileNames = {dFN.name}';

% fileNames = {'12-09-2010_Marc_Loc1.mat'
%              '12-09-2010_Marc_Loc2.mat' 
%              '12-09-2010_Marc_Loc3.mat'};


countdown = 12; %length of countdown period (in secs) at the beginning of each run 
         
[res idx] = fMRIBehAnalysis_Loc(par);        
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
%         %load behavioral data files and generate 
% trial_data = combineDataFile(fileNames, par.behavdir);
% trial_num = length(trial_data);
% 
% 
%     %extract stim presentation and behavioral variables
% if isfield(trial_data,'response')
%     resp_ = cat(1, trial_data.response);
% else
%     resp_ = nan(trial_num, 1);
% end
% if isfield(trial_data, 'cresponse')
%     cresp_ = cat(1, trial_data.cresponse);
% else
%     cresp_ = nan(trial_num, 1);
% end
% 
% % if all(isnan(resp_))
% %         %accept all the trials if we did not collect response
% %     valid_trials = 1 : trial_num;
% % elseif ~all(isnan(resp_)) && all(isnan(cresp_))
% %         %if we collected direction choices but not certainty responses
% %     valid_trials = find(~isnan(resp_));
% % else
% %         %if we collected both direction choices and certainty responses 
% %     valid_trials = find(~isnan(resp_) & ~isnan(cresp_));
% % end
% 
% valid_trials = 1 : trial_num;
% 
% 
% resp_ = resp_(valid_trials);
% cresp_ = cresp_(valid_trials);
% trial_data = trial_data(valid_trials);
% time0 = cat(1, trial_data.time0);
% start_t = cat(1, trial_data.start_t);
% event_name = cat(1, trial_data.event_name);
% event_time = cat(1, trial_data.event_time);
% event_time = event_time + repmat( start_t-time0-countdown ,[1 size(event_time,2)]);
% stim_on = event_time(strmatch('stim_on',event_name));
% stim_off = event_time(strmatch('stim_off',event_name));
% dur_ = stim_off - stim_on;
% stim_ = cat(1, trial_data.stim_group);
% coh_ = cat(1, trial_data.stim_coh_seq);
% scan_ = cat(1, trial_data.ownership);       %which scan each trial belongs to
% coh_set = unique(coh_);
% cor = cat(1,trial_data.result);
% 
% alltrials =  (scan_-1).*par.numvols(scan_)' * par.TR + stim_on; 
% coh_set_pct = coh_set*100;
% 
% 
% idx.cor = strcmp(cor, 'CORRECT');
% idx.inc = ~idx.cor;
% 
% idx.face = stim_==1;
% idx.house = stim_==2;
% 
% 
% idx.faceCor = idx.cor .* idx.face;
% idx.houseCor = idx.cor .* idx.house;
% 
% %idx.faceCor = 
% %idx.houseCor = 


perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
resp_set = {'respFace' 'respHouse'};


i = 0;

for l=1:length(class_set)
          
    idx.thisOns = idx.(class_set{l}) .* idx.coh;
    
    if ~isempty(idx.alltrials(find(idx.thisOns)))
        
        i = i+1;
        fName = sprintf('%s', class_set{l});
        

        
        stimOnsets{i}= idx.alltrials(find(idx.thisOns))';
        stimNames{i} = fName;
        stimDurations{i} = 0; 
    end
end

if ~isempty(idx.alltrials(find(~idx.coh)))
    i = i+1;
    stimOnsets{i} = idx.alltrials(find(~idx.coh))';
    stimNames{i} = 'noise';
    stimDurations{i} = 0;
end

% if ~isempty(idx.alltrials(find(idx.inc)))
%     i = i+1;
%     stimOnsets{i} = idx.alltrials(find(idx.inc))';
%     stimNames{i} = 'junk';
%     stimDurations{i} = 0;
% end

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
save mvpa_ons.mat onsets durations names;
%save regs.mat R