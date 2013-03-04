function PM_MakeRegs_PercMvpa(par)

cd (par.behavdir);

dFN = dir('*Perc*');
fileNames = {dFN.name}';

countdown = 12;

[res psy idx] = Perceptual_fMRIBehAnalysis(par);
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
% if all(isnan(resp_))
%         %accept all the trials if we did not collect response
%     valid_trials = 1 : trial_num;
% elseif ~all(isnan(resp_)) && all(isnan(cresp_))
%         %if we collected direction choices but not certainty responses
%     valid_trials = find(~isnan(resp_));
% else
%         %if we collected both direction choices and certainty responses 
%     valid_trials = find(~isnan(resp_) & ~isnan(cresp_));
% end
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
% 
% 
% %alltrials =  (scan_-1)*318 + stim_on; 
% alltrials =  (scan_-1).*par.numvols(scan_)' * par.TR + stim_on;
% coh_set_pct = coh_set*100;
% 
% 
% % idx.cor = (resp_ == stim_ );
% % idx.inc = ~(resp_ == stim_ );
% 
% idx.face = stim_==1;
% idx.house = stim_==2;
% 
% idx.resp_face = ismember(resp_ , 5:8);
% idx.resp_house = ismember(resp_ , 1:4);
% 
% idx.cor = idx.face .* idx.resp_face + idx.house .* idx.resp_house;
% idx.inc = idx.face .* idx.resp_house + idx.house .* idx.resp_face;
% 
% idx.all  = idx.cor + idx.inc;

perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
resp_set = {'resp_face' 'resp_house'};
conf_set = 1:4;


i = 0;
for c=1:length(idx.coh_set_pct)
    
    if idx.coh_set(c)~=0
        for l=1:length(class_set)
            
            
            
            for p=1:length(perf_set)
                
                
                idx.thisOns = (idx.coh_==idx.coh_set(c)) .* idx.(perf_set{p}) .* idx.(class_set{l});
                
                if ~isempty(idx.alltrials(find(idx.thisOns)))
                    
                    i = i+1;
                    fName = sprintf('%s_coh%s_%s', class_set{l}, num2str(idx.coh_set_pct(c)), perf_set{p});
                    
                    
                    
                    %ons.(fName) =
                    
                    stimOnsets{i}= idx.alltrials(find(idx.thisOns));
                    stimNames{i} = fName;
                    stimDurations{i} = 0;
                    stimRTs{i} = idx.rt(find(idx.thisOns));
                end
                
            end
            
            
        end
        
    end
    
end





% for f=1:length(conf_set)
%     
%     idx.thisConf = ismember(resp_, [f, f+4]);
%     %idx.thisOns = (coh_==coh_set(c)) .* idx.thisConf .* idx.(class_set{l});
%     
%     if ~isempty(alltrials(find(idx.thisConf)))
%         
%         i = i+1;
%         fName = sprintf('conf%g',  conf_set(f));
%         
%         
%         
%         %ons.(fName) =
%         
%         stimOnsets{i}= alltrials(find(idx.thisConf));
%         stimNames{i} = fName;
%         stimDurations{i} = 0;
%         
%     end
%     
% end






for r=1:length(resp_set)
    idx.thisOns = (idx.coh_==0) .* idx.(resp_set{r}) ;
    if ~isempty(idx.alltrials(find(idx.thisOns)))
        i = i+1;
        fName  = sprintf('coh0_%s', resp_set{r});
        stimOnsets{i}= idx.alltrials(find(idx.thisOns));
        stimNames{i} = fName;
        stimDurations{i} = 0;
        stimRTs{i} = idx.rt(find(idx.thisOns));
    end
end

for r=1:length(resp_set)
    idx.thisOns = idx.(resp_set{r}) ;
    if ~isempty(idx.alltrials(find(idx.thisOns)))
        i = i+1;
        fName  = sprintf('%s', resp_set{r});
        stimOnsets{i}= idx.alltrials(find(idx.thisOns));
        stimNames{i} = fName;
        stimDurations{i} = 0;
        stimRTs{i} = idx.rt(find(idx.thisOns));
    end
end

            

for l=1:length(class_set)
     for p=1:length(perf_set)       
            idx.thisOns = (idx.coh_~=0) .* idx.(perf_set{p}) .* idx.(class_set{l});
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                i = i+1;
                fName = sprintf('%s_%s', class_set{l}, perf_set{p});
                stimOnsets{i}= idx.alltrials(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = 0;
                stimRTs{i} = idx.rt(find(idx.thisOns));
            end
     end
end



onsets = stimOnsets;
names = stimNames;
durations = stimDurations;
rts = stimRTs;

if ~exist(par.analysisdir)
    mkdir (par.analysisdir);
end


sessReg = zeros(sum(par.numvols),length(par.numvols) - 1);
for i = 1:(length(par.numvols) - 1)
    sessReg((i-1)*par.numvols(i)+1:i*par.numvols(i),i) = ones(par.numvols(i),1);
end

R = horzcat(sessReg);



% for i = 1:par.numscans
%     cd(fullfile(par.subdir, 'functional', ['scan' (prepend(num2str(par.scans_to_include(i)),2))]   ));
%     motTxt = dir('rp*');
%     motRegs_h{i} = textread(motTxt.name);
% end
% 
% motRegs = vertcat(motRegs_h{:});

% R = horzcat(sessReg, motRegs);




cd (par.analysisdir);
save mvpa_ons.mat onsets durations names rts;
%save regs.mat R