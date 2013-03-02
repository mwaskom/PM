function [res psyphys idx] = behavAnalysisPerc(par)


cd (par.behavdir)
dFN = dir('*Perc*');
fileNames = {dFN.name};

localizer = 0;
countdown = 12;

trial_data = combineDataFile(fileNames, par.behavdir);
trial_num = length(trial_data);


%extract stim presentation and behavioral variables
if isfield(trial_data,'response')
    resp = cat(1, trial_data.response);
    
    if strcmp(par.substr, 'pm_031711')
        resp(isnan(resp))=9; %'5' button presses were not recorded for this subject.  Assume that Nans reflect '9' button presses.
        resp(resp>4) = resp(resp>4)-1; %change '6 7 8 9' to '5 6 7 8';
    end

    resp_ = 3-ceil(resp/4);
    cresp_ = mod(resp-1,4)+1;
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

%********************** XXXXXXXXXXXXXXX    



    %probability correct
[pc, pc_se] = calcGroupMean(resp_==stim_, coh_, coh_set, 'binary');
pc(coh_set==0) = NaN;
    %probability face
[pf, pf_se] = calcGroupMean(resp_==1, coh_signed, coh_set_signed, 'binary');
    %certainty
[cert, cert_se] = calcGroupMean(cresp_, coh_, coh_set);
[cert_signed, cert_se_signed] = calcGroupMean(cresp_, coh_signed, coh_set_signed);
    %certainty for correct responses
L = (resp_==stim_) | (coh_==0);
[cert_cor, cert_cor_se] = calcGroupMean(cresp_(L), coh_(L), coh_set);
[cert_cor_signed, cert_cor_se_signed] = calcGroupMean(cresp_(L), coh_signed(L), coh_set_signed);
L = (resp_~=stim_) | (coh_==0);
[cert_wrg, cert_wrg_se] = calcGroupMean(cresp_(L), coh_(L), coh_set);
[cert_wrg_signed, cert_wrg_se_signed] = calcGroupMean(cresp_(L), coh_signed(L), coh_set_signed);

logscale = 0;

if any(~isnan(pc))
    g_coh = 1e-4*2.^(0:0.1:14);
    beta = logistfit([ones(length(coh_),1) coh_ resp_==stim_]);
    g_pc = 1./(1+exp(-(beta(1)+beta(2)*g_coh)));
    figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
    hold on;
    if logscale
        L = coh_set~=0;
        errorbar(log(coh_set(L)), pc(L), pc_se(L), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        plot(log(g_coh), g_pc, 'r');
        set(gca, 'XLim', log([0.9 60]*1e-2), 'XTick', log([1 2.5 5 10 20 40]*1e-2), 'XTickLabel', [0 2.5 5 10 20 40], 'TickDir', 'out');
    else
        errorbar(coh_set, pc, pc_se, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        plot(g_coh, g_pc, 'r');
        set(gca, 'XLim', [-0.001 1.001], 'XTick', 0:0.1:1, 'XTickLabel', makeTickLabel(0:0.1:1,0.2), 'TickDir', 'out');
    end
    xlabel('Stimulus strength (%coh)');
    ylabel('Probability correct');
    box off;
    
    
    g_coh_signed = -1:0.002:1;
    beta = logistfit([ones(length(coh_signed),1) coh_signed resp_==1]);
    g_fc = 1./(1+exp(-(beta(1)+beta(2)*g_coh_signed)));
    figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
    hold on;
    errorbar(coh_set_signed, pf, pf_se, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    plot(g_coh_signed, g_fc, 'r');
    set(gca, 'XLim', [-1.001 1.001], 'XTick', -1:0.1:1, 'XTickLabel', makeTickLabel(-1:0.1:1,0.2), 'TickDir', 'out');
    xlabel('Stimulus strength (%coh)');
    ylabel('Probability face response');
    grid on;
end
if any(~isnan(cert))
    C = coh_set;
    C(1) = 0.01;
    figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
    hold on;
    if logscale
%         errorbar(log(C), cert_cor, cert_cor_se, 'o-', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
%         L = cert_wrg_se~=0 & ~isnan(cert_wrg_se);
%         errorbar(log(C(L)), cert_wrg(L), cert_wrg_se(L), 'o-', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
        errorbar(log(C), cert, cert_se, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        set(gca, 'XLim', log([0.9 110]*1e-2), 'XTick', log([1 2.5 5 10 20 40 80]*1e-2), 'XTickLabel', [0 2.5 5 10 20 40 80], ...
                 'YLim', [1 5], 'TickDir', 'out');
    else
%         errorbar(coh_set, cert_cor, cert_cor_se, 'o-', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
%         L = cert_wrg_se~=0 & ~isnan(cert_wrg_se);
%         errorbar(coh_set(L), cert_wrg(L), cert_wrg_se(L), 'o-', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
        errorbar(coh_set, cert, cert_se, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        set(gca, 'XLim', [-0.001 1.001], 'XTick', 0:0.1:1, 'XTickLabel', makeTickLabel(0:0.1:1,0.2), ...
                 'YLim', [1 5], 'TickDir', 'out');
    end
    xlabel('Stimulus strength (%coh)');
    ylabel('Certainty');
%     h = legend({'correct','error'}, 'Location', 'NorthWest');
%     set(h, 'FontSize', 9);
    box off;
    
    figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
    hold on;
    errorbar(coh_set_signed, cert_cor_signed, cert_cor_se_signed, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    set(gca, 'XLim', [-1.001 1.001], 'XTick', -1:0.1:1, 'XTickLabel', makeTickLabel(-1:0.1:1,0.2), 'TickDir', 'out');
    xlabel('Stimulus strength (%coh)');
    ylabel('Certainty');
    grid on;    
end


%%

cor_idx = find(resp_==stim_);
inc_idx = find((resp_~=stim_) & ~isnan(resp_));

[rt_signed, rt_signed_se] = calcGroupMean(rt(cor_idx), coh_signed(cor_idx), coh_set_signed);
[rt_cor, rt_cor_se] = calcGroupMean(rt(cor_idx), coh_(cor_idx), coh_set);
[rt_inc, rt_inc_se] = calcGroupMean(rt(inc_idx), coh_(inc_idx), coh_set);


figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
errorbar(coh_set_signed, rt_signed, rt_signed_se, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
set(gca, 'XLim', [-1.001 1.001], 'XTick', -1:0.1:1, 'XTickLabel', makeTickLabel(-1:0.1:1,0.2), ...
    'YLim', [0 3], 'TickDir', 'out');
grid on;
% hold on
% 
% errorbar(coh_set, rt_inc, rt_inc_se, 'o-', 'Color', 'r',  'MarkerSize', 4, 'LineStyle', 'none');
% set(gca, 'XLim', [-0.001 1.001], 'XTick', 0:0.1:1, 'XTickLabel', makeTickLabel(0:0.1:1,0.2), ...
%     'YLim', [0 2.5], 'TickDir', 'out')

xlabel('Stimulus strength (%coh)');
ylabel('RT');


%%

% cor_idx = find(resp_==stim_);
% inc_idx = find((resp_~=stim_) & ~isnan(resp_));
% [rt_cor, rt_cor_se] = calcGroupMean(rt(cor_idx), coh_(cor_idx), coh_set);
% [rt_inc, rt_inc_se] = calcGroupMean(rt(inc_idx), coh_(inc_idx), coh_set);
% 
% 
% figure('Color', 'w', 'Position', [100 100 300 250], 'PaperPositionMode', 'auto');
% errorbar(coh_set_signed, rt_cor, rt_cor_se, 'o-', 'Color', 'b',  'MarkerSize', 4, 'LineStyle', 'none');
% set(gca, 'XLim', [-0.001 1.001], 'XTick', 0:0.1:1, 'XTickLabel', makeTickLabel(0:0.1:1,0.2), ...
%     'YLim', [0 3], 'TickDir', 'out');
% 
% % hold on
% % 
% % errorbar(coh_set, rt_inc, rt_inc_se, 'o-', 'Color', 'r',  'MarkerSize', 4, 'LineStyle', 'none');
% % set(gca, 'XLim', [-0.001 1.001], 'XTick', 0:0.1:1, 'XTickLabel', makeTickLabel(0:0.1:1,0.2), ...
% %     'YLim', [0 2.5], 'TickDir', 'out')
% 
% xlabel('Stimulus strength (%coh)');
% ylabel('RT');
