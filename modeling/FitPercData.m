function output = FitPercData(dataFile)

% Loads perceptual task data, fits it to Diffusion Decision Models, 
% prints model fits/parameters and displays plots for each subject. 

% Load file containing perceptual task data
y = load(dataFile); 

% Store the totals for each parameter from all subjects,
% when the model is fitted based on the stimulus coherence data
cohMeanK = 0.0;
cohMeanB = 0.0;
cohMeanTND = 0.0;

% Store the totals for each parameter from all subjects,
% when the model is fitted based on the response certainty data
confMeanK = 0.0;
confMeanB = 0.0;
confMeanTND = 0.0;

% Responses made faster than 300 ms do not exceed the
% cutoff threshold, and are not considered in analysis
RTThresh = 0.3;
minBinSize = 10;
sortByConf = 1;
sortByCoh = 0;
fitType = 'dmat';
% Create a matrix to store the number of times that each level of response
% confidence was chosen by each of the subjects
confMat = zeros(length(y.PercBehData.sub), 8);

% Loop through each subject in the data file
for subj = find(~cellfun('isempty', y.PercBehData.sub))
    if ~isempty(y.PercBehData.sub)
        % Create a logical matrix indicating which trials contain valid
        % response time data (i.e. not NaN, and exceeds the RT threshold).
        % These trial indices are the only ones we will consider in our analysis.
        indexRealRTs = (y.PercBehData.sub{subj}.rt > RTThresh)';
        indexLegitResp = (y.PercBehData.sub{subj}.resp_face + y.PercBehData.sub{subj}.resp_house)==1;
        
        indexToInclude = indexRealRTs' & indexLegitResp;
        
        % Include valid response times in the trialData.RT field
        trialData.RT = y.PercBehData.sub{subj}.rt(indexToInclude)';
        
        % coh_signed is signed from -1 (full coherence house) to 1 (full coherence face);
        % the absolute value is taken, so trialData.Coh ranges only from 0-1
        trialData.Coh = abs(y.PercBehData.sub{subj}.coh_signed(indexToInclude))';
        
        % 0 = trial where subject responded "house,"
        % 1 = trial where subject responded "face."
        trialData.Resp = y.PercBehData.sub{subj}.resp_face(indexToInclude)';
        
        % 0 = trial with a house stimulus, 1 = trial with a face stimulus
        trialData.Stim = y.PercBehData.sub{subj}.face(indexToInclude)';
        
        % Create a logical matrix indicating which trials contained a stimulus
        % with a 0 coherence value (neither a face nor a house)
        %indexNoCoherence = (trialData.Coh == 0);
        
        % Because the 0-coherence trials represent neither a face nor a house,
        % the stimuli for those trials are divided evenly between the two
        % conditions (still trying to figure out how to do this)
        %         conditionVector = zeros(1, sum(indexNoCoherence))';
        %         conditionVector(1:(sum(indexNoCoherence) / 2)) = 1;
        %         randVector = rand(1, length(conditionVector))';
        %         [sortedRandVector, indices] = sort(randVector);
        %         conditionVectorSorted = conditionVector(indices);
        %         trialData.Stim(indexNoCoherence) = conditionVectorSorted;
        
        % Diffusion model parameters
        k = .5;
        B = 1;
        tnd = .3 ;
        
        if sortByCoh
            if strcmp(fitType, 'dmat')
                data = [trialData.Coh', trialData.Resp'==trialData.Stim', trialData.RT'];
                options = multiestv4;
                nCohs = length(unique(trialData.Coh));
                O = ones((nCohs),1);
                I = eye(nCohs);
                design_matrix = {O,O,O,O,O,O,I,O,O};
                options.DesignMatrix = design_matrix;
                output.sub{subj}.mod = multiestv4(data,options);
                output.sub{subj}.data = multiestv4(data,options);
                qpplot(data, [10:20:90], output.sub{subj}.mod.Minimum)
            elseif strcmp(fitType, 'roozbeh')
                % Fit the diffusion model to the stimulus coherence data
                % and print the parameters
                [modelParam, modelLL, trialData] = AGfitDiffusion(trialData);
                fprintf('\n\nModel parameters for subject %g (stimulus coherence data):\n', subj);
                fprintf('\tk = %g, actual %g\n\tB = %g, actual %g\n\ttnd = %g, actual %g\n\n', ...
                    modelParam.final(1), k, modelParam.final(2), B, modelParam.final(3), tnd);
                
                % Add the coherence model parameters from this participant to the total sum of
                % each parameter for all participants, from which the means will be found
                cohMeanK = cohMeanK + modelParam.final(1);
                cohMeanB = cohMeanB + modelParam.final(2);
                cohMeanTND = cohMeanTND + modelParam.final(3);
                
                % Create a set of the coherence levels that were used, by finding the
                % unique values amongst the coherence ratings from all trials
                coh_set = unique(trialData.Coh);
                
                % What are the mean probability correct and reaction time values
                % for each level of stimulus coherence?
                [pc, pc_se] = calcGroupMean(trialData.Resp==trialData.Stim, trialData.Coh, coh_set, 'binary');
                [rt, rt_se] = calcGroupMean(trialData.RT, trialData.Coh, coh_set);
                
                % Create the smooth fit curves for the stimulus coherence data
                g_coh = (0:1e-4:1);
                D.Coh = g_coh;
                D.Stim = nan(size(g_coh));
                D.Resp = nan(size(g_coh));
                D.RT = nan(size(g_coh));
                [~,~,D] = AGfitDiffusion(D, modelParam.final, ones(size(modelParam.final)));
                g_pc = D.expectedPC;
                g_rt = D.expectedRT;
                
                % Plot the probability correct data and the fitted curves
                % against the varying levels of stimulus coherence
                figure;
                hold on;
                Ic0 = coh_set==0;
                errorbar(coh_set(~Ic0), pc(~Ic0), pc_se(~Ic0), 'ok', 'MarkerFaceColor', 'k');
                plot(g_coh, g_pc, 'k');
                xlabel('Stimulus strength');
                ylabel('Probability correct');
                
                % Plot the reaction time data and the fitted curves
                % against the varying levels of stimulus coherence
                figure;
                hold on;
                errorbar(coh_set, rt, rt_se, 'ok', 'MarkerFaceColor', 'k');
                plot(g_coh, g_rt, 'k');
                xlabel('Stimulus strength');
                ylabel('Reaction time');
            end
        end
        
        if sortByConf
            conf = y.PercBehData.sub{subj}.conf_unsigned(indexRealRTs)';
            conf_set_h = unique(conf);
            binCount = hist(conf, conf_set_h);
            goodBins = binCount>minBinSize;
            idxToInclude2 = ismember(conf, conf_set_h(goodBins));
            
            trialData.Coh = (conf(idxToInclude2)-.5) / 4;
            trialData.RT = trialData.RT(idxToInclude2);
            trialData.Stim = trialData.Stim(idxToInclude2);
            trialData.Resp = trialData.Resp(idxToInclude2);
            
            if strcmp(fitType, 'dmat')
                data = [trialData.Coh', trialData.Resp'==trialData.Stim', trialData.RT'];
                options = multiestv4;
                nCohs = length(unique(trialData.Coh));
                O = ones((nCohs),1);
                I = eye(nCohs);
                design_matrix = {O,O,O,O,O,O,I,O,O};
                options.DesignMatrix = design_matrix;
                output.sub{subj}.mod = multiestv4(data,options);
                output.sub{subj}.data = multiestv4(data,options);
                qpplot(data, [10:20:90], output.sub{subj}.mod.Minimum)
            elseif strcmp(fitType, 'roozbeh')
                trialData.RT = trialData.RT*1000;
                
                % Reassign trialData.Coh to contain the data for this subject's
                % level of confidence in their responseson each trial,
                % so that the model can be re-fit

                
                conf_set = unique(trialData.Coh);
                
                % Fit the diffusion model to the response confidence data
                % and print the parameters
                [modelParam, modelLL, trialData] = AGfitDiffusion(trialData);
                fprintf('\n\nModel parameters for subject %g (response certainty data):\n', subj);
                fprintf('\tk = %g, actual %g\n\tB = %g, actual %g\n\ttnd = %g, actual %g\n\n', ...
                    modelParam.final(1), k, modelParam.final(2), B, modelParam.final(3), tnd);
                
                % Add the model parameters from this participant to the total sum of
                % each parameter for all participants, from which the means will be found
                model(subj) = modelParam;
                
                % What are the mean probability correct and reaction time values
                % for each level of response confidence?
                [pc, pc_se] = calcGroupMean(trialData.Resp==trialData.Stim, trialData.Coh, conf_set);
                [rt, rt_se] = calcGroupMean(trialData.RT, trialData.Coh, conf_set);
                
                % Create the smooth fit curves
                g_coh = (0:1e-4:1);
                D.Coh = g_coh;
                D.Stim = nan(size(g_coh));
                D.Resp = nan(size(g_coh));
                D.RT = nan(size(g_coh));
                [~,~,D] = AGfitDiffusion(D, modelParam.final, ones(size(modelParam.final)));
                g_pc = D.expectedPC;
                g_rt = D.expectedRT;
                
                % Plot the probability correct data and the fitted curves
                % against the varying levels of response confidence
                figure;
                hold on;
                Ic0 = conf_set==0;
                errorbar(conf_set(~Ic0), pc(~Ic0), pc_se(~Ic0), 'ok', 'MarkerFaceColor', 'k');
                plot(g_coh, g_pc, 'k');
                xlabel('Confidence level of response');
                ylabel('Probability correct');
                
                % Plot the reaction time data and the fitted curves
                % against the varying levels of response confidence
                figure;
                hold on;
                errorbar(conf_set, rt, rt_se, 'ok', 'MarkerFaceColor', 'k');
                plot(g_coh, g_rt, 'k');
                xlabel('Confidence level of response');
                ylabel('Reaction time');
                
            end
        end
        % Clear trialData structure for the next subject
        clear trialData;
    end
end

output.mod = model;
% use a different non-linear transform function, one that increases
% monotonically.

% exclude bins with too few data points.