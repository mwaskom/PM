function output = FitMnemData(dataFile)

% Loads mnemonic task data, fits it to Diffusion Decision Models, 
% prints model fits/parameters and displays plots for each subject. 

% Load file containing mnemonic task data
y = load(dataFile); 

% Store the totals for each parameter from all subjects
meanK = 0.0; 
meanB = 0.0; 
meanTND = 0.0; 

% Responses made faster than 300 ms do not exceed the 
% cutoff threshold, and are not considered in analysis
RTThresh = 0.3; 
minBinSize = 10;
% Stores the total number of subjects
numSubs = 0; 

% Loop through each subject in the data file
for subj = find(~cellfun('isempty', y.MnemBehData.sub))
    
    % Some cells in the data structure contained invalid data and have been
    % cleared. We do not wish to include these cells in our analysis, so we
    % only proceed if the cell actually contains valid data from a subject.
    if length(y.MnemBehData.sub{subj}) > 0
        
        % Increase the subject count by 1
        numSubs = numSubs + 1;
        
        conf = y.MnemBehData.sub{subj}.cresp;
        conf_set_h = unique(conf);
        binCount = hist(conf, conf_set_h);
        goodBins = binCount>minBinSize;
        idxBinLargeEnough = ismember(conf, conf_set_h(goodBins));
        
    % Create a logical matrix indicating which trials contain valid
    % response time data (i.e. not NaN, and exceeds the RT threshold).
    % These trial indices are the only ones we will consider in our analysis. 
    idxRealRTs = (y.MnemBehData.sub{subj}.rt > RTThresh);
    idxToInclude = idxRealRTs & idxBinLargeEnough;
    % Include valid response times in the trialData.RT field
    trialData.RT = 1000*y.MnemBehData.sub{subj}.rt(idxToInclude);
    
    % Top lines bin subject responses into "high certainty" (both 1: high
    % confidence house, and 8: high confidence face) and "low certainty"
    % (all other intermediate confidence levels).   
    % The line at the bottom bins subject responses into 4
    % categories that indicate unsigned gradations of confidence (i.e. both
    % 1 and 8 become a 4, 2 and 7 are a 3, and 4 and 5 are both 1). 
    %indexHighCertHouse = (y.MnemBehData.sub{subj}.certRat(indexRealRTs) == 1);
    %indexHighCertFace = (y.MnemBehData.sub{subj}.certRat(indexRealRTs) == 8); 
    %indexHighCert = (indexHighCertHouse + indexHighCertFace); 
    %trialData.Coh = indexHighCert;
    trialData.Coh = (conf(idxToInclude)-.5) / 4; 
    
    % 0 = trial where subject responded "house,"
    % 1 = trial where subject responded "face." 
    trialData.Resp = y.MnemBehData.sub{subj}.respFace(idxToInclude); 
    
    % 0 = trial with a house stimulus, 1 = trial with a face stimulus
    trialData.Stim = y.MnemBehData.sub{subj}.face(idxToInclude);
    
    % Diffusion model parameters
    k = 1; 
    B = 1.5;
    tnd = .3; 
    
    % Fit the diffusion model and print the parameters  
    [modelParam, modelLL, trialData] = AGfitDiffusion(trialData);
    fprintf('\n\nmodel parameters:\n\tk = %g, actual %g\n\tB = %g, actual %g\n\ttnd = %g, actual %g\n\n', ...
            modelParam.final(1), k, modelParam.final(2), B, modelParam.final(3), tnd);
        
    % Add the model parameters from this participant to the total sum of
    % each parameter for all participants, from which the means will be found
    model(subj) = modelParam;
    
    % Create a set of the certainty levels by finding the unique
    % values amongst the participant certainty responses from all trials
    cert_set = unique(trialData.Coh);     
   
    % What are the mean probability correct and reaction time values
    % for each level of stimulus coherence? 
    [pc, pc_se] = calcGroupMean(trialData.Resp==trialData.Stim, trialData.Coh, cert_set, 'binary'); 
    [rt, rt_se] = calcGroupMean(trialData.RT, trialData.Coh, cert_set); 
    
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
    % against the varying levels of response certainty
    figure; 
    hold on; 
    Ic0 = cert_set==0; 
    errorbar(cert_set(~Ic0), pc(~Ic0), pc_se(~Ic0), 'ok', 'MarkerFaceColor', 'k'); 
    plot(g_coh, g_pc, 'k');
    xlabel('Response certainty');
    ylabel('Probability correct');
    
    % Plot the reaction time data and the fitted curves 
    % against the varying levels of response certainty
    figure;
    hold on;
    errorbar(cert_set, rt, rt_se, 'ok', 'MarkerFaceColor', 'k'); 
    plot(g_coh, g_rt, 'k');
    xlabel('Response certainty');
    ylabel('Reaction time');
    
    % Clear trialData structure for the next subject
    clear trialData; 
    
    end
    
end

output.mod = model;

end
