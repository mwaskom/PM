
function [modelParam, modelLL, trialData, exitflag] = AGfitDiffusion(trialData, guess, fixed, fitwhat) %another argument for task? 
% 
% [modelParam, modelLL, trialData] = fitDiffusion(trialData, guess, fitwhat)
% 

% 
% RK, 7/8/2010
% 

if ~exist('guess')
    guess = [0.5 20 300 10 0.5];
end

if ~exist('fixed')
    fixed = zeros(size(guess));
end

if ~exist('fitwhat')
    fitwhat = 'BOTH';
elseif strcmp(fitwhat,'PC')
    fixed(end) = 1;
end
%h = figure;
%ax = subplot(1,1,1);


coh_set = unique(trialData.Coh);
[rtcor rtcor_se] = calcGroupMean(trialData.RT, trialData.Coh, coh_set); % bin coherence levels, calculates mean and STD error for each

modelParam = struct('init', guess, 'fixed', fixed, 'hessian', [], 'final', [], 'se', []);
guess(3) = log(guess(3));
if all(fixed==1)
    modelLL = -fitmodel_MLEerr([]);
    disp('All parameters are fixed. No optimization performed!');
    return;
end

options = optimset('Display', 'final', 'MaxFunEvals', 1000*sum(fixed==0), 'MaxIter', 1000*sum(fixed==0)); % optimization function
[p, fval, exitflag] = fminsearch(@fitmodel_MLEerr, guess(fixed==0), options); % review fminsearch

modelParam.final = getParam(p);
modelParam.final(3) = exp(modelParam.final(3));
modelLL = -fval;

guess = modelParam.final;
fixed = ones(size(guess));
fitmodel_MLEerr([]);

%close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function err = fitmodel_MLEerr(p)
        
        param = getParam(p);
        
        K = param(1); % how units of evidence relate to increase/decreases in decision variable 
        B = param(2); % difference between decision thresholds/bounds (big bounds = conservative) 
        T0 = exp(param(3)); % time for non-decision processes
        a1 = param(4);
        a2 = param(5);
        
        
        %coh = ((1-exp(-a1*(trialData.Coh-a2)))./(1+exp(-a1*(trialData.Coh-a2)))+1)/2;
        coh = (1./(1+exp(-a1*(trialData.Coh-a2)))+1)/2;
        coh = (coh-min(coh))/range(coh);

        %y =  1./(1+exp(-10*x + 10));
%        cla(ax);
%        plot(ax, trialData.Coh, coh, '.');
        %drawnow;
        
        s = sqrt(2*0.5);

        trialData.expectedPC = 1./(1+exp(-2.*K.*abs(coh).*B./(s^2)));
        Icor = trialData.Resp==trialData.Stim & coh~=0;
        Iwrg = trialData.Resp~=trialData.Stim & coh~=0;
        
        
        %Here, the best models have the least negative log likelihood.  
           % log likelihood = how good model fits based on given
           % parameters; likelihood that the paramters fit
           
        %log of each expectedPC for correct retrieval trials.  When
        %high-coherence trials are correct, this value tends towards zero.
        %When low-coherence trials are correct, this value is greater in the negative direction.  
        
        %log of 1-expcted PC for incorrect retrieval trials.  When
        %high-coherence trials are incorrect, this value is greater in the negative direction.  When
        %low-coherence trials are incorrect, this value tends towards zero 
        LL_PC = sum(log(trialData.expectedPC(Icor))) + sum(log(max(1-trialData.expectedPC(Iwrg),0.001)));

        Ic0 = coh==0;
        trialData.expectedRT(Ic0) = T0 + B^2/s^2;
        trialData.expectedRT(~Ic0) = T0 + (B./(K*coh(~Ic0))).*tanh(K.*coh(~Ic0).*B./(s^2));

        if ~isequal(fitwhat,'PC')
            [~,g] = ismember(trialData.Coh, coh_set);
            
            %when expected RT is close to mean RT within relevant bin, this
            %value is less negative.  when expected RT is farther from mean
            %RT within relevant bin, this value is more negative.  
            LL_RT = sum(lognormpdf(trialData.expectedRT, rtcor(g), rtcor_se(g))); % (expectedRT vs. expectedRT') transposition
        end

        switch fitwhat
            case 'PC'
                err = -LL_PC;
            case 'RT'
                err = -LL_RT;
            case 'BOTH'
                err = -(LL_PC +LL_RT);
        end
    end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %this function retrieves the full parameter set given the adjustable and
        %fixed parameters
    function param = getParam(p)
        param(fixed==0) = p;                %get adjustable parameters from param1
        param(fixed==1) = guess(fixed==1);  %get fixed parameters from guess
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %returns log of normpdf, this function helps avoinding round off errors for very small probabilities.
        %if you use log(normpdf()) instead of lognormpdf you can get 0 for small probabilities and it's 
        %detrimental to log-likelihood fitting
    function l = lognormpdf(x, mu, sigma)
        d = sqrt(2*pi);
        l = -(x-mu).^2./(2*sigma.^2)+log(1./(d*sigma));
    end

end

    
    




