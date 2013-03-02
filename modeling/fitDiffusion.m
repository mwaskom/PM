
function [modelParam, modelLL, trialData] = fitDiffusion(trialData, guess, fixed, fitwhat)
% 
% [modelParam, modelLL, trialData] = fitDiffusion(trialData, guess, fitwhat)
% 

% 
% RK, 7/8/2010
% RK, modified on 9/20/2010 to include C0 in the parameters list 
% 

if notDefined('guess')
%     guess = [0.5 20 0.1 300];
    guess = [0.5 20 0.5 300];
end

if notDefined('fixed')
    fixed = zeros(size(guess));
end

if notDefined('fitwhat')
    fitwhat = 'BOTH';
elseif strcmp(fitwhat,'PC')
    fixed(end) = 1;
end


coh_set = unique(trialData.Coh);
[rtcor rtcor_se] = calcGroupMean(trialData.RT, trialData.Coh, coh_set);

modelParam = struct('init', guess, 'fixed', fixed, 'hessian', [], 'final', [], 'se', []);
if all(fixed==1)
    modelLL = -fitmodel_MLEerr([]);
    disp('All parameters are fixed. No optimization performed!');
    return;
end
    

options = optimset('Display', 'final', 'MaxFunEvals', 500*sum(fixed==0), 'MaxIter', 500*sum(fixed==0));
[p, fval] = fminsearch(@fitmodel_MLEerr, guess(fixed==0), options); 

modelParam.final = getParam(p);
modelLL = -fval;

guess = modelParam.final;
fixed = ones(size(guess));
fitmodel_MLEerr([]);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    function err = fitmodel_MLEerr(p)

        param = getParam(p);

        K = param(1);
        B = param(2);
        C0 = param(3);
        T0 = param(4);
        Coh = trialData.Coh+C0;

        s = sqrt(2*0.5);

        trialData.expectedPC = 1./(1+exp(-2.*K.*Coh.*B./(s^2)));
        Icor = trialData.Resp==trialData.Stim & Coh~=0;
        Iwrg = trialData.Resp~=trialData.Stim & Coh~=0;
        LL_PC = sum(log(trialData.expectedPC(Icor))) + sum(log(max(1-trialData.expectedPC(Iwrg),0.001)));

        Ic0 = Coh==0;
        trialData.expectedRT(Ic0) = T0 + B^2/s^2;
        trialData.expectedRT(~Ic0) = T0 + (B./(K*Coh(~Ic0))).*tanh(K.*Coh(~Ic0).*B./(s^2));

        if ~isequal(fitwhat,'PC')
            [~,g] = ismember(trialData.Coh, coh_set);
            LL_RT = sum(lognormpdf(trialData.expectedRT', rtcor(g), rtcor_se(g)));
        end

        switch fitwhat
            case 'PC'
                err = -LL_PC;
            case 'RT'
                err = -LL_RT;
            case 'BOTH'
                err = -(LL_PC+LL_RT);
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

    
    




