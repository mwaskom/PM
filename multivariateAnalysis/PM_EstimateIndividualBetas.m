function [] = PM_EstimateIndividualBetas(par)

orig_par = par;
if strcmp(par.subTask, 'loc')
    [~, idx] = fMRIBehAnalysis_Loc(par);
elseif strcmp(par.subTask, 'DM')
    [~, ~, idx] = Mnemonic_fMRIBehAnalysis_Retrieval(par);
end

fullOnsetsStruct = load(fullfile(par.analysisdir, 'ons.mat'));
fullOnsets = fullOnsetsStruct.onsets;

fullDesign{1} = zeros(size(par.numvols,1), 2);

clear data


if ~exist(par.denoisingBetaDir)
    mkdir(par.denoisingBetaDir)
end

if strcmp(par.betaEstimationAlg, 'glmEst')
    
    
% define noise regressors
for j=1:length(par.numvols)
    thisSess = par.scans_to_include(j);
    idxScansThisSess{j} = 1+sum(par.numvols(1:(j-1))):sum(par.numvols(1:j));
    
    theseScans{j} = par.rascanfiles.(par.subTask)(idxScansThisSess{j},:);
    thisDesign{j} = zeros(par.numvols(j), length(par.denoisingConds));
    
    sessBaselineScans = sum(par.numvols(1:(j-1)));
    allTrialsThisSess{j} = idx.alltrials(idx.sess==j);
    onsThisSess{j} = allTrialsThisSess{j};
    
    r = 0;
    for k=1:length(par.denoisingConds)
        idxTheseOnsets = (ismember(fullOnsetsStruct.names, par.denoisingConds{k}));
        theseOnsets = intersect(onsThisSess{j}, [fullOnsets{idxTheseOnsets}]);
        theseOnsetsInTRs = 1 + round((theseOnsets)/par.TR) - sessBaselineScans;
        if ~isempty(theseOnsetsInTRs)
            r = r+1;
            thisDesign{j}(theseOnsetsInTRs,r) = 1;
        end
    end
    
    clear v
    for s=1:length(theseScans{j})
        v_h = spm_vol(theseScans{j}(s,:));
        v{s} = spm_read_vols(v_h);
        %v_h = load_nii(theseScans{j}(s,:));
        %v{s} = single(v_h.img);
    end
    
    data{j} = single(cat(4,v{:}));
    
    
end

opt = par.denoiseOpt;
%[results,denoiseddata] = GLMdenoisedata(thisDesign,data,par.denoiseDur,par.TR,[],[],opt,[]);




    % estimate and write out betas
    for i=1:length(idx.alltrials)
        
        thisSess = idx.sess(i);
        
        thisSessString = ['run' prepend(num2str(par.scans_to_include(thisSess)))];
        thisSessDir = fullfile(par.denoisingBetaDir, thisSessString);
        opt.numpcstotry = 0;
        
        thisTrialDesign =  zeros(size(thisDesign{thisSess},1), 2);
        
        sessBaselineScans = sum(par.numvols(1:(thisSess-1)));
        
        idxThisTrial = 1 + round(idx.alltrials(i)/par.TR) - sessBaselineScans;
        
        idxNotThisTrial_h = setdiff(find(idx.sess==thisSess),i);
        idxNotThisTrial = 1 + round(idx.alltrials(idxNotThisTrial_h)/par.TR) - sessBaselineScans;
        
        thisTrialDesign(idxThisTrial,1) = 1;
        thisTrialDesign(idxNotThisTrial,2) = 1;
        
        opt.wantpercentbold = 0; %this is crucial! having patterns of % bold messes up classification.
        results2 = GLMestimatemodel({thisTrialDesign},data(thisSess),par.denoiseDur,par.TR,'assume',[],0,opt,[]);
        
        if ~exist(thisSessDir)
            mkdir(thisSessDir)
        end
        
        thisImgMat = (squeeze(results2.modelmd{2}(:,:,:,1)));
        thisVol = v_h;
        thisVol.fname = fullfile(thisSessDir, [par.denoisingBetaPrefix '_' prepend(num2str(i),3) '.nii']);
        thisVol.dt(1) = 16; %make single, instead of int
        
        spm_write_vol(thisVol, thisImgMat);
        %this_nii = v_h;
        %this_nii.fileprefix = fullfile(thisSessDir, [par.denoisingBetaPrefix '_' prepend(num2str(i),3) '.nii']);
        %this_nii.img = int16(squeeze(results2.modelmd{2}(:,:,:,1)));
        %save_nii(this_nii, this_nii.fileprefix)
        
    end
    
elseif strcmp(par.betaEstimationAlg, 'spm')
   
    
    for i=1:length(idx.alltrials)
        new_par = par;
        
        j = idx.sess(i);
        
        thisSess = par.scans_to_include(j);
           
        idxScansThisSess{thisSess} = 1+sum(par.numvols(1:(j-1))):sum(par.numvols(1:j));
        theseScans = orig_par.rascanfilesByRun{thisSess};
         
        sessBaselineScans = sum(par.numvols(1:(j-1)));
        
        allTrialsThisSess = idx.alltrials(idx.sess==j);
        allTrialsThisSesBaselined = allTrialsThisSess - par.TR*sessBaselineScans;
        
        
        onsets{1} = idx.alltrials(i) - par.TR*sessBaselineScans;
        onsets{2} = setdiff(allTrialsThisSesBaselined, onsets{1});
        
        names = {'trial1' 'allOtherTrials'};
        
        durations = {0 0};
        
        if ~exist(par.betaAnalysisDir)
            mkdir(par.betaAnalysisDir)
        end
        
        cd (par.betaAnalysisDir)
        R = [];
        
        save ons onsets names durations
        save regs.mat R
        
        new_par.analysisdir = par.betaAnalysisDir;
        new_par.rascanfiles.(par.subTask) = theseScans;
        new_par.numvols = par.numvols(j);
        new_par.numscans = sum(new_par.numvols);
        
        
        PM_mod_spec(new_par);
        PM_mod_est(new_par);
        
        thisSessString = ['run' prepend(num2str(thisSess))];
        thisSessDir = fullfile(par.denoisingBetaDir, thisSessString);
        
        if ~exist(thisSessDir)
            mkdir(thisSessDir)
        end
        
        betaNameImg = fullfile(thisSessDir, [par.denoisingBetaPrefix '_' prepend(num2str(i),3) '.img']);
        betaNameHdr = fullfile(thisSessDir, [par.denoisingBetaPrefix '_' prepend(num2str(i),3) '.hdr']);
         
        movefile('beta_0001.img', betaNameImg)
        movefile('beta_0001.hdr', betaNameHdr)
        
        rmdir(par.betaAnalysisDir,'s')
    end
    
    
end

