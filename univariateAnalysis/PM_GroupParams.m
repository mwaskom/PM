function gpar = PM_GroupParams(subjArray, leftOutPar)

if nargin<2
   leftOutPar = []; 
end

[~, SA] = PM_SA;
gpar.subjArray = subjArray;

gpar.conGroup = { 'analysis_mnem_parmodByConfAndRT' };

if ~isempty(leftOutPar)
    gpar.tasks = {'analysis_mnem_parmodByConfAndRT_16Subs_leftOut'};
else
    gpar.tasks = {'analysis_mnem_parmodByConfAndRT_16Subs'};
end

gpar.expt_dir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/';
gpar.modelTemplate = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/scripts/analysis/univariateAnalysis/GroupTemplate.mat';
gpar.constat = 'T';
gpar.exMask = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/groupMask/inclusive_mask.img';
gpar.stat = 't1';
gpar.behTask = 'mnem';

for t = 1:length(gpar.tasks)
    
    %gpar.task{t}.conTemplate =     fullfile(gpar.expt_dir, 'pm_052311_2', gpar.conGroup{t}, 'SPM');
    gpar.task{t}.conTemplate =  fullfile(gpar.expt_dir, subjArray{1}, gpar.conGroup{t}, 'SPM');
    
    ldTemp = load(gpar.task{t}.conTemplate);
    
    gpar.task{t}.SPMcons = ldTemp.SPM.xCon;
    
    gpar.nCovs = 0;
    
    %erp_cov = load('/Users/alangordon/mounts/w5/alan/eegfmri/fmri_data/erp_data/PMMeanZscore_0p4to0p8sec_RHCH_LCH_LCCR_HCCR.mat');
 
    %for i=1:length(subjArray), subIDs(i) = EF_num2Sub(subjArray{i}); end
    %idxSubsToInclude = ismember(erp_cov.Y(:,1), subIDs);
    
    %gpar.covVec = erp_cov.Y(idxSubsToInclude,3) - erp_cov.Y(idxSubsToInclude,2);
    %gpar.covName = 'subjectWise_LCHits_vs_HCHitsAndRHits_MedParietal_400To800';
    
    %covVec = [0.8993 0.7149 0.7385 0.7512 0.6734 0.8600 0.8104 0.7800 ...
    %    0.8241 0.7067 0.8400 0.8387 0.9567 0.7179 0.8067 0.7400];
    %covVec = [0.8919    0.6452    0.7703    0.6757   0.6892    0.8267    0.8000    0.8000 ...
    %    0.8082    0.6400    0.8933    0.8108    0.9322    0.6892    0.7867    0.7333];
    %covVec = [0.9067    0.7846    0.7067    0.8267    0.6575    0.8933    0.8209    0.7600 ...
    %    0.8400    0.7733    0.7867    0.8667    0.9811    0.7467    0.8267    0.7467];

    
    %covName = 'pctCor';
    gpar.covVec = [];
    gpar.covName = [];
     
    for c= 1:length(gpar.task{t}.SPMcons);
        gpar.exMask = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/groupMask/inclusive_mask.img';
        
        if ~isempty(leftOutPar)
            gpar.task{t}.cons{c}.dir = {fullfile(leftOutPar.subdir, gpar.tasks{t}, gpar.task{t}.SPMcons(c).name)};
        else
            gpar.task{t}.cons{c}.dir = {fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t},gpar.task{t}.SPMcons(c).name)};
        end
        gpar.task{t}.cons{c}.name = gpar.task{t}.SPMcons(c).name;
        
        if strcmp(gpar.behTask,'mnem');
            if (~isempty(strfind(gpar.task{t}.cons{c}.name, 'Conf')) || ~isempty(strfind(gpar.task{t}.cons{c}.name, 'conf')))
                if ~isempty(strfind(gpar.task{t}.cons{c}.name, 'face_conf1'))
                    thisSubSet = SA.(gpar.behTask).sa16_mnemFaceConf1;
                else
                    thisSubSet = SA.(gpar.behTask).sa16_Conf;
                end
            elseif (~isempty(strfind(gpar.task{t}.cons{c}.name, 'Inc')) || ~isempty(strfind(gpar.task{t}.cons{c}.name, 'acc')))
                thisSubSet = SA.(gpar.behTask).sa16_CorVsInc;
            else
                thisSubSet = gpar.subjArray;
            end
        else
            thisSubSet = gpar.subjArray;
        end
        
        idxThisSubSet = ismember(gpar.subjArray,thisSubSet);
        subsToInclude = gpar.subjArray(idxThisSubSet);
        %gpar.covVec{c} = covVec(idxThisSubSet);
        %gpar.covName{c} = covName;
        nSubs = length(subsToInclude);
        for s = 1:length(subsToInclude)
            gpar.task{t}.cons{c}.scans{s} = fullfile(gpar.expt_dir, subsToInclude{s}, gpar.conGroup{t}, ['con_' prepend(num2str(c), 4) '.img']);
        end
        gpar.task{t}.cons{c}.groupContrast = {[1] [-1]};

    end
end