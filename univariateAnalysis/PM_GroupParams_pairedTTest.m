function gpar = PM_GroupParams_pairedTTest(SA1,SA2)

[~, SA] = PM_SA;
gpar.subjArray1 = SA1;
gpar.subjArray2 = SA2;

gpar.conGroup1 = {'analysisPercByAllConf'};
gpar.conGroup2 = {'analysisRetByAllConf'};
gpar.tasks = {'analysis_percVsMnem_Conf'};
gpar.expt_dir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/';
gpar.modelTemplate = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/templates/pairedTTemplate.mat';

gpar.constat = 'T';
gpar.exMask = [];
gpar.stat = 'pt';
for t = 1:length(gpar.tasks)
    
    %gpar.task{t}.conTemplate =     fullfile(gpar.expt_dir, 'pm_052311_2', gpar.conGroup{t}, 'SPM');
    gpar.task{t}.conTemplate1 =  fullfile(gpar.expt_dir, SA1{end-1}, gpar.conGroup1{t}, 'SPM');
    gpar.task{t}.conTemplate2 =  fullfile(gpar.expt_dir, SA2{end-1}, gpar.conGroup2{t}, 'SPM');
    
    ldTemp1 = load(gpar.task{t}.conTemplate1);
    ldTemp2 = load(gpar.task{t}.conTemplate2);
    
    gpar.task{t}.SPMcons1 = ldTemp1.SPM.xCon;
    gpar.task{t}.SPMcons2 = ldTemp2.SPM.xCon;
    
    gpar.nCovs = 0;
    
    gpar.covVec = [];
    gpar.covName = [];
     
    for c= 1:length(gpar.task{t}.SPMcons1);
        gpar.exMask = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/groupMask/inclusive_mask.img';
        
        gpar.task{t}.cons{c}.dir(1) = {fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t},gpar.task{t}.SPMcons1(c).name)};
        gpar.task{t}.cons{c}.dir(2) = {fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t},gpar.task{t}.SPMcons2(c).name)};

        gpar.task{t}.cons{c}.name = gpar.task{t}.SPMcons1(c).name;
        
        if (~isempty(strfind(gpar.task{t}.cons{c}.name, 'Conf')) || ~isempty(strfind(gpar.task{t}.cons{c}.name, 'conf')))
            thisSubSet1 = SA.perc.sa16_Conf;
            thisSubSet2 = SA.mnem.sa16_Conf;
        elseif (~isempty(strfind(gpar.task{t}.cons{c}.name, 'Inc')) || ~isempty(strfind(gpar.task{t}.cons{c}.name, 'acc')))
            thisSubSet1 = SA.perc.sa16_CorVsInc;
            thisSubSet2 = SA.mnem.sa16_CorVsInc;
        else
            thisSubSet1 = gpar.subjArray1;
            thisSubSet2 = gpar.subjArray2;
        end
        
        idxThisSubSet = ismember(gpar.subjArray1,thisSubSet1);
        %gpar.covVec{c} = covVec(idxThisSubSet);
        %gpar.covName{c} = covName;
        nSubs = length(thisSubSet1);
        for s = 1:length(thisSubSet1)
            gpar.task{t}.cons{c}.scans{s,1} = fullfile(gpar.expt_dir, thisSubSet1{s}, gpar.conGroup1{t}, ['con_' prepend(num2str(c), 4) '.img']);
            gpar.task{t}.cons{c}.scans{s,2} = fullfile(gpar.expt_dir, thisSubSet2{s}, gpar.conGroup2{t}, ['con_' prepend(num2str(c), 4) '.img']);
        end
        gpar.task{t}.cons{c}.groupContrast = {[zeros(1,nSubs) 1 -1] [zeros(1,nSubs) -1 1]};
    end
end