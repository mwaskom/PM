function [] = PM_PatLinRegGroupModel( subjArray, S )

%% gpar definition
gpar.subjArray = subjArray;



gpar.expt_dir = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/';
gpar.modelTemplate = '/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/scripts/GroupTemplate.mat';
gpar.constat = 'T';

gpar.exMask = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/groupMask/inclusive_mask.img';

for s = 1:length(S.patRegSet)
    thisPatReg = S.patRegSet{s};
    thisPatRegDir = [S.expt_dir S.subj_id '/' thisPatReg];
    
    gpar.conGroup = thisPatReg;
    gpar.tasks = {[thisPatReg '_13Subs']};
    
    for t = 1:length(gpar.tasks)
        
        d = dir([thisPatRegDir '/*.img']);
        
        for c= 1:length(d);
            
            gpar.task{t}.cons{c}.name = d(c).name(1:end-4);
            gpar.task{t}.cons{c}.dir = {fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t}, gpar.task{t}.cons{c}.name)};
            
            for s = 1:length(gpar.subjArray)
                gpar.task{t}.cons{c}.scans{s} =  fullfile(S.expt_dir, gpar.subjArray{s}, thisPatReg, d(c).name);
            end
            
        end
    end
    
    %% group model spec
    
    for t = 1:length(gpar.tasks)
        clear jobs
        load (gpar.modelTemplate)   %%/Users/alangordon/Accumulator_fMRI/AG1/scripts/GroupTemplate.mat;
        
        
        
        for c= 1:length(gpar.task{t}.cons);
            jobs{c}.stats{1}.factorial_design.des.t1.scans = [];
            
            if ~exist(fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t}, gpar.task{t}.cons{c}.name));
                mkdir(fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t}, gpar.task{t}.cons{c}.name));
            end
            jobs{c}.stats{1}.factorial_design.dir = {fullfile(gpar.expt_dir, 'group_analyses', gpar.tasks{t}, gpar.task{t}.cons{c}.name)};
            jobs{c}.stats{1}.factorial_design.masking.em{1} = gpar.exMask;
            
            jobs{c}.stats{1}.factorial_design.des.t1.scans = gpar.task{t}.cons{c}.scans;
            
        end
        
        spm_jobman('run',jobs)
        
    end
    
    %% group model estimate and set contrasts
    
    for tsk = 1:length(gpar.tasks)
        for cnd = 1:length(gpar.task{tsk}.cons)
            AG1_group_mod_est(gpar, tsk, cnd);
        end
    end
    
    for tsk = 1:length(gpar.tasks)
        for cnd = 1:length(gpar.task{tsk}.cons)
            AG1_groupsetcontrasts(gpar, tsk, cnd);
        end
    end
    
    
end
end

