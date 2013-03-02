function [S par] = PMPatternParams( subidx, task )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Params
par = PM_Params(subidx, task);
subj_id = subidx;
S.sub_no = subj_id;
S.subj_id = par.substr;
S.expt_dir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/';
S.univar_dir = [S.expt_dir S.subj_id '/' 'analysis_loc_perc_3d'];
S.workspace_dir = [S.expt_dir S.subj_id '/' 'mvpa'];
S.analysisName = 'Perc_byEv';
S.exp_name = 'PM';

%% Analysis types
S.MeanROIAnalysis = 0;
S.WriteRegressionMaps = 1;

%% Condition Parameters
S.trainTask = 'perc';
S.testTask = 'perc';

if strcmp(S.trainTask, S.testTask)
    S.xvalTrainData = 1;
else
    S.xvalTrainData = 0;
end

if strcmp(S.trainTask,'mnemonicloc')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d/'];
    S.condsTrain = {'faceCor'  'houseCor'} ;
    S.TrainRuns = par.scansSelect.(par.task).loc;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swascanfilesByRun{S.TrainRuns});
elseif strcmp(S.trainTask,'percloc_2sess')
    par2 = PM_Params2(subj_id, 'mnem');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTrain = {'faceCor'  'houseCor'} ;
    S.TrainRuns1 = par.scansSelect.(par.task).loc;
    S.TrainRuns2 = par2.scansSelect.(par2.task).loc;
    S.TrainRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.perc.all)];
    S.durTrain = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;
    S.filenames_train_h{1} = char(par.swascanfilesByRun{S.TrainRuns1});
    S.filenames_train_h{2} = char(par2.swascanfilesByRun{S.TrainRuns2});
    S.filenames_train = char(S.filenames_train_h);
elseif strcmp(S.trainTask,'mnemonicloc_2sess')
    par2 = PM_Params(subj_id, 'perc');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTrain = {'faceCor'  'houseCor' } ;
    S.TrainRuns1 = par.scansSelect.(par.task).loc;
    S.TrainRuns2 = par2.scansSelect.(par2.task).loc;
    S.TrainRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.mnem.all)];
    S.durTrain = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;
    S.filenames_train_h{1} = char(par.swascanfilesByRun{S.TrainRuns1});
    S.filenames_train_h{2} = char(par2.swascanfilesByRun{S.TrainRuns2});
    S.filenames_train = char(S.filenames_train_h);
elseif strcmp(S.trainTask, 'ret')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_mnemDM_3d/'];
    S.condsTrain = {'face_cor' 'house_cor'};
    S.TrainRuns = par.scansSelect.(par.task).DM;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swascanfilesByRun{S.TrainRuns});
elseif strcmp(S.trainTask, 'retIntermixed')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/mvpa_ret_intermixed/'];
    S.condsTrain = {'face_cor' 'house_cor'};
    S.TrainRuns = [4 5 6 7 ];
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swascanfilesByRun{S.TrainRuns});
elseif strcmp(S.trainTask,'percloc')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d/'];
    S.condsTrain = {'faceCor'  'houseCor'} ;
    S.TrainRuns = par.scansSelect.(par.task).loc;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swascanfilesByRun{S.TrainRuns});
elseif strcmp(S.trainTask, 'perc')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_percDMByCoh_3d/'];
    S.condsTrain = {'resp_face' 'resp_house'};
    S.TrainRuns = par.scansSelect.(par.task).DM;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});
end

if strcmp(S.testTask,'mnemonicloc')
    S.onsetsTestDir =[S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d/'];
    S.condsTest = {'faceCor'  'houseCor'} ;
    S.TestRuns = par.scansSelect.(par.task).loc;
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
elseif strcmp(S.testTask,'mnemonicloc_2sess') 
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTest = {'faceCor'  'houseCor'} ;
    S.TestRuns1 = par.scansSelect.(par.task).loc;
    S.TestRuns2 = par2.scansSelect.(par2.task).loc;
    S.TestRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.mnem.all)];
    S.durTest = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;  
elseif strcmp(S.testTask, 'percloc_2sess')
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTest = {'faceCor'  'houseCor'} ;
    S.TestRuns1 = par.scansSelect.(par.task).loc;
    S.TestRuns2 = par2.scansSelect.(par2.task).loc;
    S.TestRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.mnem.all)];
    S.durTest = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;
elseif strcmp(S.testTask, 'ret')
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_mnemDM_3d/'];
    S.condsTest = {'face_cor' 'house_cor'};
    S.TestRuns = par.scansSelect.(par.task).DM;
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
elseif strcmp(S.testTask,'percloc')
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d/'];
    S.condsTest = {'faceCor'  'houseCor'} ;
    S.TestRuns = par.scansSelect.(par.task).loc;
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
elseif strcmp(S.testTask, 'perc')
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_percDMByCoh_3d'];
    S.condsTest = {'resp_face' 'resp_house'};
    S.TestRuns = par.scansSelect.(par.task).DM;
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
end

%% Smoothing Parameters
S.smoothTxt = {'smoothed' 'unsmoothed'};
S.use_unsmoothed = false;

if S.use_unsmoothed
    par.filesForPatterns = par.wascanfiles.all;
else
    par.filesForPatterns = par.swrascanfiles.all;
end

if S.xvalTrainData
    S.filenames = S.filenames_train;
else
    S.filenames_test = vertcat(par.swascanfilesByRun{S.TestRuns});
    S.filenames_h{1} = S.filenames_train;
    S.filenames_h{2} = S.filenames_test;
    S.filenames = char(S.filenames_h);
end

S.img_files =  mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);


%% Volume Params
S.roi_name = 'inclusive_mask.img';
S.roi_file = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/groupMask/' S.roi_name];

S.vol_info = spm_vol(fullfile(S.univar_dir, 'spmT_0001.img'));

%% Workspace Parameters
S.use_premade_workspace = 1;
S.preprocType = 'spm_preproc';
S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.use_unsmoothed + 1} '_train_' S.trainTask '_test_' S.testTask S.preprocType '.mat']);
S.preprocPatName = 'spiral_d_z';
S.preprocPatCondensedName = 'spiral_d_z_condensed';


%% Runs Parameters
if ismember(S.trainTask,{'percloc_2sess' 'mnemonicloc_2sess'})  
    if S.xvalTrainData
        S.runs_vector = [par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)];
    else
        S.runs_vector = [par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2) par.(par.task).numvols(S.TestRuns)];
    end
else
    if S.xvalTrainData
        S.runs_vector =  [par.(par.task).numvols(S.TrainRuns)];
    else
        S.runs_vector =  [par.(par.task).numvols(S.TrainRuns) par.(par.task).numvols(S.TestRuns)];
    end
end

S.meta_runs = S.runs_vector;
S.num_runs = length(S.runs_vector);
S.num_vols = sum(S.runs_vector);
S.TR = 2;

%% TRs
S.TR_weights = [0 0 .5 .5];
S.TRs_to_average_over = 1:length(S.TR_weights);

%% Artifacts
S.artStruct = load(fullfile(par.artrepdir, ['art_global_modified_' par.substr]));

%% classFile
S.classMatDescriptor = 'trainTestPerc_forImpmaps.mat';
S.classMatFile = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/mvpa_files/' S.classMatDescriptor ];

S.classMatDescriptor2 = 'trainTestPerc_forImpmaps.mat';
S.classMatFile2 = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/mvpa_files/' S.classMatDescriptor2 ];
end

