function [S idxTr idxTe par]= PM_mvpa_params(subj_id, task)

% establish parameters for mvpa analysis
% <subj_id> - identifier for the given subject. Can be numerical or a
% string
% <task> 'perc' or 'mnem'

%% establish general parameters
idxTr = [];
idxTe = [];
par = PM_Params(subj_id, task);
S.exp_name = 'PM';

%% directories
S.subj_id = par.substr;
S.expt_dir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/';
S.mvpa_dir = [S.expt_dir  S.subj_id '/analysis_mvpa_Loc_noPerf'];
S.anat_dir = [S.expt_dir  S.subj_id '/anat'];
S.importance_maps_dir=[S.expt_dir 'mvpa_results/ImpMaps_' date  ];
S.group_mvpa_dir = [S.expt_dir 'mvpa_files'];

if strcmp(task, 'perc')
    S.univar_dir = [par.subdir '/' 'analysis_loc_perc_3d'];
elseif strcmp(task, 'mnem')
    S.univar_dir = [par.subdir '/' 'analysis_loc_mnem'];
end
S.workspace_dir = [par.subdir '/' 'mvpa_workspace'];

%% preprocessing
S.preprocType = 'knk'; % 'spm' for spm preprocessing, 'knk' for kendrick preprocessing

%% tasks
S.trainTask = 'percloc';
S.testTask = 'perc';

%% cross-validation scheme
if strcmp(S.trainTask, S.testTask)
    S.xval = 1;
    S.thisSelector =  'randomNFold_xval'; % cross validation
else
    S.xval = 0;
    S.thisSelector = 'TrainTestOneIterGroup'; % train on one group, test on another
end

%% information specific to the training and testing of specific sets of data

% parTr = par for the training data
% S.onsetsTrainDir = location of training onsets
% S.condsTrain = conditions on which to train
% S.dnCondsTrain = conditions which which to denoise, if denoising is used
% S.TrainRuns = runs of data on which to train
% S.durTrain = duration of training
% S.filenames_train = names of data images to use for training
% idxTr = behavioral information for training task
    
if strcmp(S.trainTask,'mnemonicloc')
    parTr = PM_Params(subj_id, 'mnem');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_noPerf/'];
    S.condsTrain = {{'face'}  {'house'}} ;
    S.dnCondsTrain = {{'face'}  {'house'} {'noise'}};
    S.TrainRuns = par.scansSelect.(par.task).loc;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});
    %S.filenames_train = vertcat(par.betasByRun{S.TrainRuns});
    [~, idxTr] = fMRIBehAnalysis_Loc(par);
elseif strcmp(S.trainTask,'percloc_2sess')
    par2 = PM_Params3(subj_id, 'mnem');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTrain = {{'noise'}  {'houseCor'}} ;
    S.TrainRuns1 = par.scansSelect.(par.task).loc;
    S.TrainRuns2 = par2.scansSelect.(par2.task).loc;
    S.TrainRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.perc.all)];
    S.durTrain = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;
    S.filenames_train_h{1} = char(par.rascanfilesByRun{S.TrainRuns1});
    S.filenames_train_h{2} = char(par2.rascanfilesByRun{S.TrainRuns2});
    S.filenames_train = char(S.filenames_train_h);
elseif strcmp(S.trainTask,'mnemonicloc_2sess')
    parTr = PM_Params(subj_id, 'perc');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTrain = {{'faceCor'} {'houseCor' 'noise'} } ;
    S.TrainRuns1 = par.scansSelect.(par.task).loc;
    S.TrainRuns2 = parTr.scansSelect.(parTr.task).loc;
    S.TrainRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.mnem.all)];
    S.durTrain = sum([par.(par.task).numvols(S.TrainRuns1) parTr.(parTr.task).numvols(S.TrainRuns2)])* par.TR;
    S.filenames_train_h{1} = char(par.swrascanfilesByRun{S.TrainRuns1});
    S.filenames_train_h{2} = char(parTr.swrascanfilesByRun{S.TrainRuns2});
    S.filenames_train = char(S.filenames_train_h);
elseif strcmp(S.trainTask, 'ret')
    parTr = PM_Params(subj_id, 'mnem');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_mnemDM_3d/'];
    S.condsTrain = {{'face_cor_high' 'face_cor_low'} {'house_cor_high' 'house_cor_low'}};
    S.dnCondsTrain = {{'face'}  {'house'} {'noise'}};
    S.TrainRuns = par.scansSelect.(par.task).DM;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});
    %S.filenames_train = vertcat(par.betasByRun{S.TrainRuns});
    [~, ~, idxTr] = Mnemonic_fMRIBehAnalysis_Retrieval(parTr);
elseif strcmp(S.trainTask, 'retIntermixed')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/mvpa_ret_intermixed/'];
    S.condsTrain = {'face_cor' 'house_cor'};
    S.TrainRuns = [4 5 6 7 ];
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.rascanfilesByRun{S.TrainRuns});
elseif strcmp(S.trainTask,'percloc')
    parTr = PM_Params(subj_id, 'perc');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d/'];
    S.condsTrain = {{'faceCor'}  {'houseCor'}} ;
    S.TrainRuns = par.scansSelect.(par.task).loc;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});
    [~, idxTr] = fMRIBehAnalysis_Loc(par);
elseif strcmp(S.trainTask, 'perc')
    parTr = PM_Params(subj_id, 'perc');
    S.onsetsTrainDir = [parTr.subdir '/analysis_mvpa_percDMByCoh_3d/'];
    S.condsTrain = {{'house_coh100_cor' 'house_coh60_cor' 'house_coh45_cor' 'house_coh35_cor' 'coh0_resp_house'} {'null'}} ;
    S.TrainRuns = parTr.scansSelect.(parTr.task).DM;
    S.durTrain = sum(parTr.(parTr.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(parTr.swrascanfilesByRun{S.TrainRuns});
    [~, ~, idxTr] = Perceptual_fMRIBehAnalysis(parTr);
end

%% testing
if strcmp(S.testTask,'mnemonicloc')
    S.onsetsTestDir =[S.expt_dir S.subj_id '/analysis_mvpa_Loc_noPerf/'];
    S.condsTest = {{'face'}  {'house'}};
    S.TestRuns = par.scansSelect.(par.task).loc;
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    S.filenames_test = vertcat(par.betasByRun{S.TestRuns});
    [~, idxTe] = fMRIBehAnalysis_Loc(par);
elseif strcmp(S.testTask,'mnemonicloc_2sess') 
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTest = {{'faceCor'}  {'houseCor'}} ;
    S.TestRuns1 = par.scansSelect.(par.task).loc;
    S.TestRuns2 = par2.scansSelect.(par2.task).loc;
    S.TestRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.mnem.all)];
    S.durTest = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;  
elseif strcmp(S.testTask, 'percloc_2sess')
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTest = {{'faceCor'}  {'houseCor'}} ;
    S.TestRuns1 = par.scansSelect.(par.task).loc;
    S.TestRuns2 = par2.scansSelect.(par2.task).loc;
    S.TestRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.mnem.all)];
    S.durTest = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;
    S.filenames_test_h{1} = char(par.swrascanfilesByRun{S.TestRuns1});
    S.filenames_test_h{2} = char(par2.swrascanfilesByRun{S.TestRuns2});
    S.filenames_test = char(S.filenames_test_h);
elseif strcmp(S.testTask, 'ret')
    S.onsetsTestDir = [S.expt_dir par.substr '/analysis_mvpa_mnemDM_3d/'];
    S.condsTest = { {'face_cor' 'face_inc'} {'house_cor' 'house_inc'}};
    S.dnCondsTest = { {'face_cor'} {'house_cor'} {'face_inc'}  {'house_inc'}};
    S.TestRuns = par.scansSelect.(par.task).DM;
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    S.filenames_test = vertcat(par.swrascanfilesByRun{S.TestRuns});
    %S.filenames_test = vertcat(par.betasByRun{S.TestRuns});
    [~, ~, idxTe] = Mnemonic_fMRIBehAnalysis_Retrieval(par);
elseif strcmp(S.testTask,'percloc')
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d/'];
    S.condsTest = {{'faceCor'}  {'houseCor'}} ;
    S.TestRuns = par.scansSelect.(par.task).loc;
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    S.filenames_test = vertcat(par.swrascanfilesByRun{S.TestRuns});
elseif strcmp(S.testTask, 'perc')
    S.onsetsTestDir = [S.expt_dir S.subj_id '/analysis_mvpa_percDMByCoh_3d'];
    S.condsTest = {{'face_cor' 'face_inc' 'coh0_resp_face'} {'house_cor' 'house_inc' 'coh0_resp_house'}};
    S.TestRuns = par.scansSelect.(par.task).DM;
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    S.filenames_test = vertcat(par.swrascanfilesByRun{S.TestRuns});
    [~, ~, idxTe] = Perceptual_fMRIBehAnalysis(par);
end

S.condnames = S.condsTrain;
S.regName = 'conds';


%% Smoothing Parameters
S.funcType = 2;
S.smoothTxt = { 'unsmoothed' 'smoothed' 'native'};
switch S.funcType
    case 1
    par.filesForPatterns = par.wrascanfiles.all;
    case 2
    par.filesForPatterns = par.swrascanfiles.all;
    case 3
    par.filesForPatterns = par.rascanfiles.all;     
end

%% specify which files to load for classification
if S.xval
    S.filenames = S.filenames_train;
else
    S.filenames_h{1} = S.filenames_train;
    S.filenames_h{2} = S.filenames_test;
    S.filenames = char(S.filenames_h);
end
S.img_files =  mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);

%% Runs Parameters

%S.runs_vector - number of volumes per each run
if ismember(S.trainTask,{'percloc_2sess' 'mnemonicloc_2sess'})  
    if S.xval
        S.runs_vector = [par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)];
    else
        S.runs_vector = [par.(par.task).numvols(S.TrainRuns1) parTr.(parTr.task).numvols(S.TrainRuns2) par.(par.task).numvols(S.TestRuns)];
    end
else
    if S.xval
        S.runs_vector =  [par.(par.task).numvols(S.TrainRuns)];
    else
        S.runs_vector =  [parTr.(parTr.task).numvols(S.TrainRuns) par.(par.task).numvols(S.TestRuns)];
        S.TrainTestOrder = [ones(size(S.TrainRuns)) 2*ones(size(S.TestRuns))];
    end
end

S.meta_runs = S.runs_vector;
S.num_runs = length(S.runs_vector);
S.num_vols = sum(S.runs_vector);
S.TR = par.TR;

%% Volume Parameters
S.vol_info = spm_vol(fullfile(par.funcdir, 'run01', 'swRrun010030.nii')); %get functional data resolution info for spm .img writing

S.roiWithNonTaksVoxels = fullfile(par.anatdir, 'tnativeOccTemp.nii');
S.roiWithNonTaksVoxelsName = 'tnativeOccTemp.nii';

%S.roi_name = ['parietalNoPostCentral.img'];
%S.roi_file = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/parietalNoPostCentral.img';
S.roi_file = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/OTnoHipp.img';
S.roi_name = 'OTnoHipp.img';

%S.roi_name = 'rc1V001.nii';
%S.roi_file = fullfile(S.anat_dir, S.roi_name);

%S.roi_name = 'tnativeOccTemp.nii';
%S.roi_file = fullfile(par.anatdir, S.roi_name);

S.noiseVoxels_file = fullfile(par.anatdir, 'rnativec2V001.nii');
S.noiseVoxels_name = 'rnativec2V001.nii';

S.sigVoxels_file = fullfile(par.anatdir, 'tnativeOccTempGrey.nii');
S.sigVoxels_name = 'tnativeOccTempGrey.nii';

S.secondaryMask = []; % secondary mask (the specific classification mask)
%masks the primary data loaded in the workspace. [] = no secondary mask.
%S.secondaryMask = fullfile(par.anatdir, 'tnativeOccTempGrey.nii');
%S.secondaryMask = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/OTnoHipp.img'];
%S.secondaryMask = ['/Users/gordonam/Studies/AG1/mvpa/OT_2012/rhipp.img'];
%S.secondaryMask = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/rparietal.img'];


%% Workspace Parameters
S.use_premade_workspace = 1;
S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.funcType} '_train_' S.trainTask '_test_' S.testTask S.preprocType '.mat']);

%% Pattern names
S.patternType = 'raw'; %'raw' or 'betas'
% S.preprocPatName = 'spiral_dn';
% S.preprocPatName = 'patsAllVox_z_dn';
S.preprocPatName = 'spiral_hp_z';
S.preprocPatCondensedName = [S.preprocPatName '_condensed'];

if isempty(S.secondaryMask)
    S.preprocPatNameFinalMask = S.preprocPatName;
else
    S.preprocPatNameFinalMask = [S.preprocPatName '_masked'];
end

%% Artifacts
S.artFile = (fullfile(par.artrepdir, ['art_global_modified_' par.substr])); %directory where artifact information is contained
S.inactivateArtifacts = 0; %remove artifact trials?

%% Iteration Parameters
S.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
S.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data

%% Balancing Parameters
S.equate_number_of_trials_in_groups = 0; % equate number of trials in conditions 1 and 2
S.numBalancedIts = 1; % number of iterations to run, with different randomization for the balancing

%% Z-Scoring and outlier detection
S.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification 
S.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
S.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers
S.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers
S.remove_outlier_trials = 0;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)

%% Importance Maps
S.generate_importance_maps = 0; %visualize classifier weights
S.generateBetaMaps = 1; %use betas, instead of importance values
S.impType = {'pos' 'neg' 'both' 'raw'}; %importance map types
S.regNames = {'face' 'house'};

%% Special types of analysis
S.searchlightAnalysis = 0; % run a searchlight analysis
S.linReg = 0; % run an analysis with a continuous outcome variable

%% Subsample
S.subsampleToMatch = 0; %subsample trials to match quantities across them.
S.balanceHiAndLowConf = 0;% match N of hi and low confidence trials?

%% voxel interactions
S.includeVoxelInteractions = 0; %include interactions among voxels?   
S.interactionType = 2;
S.intEst = 1;
S.intConcat = 0;
S.intPThresh = .001;
S.intReportIncrement = 100;
S.intFilePrefix = 'intVox';
S.intMaskName = 'interactionMask';
S.intPatName = 'interactions';
S.intGroupName = 'interactionsGroup';
S.intUseIntsWithUniqueInfo = 1;

%% Signal intensity analysis
S.thisSigIntenseSelector = 'randomNFold_xval'; %which selector to use for signal intensity analysis
S.zscoreIntensityVals = 1; % zscore the intensity values?

%% Denoising
S.denoise = 0; %undergo denoising?
S.denoiseOpt.denoisespec = '10001'; %which parts of the glm output do we want to save?

%% Mean Signal Extraction Params
% parameters for selecting the mean signal from a class-specific ROI for each pattern.

S.extractMeanSignal = 0; %1 - do signal extraction. 0 - don't do this. 
S.defineROIsFromANOVAFS = 0; % define ROIs using ANOVA-based feature selection, instead of pre-defining them. 
S.logreg_2Features = 0; %perform a logistic regression, using the two extracted intensity vectors

S.ROI1PatName = [S.preprocPatCondensedName '_ROI1'];
S.ROI1_name = [ 'occipitoTemporal_faceVsScene_500vox.img'];
S.ROI1_file  = [par.subdir '/analysis_loc_mnem/' S.ROI1_name];

S.ROI2PatName = [S.preprocPatCondensedName '_ROI2'];
S.ROI2_name  = ['occipitoTemporal_sceneVsFace_500vox.img'];
S.ROI2_file   = [par.subdir '/analysis_loc_mnem/' S.ROI2_name];

%% TR Weighting
%which post-stimulus TRs should be used (and if more than one, averaged
%across) before feeding data to the classifier?  S.TR_weights_set{1} -
%training weights, S.TR_weights_set{2} = testing weights
S.TR_weights_set = {[0 0 0 .5 .5] [0 0 0 .5 .5]}; 

%% classifier parameters
S.class_args.train_funct_name = 'train_liblinear'; %training function
S.class_args.test_funct_name = 'test_liblinear'; %testing function
S.class_args.classType = 'libLin';
S.perfmet_functs = 'perfmet_maxclass'; % performance metric
S.statmap_funct = 'AG_statmap_anova'; % performance metric
S.nPlsCompsSet = 0; % number of pls components to include. 0 = do not use pls components.
S.nFolds = 10; % number of cross validation iterations

S.class_args.nVox = 0; % number of voxels to select with feature selection e.g. [1000 5000 10000]
S.class_args.libLin = '-q -s 0 -B 1'; %arguments for liblinear
S.class_args.libsvm = '-q -s 0 -t 2 -d 3'; % arguments for libsvm
S.class_args.constant = true; % include a constant term?
S.class_args.prefitWeights = true; 
S.class_args.chooseOptimalPenalty = 1; % cycle through cost parameters in the training set, and chose the optimal one?
S.class_args.penaltyRange = [.001 .005 .01 .05 .1 .5 1 5 10 50 100 500 1000 50000]; % cost parameters to cycle through
S.class_args.radialBasisSelection = [];%[.00001 .0001 .001 .01 .1 1 10];
S.class_args.nFoldsPenaltySelection = 10; % number of cross validation folds for penalty parameter selection. 

end