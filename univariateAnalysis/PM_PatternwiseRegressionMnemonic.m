function S = PMPatternwiseRegression_Mnemonic(subidx, task)


%% Params
par = PM_Params(subidx, task);
subj_id = subidx;
S.subj_id = par.substr;
S.expt_dir = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/';
S.univar_dir = [S.expt_dir S.subj_id '/' 'analysis_loc_mnem'];
S.workspace_dir = [S.expt_dir S.subj_id '/' 'mvpa'];
S.analysisName = 'Ret';
S.exp_name = 'PM';


%% Condition Parameters
S.trainTask = 'ret';
S.testTask = 'ret';

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
    S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});
elseif strcmp(S.trainTask,'percloc_2sess')
    par2 = PM_Params(subj_id, 'mnem');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTrain = {'faceCor'  'houseCor'} ;
    S.TrainRuns1 = par.scansSelect.(par.task).loc;
    S.TrainRuns2 = par2.scansSelect.(par2.task).loc;
    S.TrainRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.perc.all)];
    S.durTrain = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;
    S.filenames_train_h{1} = char(par.swrascanfilesByRun{S.TrainRuns1});
    S.filenames_train_h{2} = char(par2.swrascanfilesByRun{S.TrainRuns2});
    S.filenames_train = char(S.filenames_train_h);
elseif strcmp(S.trainTask,'mnemonicloc_2sess')
    par2 = PM_Params(subj_id, 'perc');
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d_2sess/'];
    S.condsTrain = {'faceCor'  'houseCor' } ;
    S.TrainRuns1 = par.scansSelect.(par.task).loc;
    S.TrainRuns2 = par2.scansSelect.(par2.task).loc;
    S.TrainRuns = [S.TrainRuns1 S.TrainRuns2 + max(par.scansSelect.mnem.all)];
    S.durTrain = sum([par.(par.task).numvols(S.TrainRuns1) par2.(par2.task).numvols(S.TrainRuns2)])* par.TR;
    S.filenames_train_h{1} = char(par.swrascanfilesByRun{S.TrainRuns1});
    S.filenames_train_h{2} = char(par2.swrascanfilesByRun{S.TrainRuns2});
    S.filenames_train = char(S.filenames_train_h);
elseif strcmp(S.trainTask, 'ret')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_mnemDM_3d/'];
    S.condsTrain = {'face_cor' 'house_cor'};
    S.TrainRuns = par.scansSelect.(par.task).DM;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});
elseif strcmp(S.trainTask, 'retIntermixed')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/mvpa_ret_intermixed/'];
    S.condsTrain = {'face_cor' 'house_cor'};
    S.TrainRuns = [4 5 6 7 ];
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});
elseif strcmp(S.trainTask,'percloc')
    S.onsetsTrainDir = [S.expt_dir S.subj_id '/analysis_mvpa_Loc_3d/'];
    S.condsTrain = {'faceCor'  'houseCor'} ;
    S.TrainRuns = par.scansSelect.(par.task).loc;
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swrascanfilesByRun{S.TrainRuns});
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
    par.filesForPatterns = par.wrascanfiles.all;
else
    par.filesForPatterns = par.swrascanfiles.all;
end

if S.xvalTrainData
    S.filenames = S.filenames_train;
else
    S.filenames_test = vertcat(par.swrascanfilesByRun{S.TestRuns});
    S.filenames_h{1} = S.filenames_train;
    S.filenames_h{2} = S.filenames_test;
    S.filenames = char(S.filenames_h);
end

S.img_files =  mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);


%% Volume Params
S.roi_name = 'mask.img';
S.roi_file = [S.univar_dir '/' S.roi_name];
S.vol_info = spm_vol(fullfile(S.univar_dir, 'beta_0001.hdr'));

%% Workspace Parameters
S.use_premade_workspace = 1;
S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.use_unsmoothed + 1} '_train_' S.trainTask '_test_' S.testTask '.mat']);
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


%% classFile
classMatDescriptor = 'linRidge_14Subs_FaceOnly.mat';
classMatFile = ['/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/mvpa_files/' classMatDescriptor] ;

classMatDescriptor2 = 'linRidge_14Subs_HouseOnly.mat';
classMatFile2 = ['/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/mvpa_files/' classMatDescriptor2 ];

%% Workspace loading / creation.
existWorkspace = exist(S.workspace);

if (S.use_premade_workspace&&existWorkspace) % if we are supposed to use premade workspace, and one with the correct name exists
    load(S.workspace, 'subj', 'condensed_regs_of_interest', 'condensed_runs');
else
    [subj] = PM_mvpa_load_and_preprocess_raw_data(S);
    
    % load onsets
    [S.onsets S.RTs S.enoughTRs] = PM_mvpa_onsets_and_images(S);
    S.num_conds = size(S.onsets,2);
    
    % turn onsets into regressors matrix
    all_regs = zeros(S.num_conds,S.num_vols); % initialize regs matrix as conditions x timepoints
    for cond = 1: S.num_conds
        for trial = 1: length(S.onsets{cond})
            time_idx = round(S.onsets{cond}(trial)/S.TR) + 1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
            all_regs(cond, round(time_idx)) = 1;
        end
    end
    
    % condense regs by removing zeros
    condensed_runs = [];
    condensed_regs_of_interest = [];
    trial_counter = 1;
    for i = 1: size(all_regs,2)
        if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
            condensed_regs_of_interest(:,trial_counter) = all_regs(:,i);
            condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
            trial_counter = trial_counter + 1;
        end
    end
    
    % create meta runs structure
    idx_condense =find(sum(all_regs));
    all_trials = sum(all_regs,1);
    meta_runs = find(all_trials);
    
    % condense data
    data_by_TR = [];
    for dt = 1:length(S.TR_weights)
        data_by_TR(dt,:,:) = S.TR_weights(dt)*subj.patterns{end}.mat(:,meta_runs+(dt-1));
    end
    temporally_condensed_data = squeeze(sum(data_by_TR(S.TRs_to_average_over,:,:),1));
    
    % add new condensed activation pattern
    subj = duplicate_object(subj,'pattern',S.preprocPatName,S.preprocPatCondensedName);
    subj = set_mat(subj,'pattern',S.preprocPatCondensedName,temporally_condensed_data,'ignore_diff_size',true);
    subj = remove_mat(subj,'pattern', S.preprocPatName);
    subj = remove_mat(subj,'pattern','spiral');
    save (S.workspace)
end

%% post workspace loading
[res psy idx] = Mnemonic_fMRIBehAnalysis_Retrieval(par);

% load classmat
classMat = load(classMatFile);
classMat2 = load(classMatFile2);
subjArrayIdx = classMat.qqq.subjArray==subj_id;

% load relevant performance vectors
%corTrials = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations.perfmet.corrects;
%allActs = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations.acts(1,:);
allPats = get_mat(subj, 'pattern', S.preprocPatCondensedName);
regs = condensed_regs_of_interest(1,:);
allConf = idx.certRat(idx.cor==1);
%allClass = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations.perfmet.desireds;



for k = 1:length(classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations)
    corTrials_h{k} = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.corrects;
    allActs_h{k} = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).acts(1,:);
    allDesireds_h{k} = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.desireds;
    allClass_h{k} = classMat.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.desireds;

end

corTrials = [corTrials_h{:}];
allActs = [ allActs_h{:}];
allClass = [ allClass_h{:}];
allActsLogit = log(allActs ./ (1- allActs)); %use logit instead of probability
allRT = idx.rt;

for k = 1:length(classMat2.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations)
    corTrials_h2{k} = classMat2.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.corrects;
    allActs_h2{k} = classMat2.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).acts(1,:);
    allDesireds_h2{k} = classMat2.qqq.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.desireds;
end

% turn this from face evidence into general evidence in the correct
% direction, by flipping evidence for scenes.
%allActs(allClass==2) = 1-allActs(allClass==2);
%allActsLogit(allClass==2) = allActsLogit(allClass==2);

allActs2 = horzcat(allActs_h2{:});
%allActsLogit2 = log(allActs2 ./ (1- allActs2));



% session
allSess = zeros(max(condensed_runs), length(regs));
for i = 1:length(regs)
   allSess(condensed_runs(i),i) = 1; 
end

% turn the final session column into a constant column
allSess(end,:) = 1;

%idxClass{1} = (idx.high(idx.cor==1) + idx.low(idx.cor==1)) == 1;
% idxClass{1} = logical((idx.high(idx.cor==1) + idx.low(idx.cor==1)) .* (idx.face==1) );
% idxClass{2} = logical((idx.high(idx.cor==1) + idx.low(idx.cor==1)) .* (idx.house==1) );

idxClass{1} = logical(idx.face(idx.cor==1));
idxClass{2} = logical(idx.house(idx.cor==1));

% allActsDirectionOfCorrectSource = allActs;
% allActsDirectionOfCorrectSource(idxClass{2}) = -allActs2(idxClass{2});
% 
% allActsLogitDirectionOfCorrectSource = allActsLogit;
% allActsLogitDirectionOfCorrectSource(idxClass{2}) = -allActsLogit2(idxClass{2});


%className{1} = 'percCombinedFaceAndSceneBySignedEv';
 className{1} = 'percFace';
 className{2} = 'percHouse';

S.patRegSet = {'mnemByConf'};

betaMapsSet = { {'conf'} };

for bS = 1:length(betaMapsSet)
    S.patReg = S.patRegSet{bS};
    S.patRegDir = [S.expt_dir S.subj_id '/' S.patReg];
    for i = 1:length(idxClass)
        pat{i} = allPats(:,idxClass{i});
        sess{i} = allSess(:,idxClass{i});
        %acts{i} = allActsDirectionOfCorrectSource(idxClass{i});
        conf{i} = allConf(idxClass{i});
        %class{i} = allClass(idxClass{i});
        rt{i} = allRT(idxClass{i});
        
        datMat = {[conf{i}'] };
        
        X{i} = [datMat{bS} sess{i}'];
        b{i} = nan(size(X{i} ,2), size(pat{i},1));
    end
    
    
    for i = 1:length(pat)
        nRegRows = size(X{i},2) - size(sess{i},1);
        XZ = [zscore(X{i}(:,1:nRegRows)) sess{i}'];
        for j = 1:size(pat{i},1)
            
            thisPat = pat{i}(j,:)';
            thisPatZ = zscore(thisPat);
            
            [b{i}(:,j)] = regress(thisPatZ, XZ);
        end
    end
    
    vol_info = S.vol_info;
    
    for p = 1:length(subj.patterns)
        patNames{p} = subj.patterns{p}.name;
    end
    
    voxel_inds = subj.masks{end}.mat;
    
    betaMaps = betaMapsSet{bS};
    
    %%
    contrastName= {'acrossFacesAndHouses' 'FacesVsHouses'};
    imcalcText = {'(i1 + i2)/2' 'i1 - i2'};
    
    for j = 1:length(betaMaps)
        for i = 1:length(pat)
            
            outMap=  zeros(vol_info.dim);
            
            outMap(voxel_inds==1) = b{i}(j,:);
            
            vol_info.dir = S.patRegDir;
            vol_info.fname = [ vol_info.dir '/' className{i} '_' betaMaps{j} '.img'];
            
            S.pat{i}.map{j}.fname = vol_info.fname;
            
            if isempty(dir([vol_info.dir]))
                mkdir(vol_info.dir);
            end
            
            spm_write_vol(vol_info,outMap);
            
            fprintf('\n wrote out %s', vol_info.fname);
            
            thisScript = which('PM_PatternwiseRegressionMnemonic.m');
            copyfile(thisScript, vol_info.dir);
            
            d{i} = vol_info.fname;
        end
        
        dirChar = char(d);
        
        % print out contrasts
        for c = 1:length(contrastName)
            outputfileName = ['con_' betaMaps{j} '_' contrastName{c} '.img'];
            outputfile = [ vol_info.dir '/' outputfileName];
            spm_imcalc_ui(dirChar,outputfile, imcalcText{c})
        end
    end
end




