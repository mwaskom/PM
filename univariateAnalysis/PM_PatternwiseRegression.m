function [S psyphys] = PM_PatternwiseRegression(subidx, task)

[S par] = PM_PatternParams(subidx,task);

psyphys = [];

%% Workspace loading / creation.
existWorkspace = exist(S.workspace);

if (S.use_premade_workspace&&existWorkspace) % if we are supposed to use premade workspace, and one with the correct name exists
    load(S.workspace, 'subj', 'condensed_regs_of_interest', 'condensed_runs', 'all_regs');
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
[res psy idx] = Perceptual_fMRIBehAnalysis(par);

% load classmat
classMat = load(S.classMatFile);
classMat2 = load(S.classMatFile2);
subjArrayIdx = classMat.res.subjArray==S.sub_no ;

% load relevant performance vectors
for k = 1:length(classMat.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations)
    corTrials_h{k} = classMat.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.corrects;
    allActs_h{k} = classMat.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).acts(1,:);
    allNActs_h{k} = classMat.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).acts(1:2,:);
    allDesireds_h{k} = classMat.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.desireds;
end

for k = 1:length(classMat2.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations)
    corTrials_h2{k} = classMat2.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.corrects;
    allActs_h2{k} = classMat2.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).acts(1,:);
    allDesireds_h2{k} = classMat2.res.subj{subjArrayIdx==1}.penalty.nVox.weights.iter{1}.iterations(k).perfmet.desireds;
end

allActs = horzcat(allActs_h{:});
allActs2 = horzcat(allActs_h2{:});
corTrials = horzcat(corTrials_h{:});
allDesidreds = horzcat(allDesireds_h{:});
allClass = idx.face + 2*idx.house;
allRespClass = idx.resp_face + 2*idx.resp_house;

allActsWithoutNoise = allActs ./ sum(horzcat(allNActs_h{:}));
allNActs_h2 = [allNActs_h{:}];
faceActs = allNActs_h2(1,:);
houseActs = allNActs_h2(2,:);
%noiseActs = horzcat(allNActs_h{:}(3,:));


allRT = idx.rt;

allPats = get_mat(subj, 'pattern', S.preprocPatCondensedName);
regs = condensed_regs_of_interest(1,:);
allCoh = idx.coh_;
allConf = idx.conf_unsigned; 


allActsLogit = log(allActs ./ (1- allActs)); %use logit instead of probability
allActsLogit2 = log(allActs2 ./ (1- allActs2));
%allActsLogitFace = log(faceActs ./ (1- faceActs));
%allActsLogitHouse = log(houseActs ./ (1- houseActs));
%allActsLogitNoise = log(noiseActs ./ (1- noiseActs));
%allActsWithoutNoiseLogit = log(allActsWithoutNoise ./ (1- allActsWithoutNoise)); 

% turn this from face evidence into general evidence in the correct
% direction, by flipping evidence for scenes.
%allActs(allClass==2) = -allActs(allClass==2);



%allActsDiff = (allActs + allActs2); 
%allActsDiff(allRespClass==2) = -allActsDiff(allRespClass==2);

% session
allSess = zeros(max(condensed_runs), length(regs));
for i = 1:length(regs)
   allSess(condensed_runs(i),i) = 1; 
end

% turn the final session column into a constant column
allSess(end,:) = 1;

%% artifact index

rawSigInts = S.artStruct.zscoreA_cell(1:length(S.TestRuns));
rawMotDerivs = S.artStruct.delta_cell(1:length(S.TestRuns));

catSigInts = vertcat(rawSigInts{:});
catMotDerivs = vertcat(rawMotDerivs{:});

Art = ((abs(catSigInts) > par.art.sigThresh) + (abs(catMotDerivs) > par.art.motThresh)) >0;

if (length(all_regs) ~= length(Art))
   error('regs do not match artefact vector') 
end

ArtWithNeighbors = unique([find(Art); find(Art)-1; find(Art)+1]);

regOnsets = find(all_regs);

idx.noArt = ~ismember(regOnsets, ArtWithNeighbors);
%% class index

%idxClass{1} = idx.cor==1;
%idxClass{1} = (idx.coh_==0);
idxClass{1} = ((idx.corWithZeros==1) .* (idx.resp_face==1) .* idx.noArt) == 1;
idxClass{2} = ((idx.corWithZeros==1) .* (idx.resp_house==1) .* idx.noArt) == 1;

mean(idx.noArt)

%allActsLogitDirectionOfCorrectSource = allActsLogit;
%allActsLogitDirectionOfCorrectSource(idxClass{2}) = -allActsLogit2(idxClass{2});

%allFaceActsDirectionOfCorrectSource = allActsLogitFace;

%allHouseActsDirectionOfCorrectSource = allActsLogitHouse;

%allNoiseActsDirectionOfCorrectSource = allActsLogitNoise;

%allClassAct = allFaceActsDirectionOfCorrectSource;
%allClassAct(idxClass{2}) = allHouseActsDirectionOfCorrectSource(idxClass{2});

%allActsWithoutNoiseLogitDirOfCorSource = allActsLogit;
%allActsWithoutNoiseLogitDirOfCorSource(idxClass{2}) = -allActsLogit2(idxClass{2});

%className{1} = 'percSeparateFaceAndSceneByUnsignedEv';
className{1} = 'percFace';
className{2} = 'percHouse';

S.patRegSet = {'percByConfAndRT'};


betaMapsSet = {  {'conf' 'RT'} };


%% Mean ROI analysis
% 
% if S.MeanROIAnalysis    
%     S.roi_name = 'sphere_5-63_-41_40.nii';
%     S.roi_file = ['/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data/group_analyses/conjAnalysis/' S.roi_name];
%     subj = load_spm_mask(subj,  S.roi_name, S.roi_file);
%     all_roi_pat = get_masked_pattern(subj, S.preprocPatCondensedName, S.roi_name);
%     
%     
%     for bS = 1:length(betaMapsSet)
%         S.patReg = S.patRegSet{bS};
%         S.patRegDir = [S.expt_dir S.subj_id '/' S.patReg];
%         for i = 1:length(idxClass)
%             pat{i} = all_roi_pat(:,idxClass{i});
%             sess{i} = allSess(:,idxClass{i});
%             acts{i} = allActsLogitDirectionOfCorrectSource(idxClass{i});
%             coh{i} = allCoh(idxClass{i});
%             conf{i} = allConf(idxClass{i});
%             class{i} = allClass(idxClass{i});
%             
%             datMat = {[actsDiff{i}' coh{i}]};
%             
%             X{i} = [datMat{bS} sess{i}'];
%             nRegRows = size(X{i},2) - size(sess{i},1);
%             XZ = [zscore(X{i}(:,1:nRegRows)) sess{i}']; %currently we just include the intercept term.
%             
%             mean_roi_pat{i} = mean(pat{i});
%             mean_roi_pat_z{i} = zscore(mean_roi_pat{i});
%             [b{i}, ~, ~, ~, stats{i}] = regress(mean_roi_pat_z{i}', XZ);
%             
%             
%             psyphys.(className{i}).y = mean_roi_pat_z{i}';
%             psyphys.(className{i}).x = XZ(:,1);
%             psyphys.(className{i}).betas = [b{i}(end) b{i}(1)];
%             psyphys.(className{i}).SE = [];
%             psyphys.(className{i}).ticks.x = [-2.5 0 2.5];
%             psyphys.(className{i}).ticks.y = [-3 3];
%             psyphys.(className{i}).xlabel = 'evidence';
%             psyphys.(className{i}).ylabel = 'activity' ;
%             %psyphys.(className{i}).xVals = -2:2;
%             psyphys.(className{i}).marker = 'o-';
%             
%             xValSet = unique(XZ(:,1));
%             for k = 1:length(xValSet)
%                 psyphys.(className{i}).xVals(k) = xValSet(k);
%                 psyphys.(className{i}).mean(k) = mean(mean_roi_pat_z{i}(XZ(:,1) == xValSet(k)));
%                 psyphys.(className{i}).SE(k) = ste( mean_roi_pat_z{i}(XZ(:,1) == xValSet(k)));
%             end
%         end
%         
%     end
% end

%% Map creation

if S.WriteRegressionMaps
    for bS = 1:length(betaMapsSet)
        S.patReg = S.patRegSet{bS};
        S.patRegDir = [S.expt_dir S.subj_id '/' S.patReg];
        for i = 1:length(idxClass)
            pat{i} = allPats(:,idxClass{i});
            sess{i} = allSess(:,idxClass{i});
            %acts{i} = allActsLogitDirectionOfCorrectSource(idxClass{i});
            coh{i} = allCoh(idxClass{i});
            conf{i} = allConf(idxClass{i});
            class{i} = allClass(idxClass{i});
            rt{i} = allRT(idxClass{i});
            respClass{i} = allRespClass(idxClass{i});
            
            %actsFace{i} = allFaceActsDirectionOfCorrectSource(idxClass{i});
            %actsHouse{i} = allHouseActsDirectionOfCorrectSource(idxClass{i});
            
            %classAct{i} = allClassAct(idxClass{i});
            %noiseAct{i} = allNoiseActsDirectionOfCorrectSource(idxClass{i});
            
            datMat = {[conf{i} rt{i}]};
            
            X{i} = [datMat{bS} sess{i}'];
            
            b{i} = nan(size(X{i} ,2), size(pat{i},1));
            
            nRegRows = size(X{i},2) - size(sess{i},1);
            
            XZ{i} = nan(size(X{i}));
            for j=1:nRegRows
                idxNotNan = ~isnan(X{i}(:,j));
                XZ{i}(idxNotNan,j) = zscore(X{i}(idxNotNan,j));
            end
            idxNonRegRows = (j+1):size(X{1},2);
            XZ{i}(:,idxNonRegRows) = X{i}(:,idxNonRegRows);
            
            %XZ{i} = [zscore(X{i}(:,1:nRegRows)) sess{i}'];
            
            %thisConf = XZ{i}(:,2);
            %thisCoh = XZ{i}(:,1);
            
            %             %% create fake coh
            %             coh_sorted = sort(thisCoh);
            %             rCohConf = corr(thisCoh,thisConf);
            %             k = .2;
            %             diffR = 1;
            %             nIts = 0;
            %             while (diffR>.025 && nIts<100)
            %                 nIts = nIts+1;
            %                 fakeCoh_h = thisConf + k*randn(size(thisConf));
            %                 [fC_sorted,fC_idx] = sort(fakeCoh_h);
            %                 [~,fC_idx_inv] = sort(fC_idx);
            %
            %                 fakeCoh = coh_sorted(fC_idx_inv);
            %
            %                 rFakeCohConf = corr(fakeCoh,thisConf);
            %
            %                 diffR = abs(rFakeCohConf-rCohConf);
            %
            %                 if diffR>0
            %                     k = k+.025;
            %                 else
            %                     k = k-.025;
            %                 end
            %             end
            %
            %             fakeCoh_z = zscore(fakeCoh);
            %%
            
            %            normal linear regression
            for j = 1:size(pat{i},1)
                thisPat = pat{i}(j,:)';
                thisPatZ = zscore(thisPat);
                [b{i}(:,j)] = regress(thisPatZ, XZ{i});
            end
        
            %% special analysis of effect of coherence + BOLD on confidence
            
            
%             for i = 1:length(pat)
%                 thisCoh = coh{i}';
%                 thisCohZ = zscore(thisCoh);
%                 
%                 thisConf = conf{i}';
%                 thisConfZ = zscore(thisConf);
%                 
%                 thisSess = sess{i}';
%                 
%                 xSize1 = size([ thisCohZ' thisSess] ,2) + 1;
%                 xSize2  = size(pat{i},1);
%                 b{i} = nan(xSize1, xSize2);
%                 stats{i}{k} = nan(xSize1, xSize2);
%                 for k = 1:size(pat{i},1)
%                     thisPat = pat{i}(k,:)';
%                     thisPatZ = zscore(thisPat);                    
%                     datMat = [thisPatZ thisCohZ' thisSess];
%                     
%                     [b{i}(:,k)] = regress(thisConfZ', datMat);
%                     
%                 end
             end
    end
        
        vol_info = S.vol_info;
        
        for p = 1:length(subj.patterns)
            patNames{p} = subj.patterns{p}.name;
        end
        
        voxel_inds = subj.masks{end}.mat;
        
        contrastName= {'acrossFacesAndHouses' 'FacesVsHouses'};
        imcalcText = {'(i1 + i2)/2' 'i1 - i2'};
        %%
        betaMaps = betaMapsSet{bS};
        for k = 1:length(betaMaps)
            for i = 1:length(pat)
                
                outMap=  zeros(vol_info.dim);
                
                outMap(voxel_inds==1) = b{i}(k,:);
                
                vol_info.dir = S.patRegDir;
                vol_info.fname = [ vol_info.dir '/' className{i} '_' betaMaps{k} '.img'];
                
                S.pat{i}.map{k}.fname = vol_info.fname;
                
                if isempty(dir([vol_info.dir]))
                    mkdir(vol_info.dir);
                end
                
                spm_write_vol(vol_info,outMap);
                
                fprintf('\n wrote out %s \n', vol_info.fname);
                
                d{i} = vol_info.fname;
            end
            
            dirChar = char(d);
            
            % print out contrasts
            for c = 1:length(contrastName)
                outputfileName = ['con_' betaMaps{k} '_' contrastName{c} '.img'];
                outputfile = [ vol_info.dir '/' outputfileName];
                spm_imcalc_ui(dirChar,outputfile, imcalcText{c})
            end
        end
end
end



