function [subj S] = PM_model_voxel_interactions(subj, S)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

origPat = single(get_masked_pattern(subj,S.classifier_pattern, S.classifier_mask));
idxOrigPat = find(get_mat(subj,'mask', S.classifier_mask));

patsToInteract = single(get_masked_pattern(subj,S.classifier_pattern, S.roi_name));
idxPatsToInteract = find(get_mat(subj,'mask', S.roi_name));

groups_h = get_mat(subj,'regressors', S.regName); 
groups = ones(1,size(groups_h,2));

idxTrain = squeeze(get_group_as_matrix(subj,'selector', S.classSelector))==1;

if ~isempty(S.portion)
    part = S.portion;
    idxThisTrain = idxTrain(part,:)';
end

for i=size(groups_h,1)
    idxTheseTrials = logical(groups_h(i,:));
    groups(idxTheseTrials) = i;
end

trainGroups = groups.*idxThisTrain';

pThresh = S.intPThresh;

a = subj.patterns{end}.mat;
a2 = sort(a);
pThresh = a2(1000)/10;

szP = size(origPat);

if S.interactionType==1 % interaction of all voxels in a pattern with a single ROI.
    newPats = cell(szP(1),1);
    for vx = 1:szP(1)
        newPats_h = repmat(origPat(vx,:),szP(1),1) .* origPat;
        newPats{vx} = newPats_h((vx+1):szP(1),:);
    end
    newPats = vertcat(newPats{:});
    newPats = vertcat(newPats,origPat);
    
    subj.patterns{end}.mat = [];
    clear origPat
    
    subj = duplicate_object(subj,'pattern',[S.preprocPatName],[S.preprocPatName '_conj']);
    subj = set_mat(subj,'pattern',[S.preprocPatName '_conj'],newPats,'ignore_diff_size',true);
    clear newPats
    
    S.preprocPatCondensedName = [S.preprocPatName '_conj'];
elseif S.interactionType==2 % interaction of all voxels among each other
    
    intVox.idx = [];
    intVox.dat = [];
    
    reverseStr = '';
    
    if S.intEst % perform feature selection on interaction features
        
        %important for dividing up the labor among different compute jobs.
        vxIts = 1:szP(1);
%         if isempty(S.portion)
%             vxIts = 1:szP(1);
%         else
%             vxInit = round((S.portion(1)-1)*szP(1)/S.portion(2)+1);
%             vxFinal = round((S.portion(1))*szP(1)/S.portion(2));
%             vxIts = vxInit:vxFinal;
%         end
        
        for vx = vxIts
            
            if mod(vx,S.intReportIncrement) == 0
                msg = sprintf('%d of %d',vx,szP(1));
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            
            %thisVol = origPat .* repmat(origPat(vx,:), szP(1),1);
            if strcmp(S.multiVoxMethod, 'interaction')
                thisVol = patsToInteract .* repmat(origPat(vx,:), size(patsToInteract,1),1);
            elseif strcmp(S.multiVoxMethod, 'corr')
                thisVol = corr(patsToInteract', origPat(vx,:)');
            end
            
            toExclude = ismember(idxPatsToInteract, idxOrigPat(1:vx));
            
            
            %ttest2 for the case in which there are two groups.  For multiple
            %groups, implement an anova call (TODO).
            [~,pMat] = ttest2(thisVol(:,trainGroups==1)',thisVol(:,trainGroups==2)');
            
            pIdx = find((pMat<=pThresh).*~toExclude');
            
            intVox.idx = [intVox.idx, [repmat(vx,size(pIdx)); pIdx]];%in coordinates of origPats x patsToInteract
            intVox.dat = [intVox.dat; thisVol(pIdx,:)] ;
            
        end
        if ~isempty(S.portion)
            saveStr = [S.intFilePrefix prepend(S.portion(1)) 'Of' prepend(S.portion(2))];
            save(fullfile(S.workspace_dir, saveStr), S.intFilePrefix);
        end
    end
    
    if S.intConcat %concatenate interaction features
        intParts = dir(fullfile(S.workspace_dir, [S.intFilePrefix '*Of*']));
        
        if ~isempty(intParts)
            for p=1:length(intParts)
                r = load(fullfile(S.workspace_dir,intParts(p).name));
                intIdx_h{p} = r.intVox.idx;
                intDat_h{p} = r.intVox.dat;
            end
            intVox.idx = horzcat(intIdx_h{:});
            intVox.dat = vertcat(intDat_h{:});
            
            clear intIdx_h intDat_h;
            delete(fullfile(S.workspace_dir, [S.intFilePrefix '*Of*'])); %Change this. we dont' want to delete everything. also, have a name in the params
            save (fullfile(S.workspace_dir, S.intFilePrefix), 'intVox');
            
        else
            ldDat = load(fullfile(S.workspace_dir, S.intFilePrefix) );
            intVox = ldDat.intVox;
        end
    end
    
    if S.intUseIntsWithUniqueInfo
        newDat = nan(size(intVox.dat));
       
        for i=1:size(intVox.dat,1)
            if intVox.idx(1,i) ~= intVox.idx(2,i)
                X(1,:) = origPat(intVox.idx(1,i),:);
                X(2,:) = patsToInteract(intVox.idx(2,i),:);
                X(3,:) = intVox.dat(i,:);
                XOrth = spm_orth(X');
                [~, ~, stats] = mnrfit(XOrth(idxTrain,3),groups(idxTrain)');
                
                P(i) = stats.p(2);
                
                newDat(i,:) =  XOrth(:,3)';
             
            end
        end
        if length(P)<1000
            intVox.dat = [newDat];
        else
            PSorted = sort(P);
            pThresh2 = PSorted(1000);
            newFeat = newDat((PSorted<=pThresh2),:);
            intVox.dat = [newFeat];
            
        end
        size(intVox.dat)
    end
    
    %put pattern in subj
    subj = duplicate_object(subj,'pattern',S.preprocPatNameFinalMask,S.intPatName);
    subj = set_mat(subj,'pattern',S.intPatName,intVox.dat,'ignore_diff_size',true);
    subj = set_objfield(subj, 'pattern', S.intPatName, 'group_name', S.intGroupName);
    subj = set_objfield(subj, 'pattern', S.intPatName, 'masked_by', S.intMaskName);
    
    %make a dummy mask
    subj = duplicate_object(subj,'mask',subj.masks{1}.name, S.intMaskName);
    subj.masks{end}.mat = ones(size(intVox.dat,1),1);
    subj.masks{end}.matsize = size(subj.masks{end}.mat);
    %subj = set_mat(subj,'mask','interactionMask',ones(size(intVox.dat,2),1), 'ignore_diff_size',true);
    
    S.classifier_mask = S.intMaskName;
    S.classifier_pattern = S.intGroupName;
    
end
end


