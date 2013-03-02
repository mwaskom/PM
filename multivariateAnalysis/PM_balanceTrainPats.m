function subj = PM_balanceTrainPats(S, subj)

sels = squeeze(get_group_as_matrix(subj, 'selector', S.thisSelector));
regs = get_mat(subj, 'regressors', 'conds');


% if strcmp(S.thisSelector, 'TrainTestOneIterGroup')
%     sels = sels';
% end
% 
if ~strcmp(S.thisSelector, 'randomNFold_xval')
    sels = sels';
end

for i = 1:S.numBalancedIts
    for s = 1:size(sels,1)
        selIt = (sels(s,:));
        trainRegIdx = find(selIt==1);
        testRegIdx = find(selIt==2);
        
        trainregs = regs(:,selIt==1);
        nTrainPats = selIt(selIt==1) * trainregs';
        
        newBinSize = min(nTrainPats(nTrainPats~=0));
        active_trials_train = zeros(1,size(trainregs,2));
        
        for t = 1:length(nTrainPats)
            if (nTrainPats(t)>0)
            theseTrials_h = shuffle(find((trainregs(t,:))==1));
            theseTrials = trainRegIdx(theseTrials_h);
            idxInclude{t} = theseTrials(1:newBinSize);
            active_trials_train(idxInclude{t}) = 1;
            end
        end
        
        allIncluded = sort([idxInclude{:}]);
        active_trials = zeros(1,size(selIt,2));
        active_trials(allIncluded) = 1;
        active_trials(testRegIdx) = 1;
%         NDiff = diff(nTrainPats);
%         active_trials = sum(regs);
%         
%         if diff(nTrainPats)<0 
%             theseTrials_h = shuffle(find((trainregs(1,:))==1));
%             theseTrials = trainRegIdx(theseTrials_h);
%             idxRemove = theseTrials(1:abs(NDiff));
%             active_trials(idxRemove) = 0;
%         elseif diff(nTrainPats)>0 
%             theseTrials_h = shuffle(find((trainregs(2,:))==1));
%             theseTrials = trainRegIdx(theseTrials_h);
%             idxRemove = theseTrials(1:abs(NDiff));
%             active_trials(idxRemove) = 0;
%         end
        
        selItName = [S.thisSelector num2str(s) '_balanced' num2str(i)];
        selItGroupName = [S.thisSelector 'balanced'];
        
        subj = init_object(subj,'selector', selItName);
        subj = set_mat(subj,'selector',selItName, selIt .* active_trials);
        subj = set_objfield(subj, 'selector', selItName, 'group_name', selItGroupName);
        
    end
end