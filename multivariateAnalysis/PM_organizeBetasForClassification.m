function [subj, S] = PM_organizeBetasForClassification(subj, S)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

            PM_organizeBetasForClassification
            conds_orig=[1*ones(size(S.onsets{1})) 2*ones(size(S.onsets{2}))];
            o=[S.onsets{:}];
            [~,ixOnsetsSort]=sort(o);
            conds_h = conds_orig(ixOnsetsSort);
            
            if strcmp(S.thisSelector, 'TrainTestOneIterGroup')
                TrainTestOneIter = [1*S.idxOnsets_train_in_classifier 2*S.idxOnsets_test_in_classifier];
                actives =  TrainTestOneIter>0;
                
                conds = zeros(2,length(TrainTestOneIter));
                
                conds_h2 = zeros(size(TrainTestOneIter));
                conds_h2(actives) = conds_h;
                
                conds(1,conds_h2==1)=1;
                conds(2,conds_h2==2)=1;
                
                subj = init_object(subj,'selector','TrainTestOneIter');
                subj = set_mat(subj,'selector','TrainTestOneIter', TrainTestOneIter);
                subj = set_objfield(subj, 'selector', 'TrainTestOneIter', 'group_name', 'TrainTestOneIterGroup');
            elseif strcmp(S.thisSelector, 'randomNFold_xval')
                actives =  (S.idxOnsets_train_in_classifier | S.idxOnsets_test_in_classifier);
                train_actives = S.idxOnsets_train_in_classifier;
                
                randomNFold_h = ceil(shuffle(1:sum(actives))/(sum(actives)/S.nFolds));
                
                randomNFold = zeros(size(actives));
                randomNFold(actives) = randomNFold_h;
                
                conds_h2 = zeros(size(actives));
                conds_h2(actives) = conds_h;
                
                conds = zeros(2,length(actives));
                conds(1,conds_h2==1)=1;
                conds(2,conds_h2==2)=1;
                
                subj = init_object(subj,'selector','trainActives');
                subj = set_mat(subj,'selector','trainActives', train_actives);
                
                subj = init_object(subj,'selector','randomNFold');
                subj = set_mat(subj,'selector','randomNFold', randomNFold);
                subj = set_objfield(subj, 'selector', 'randomNFold', 'group_name', 'randomNFoldGroup');
                subj = PM_create_xvalid_indices_trainActivesOnly(subj,'randomNFold', 'actives_selname', 'trainActives');
            end
            
            
            
            
            subj = init_object(subj,'regressors','conds');
            subj = set_mat(subj,'regressors','conds',conds);
            subj = set_objfield(subj,'regressors','conds','condnames',S.condnames);
            
            S.classSelector = S.thisSelector;

end

