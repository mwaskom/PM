function [S] = PM_subsample(S)
 allTrainOns = sort([S.onsets_train_in_classifier{:}]);
        allOns = sort([S.onsets{:}]);
        allTrainOnsFace = S.onsets_train_in_classifier{1};
        allTrainOnsHouse = S.onsets_train_in_classifier{2};
        
        idxTrialSubset = ismember(S.idxTr.alltrials, allTrainOns);
        idxTrialSubsetClass1 = ismember(S.idxTr.alltrials, S.onsets_train_in_classifier{1});
        idxTrialSubsetClass2 = ismember(S.idxTr.alltrials, S.onsets_train_in_classifier{2});
        
        idx = S.idxTr;

        idx.binaryConf = idx.unsignedConf==4;
        idx.toRemoveFaces = zeros(size(idx.face));
        idx.toRemoveHouses = zeros(size(idx.house));
            
        for j=1:4

            idxTheseFaces = (idx.face==1 & idx.unsignedConf==j & idxTrialSubset);
            idxTheseHouses = (idx.house==1 & idx.unsignedConf==j & idxTrialSubset);

            %idxTheseFaces = (idx.face==1 & idx.binaryConf==j & idxTrialSubset);
            %idxTheseHouses = (idx.house==1 & idx.binaryConf==j & idxTrialSubset);

            diffNTrials = sum(idxTheseFaces) - sum(idxTheseHouses);
            
            if diffNTrials>0
                theseFaceTrials = find(idxTheseFaces);
                toRemoveFaces_h = randperm(length(theseFaceTrials), diffNTrials);
                toRemoveFaces_h2 = theseFaceTrials(toRemoveFaces_h);
                idx.toRemoveFaces(toRemoveFaces_h2) = 1;
            elseif diffNTrials<0
                theseHouseTrials = find(idxTheseHouses);
                toRemoveHouses_h = randperm(length(theseHouseTrials), abs(diffNTrials));
                toRemoveHouses_h2 = theseHouseTrials(toRemoveHouses_h);
                idx.toRemoveHouses(toRemoveHouses_h2) = 1;
            end
        end

        
        
        if S.balanceHiAndLowConf
            NFaceHighConf = sum(idx.binaryConf .* ~idx.toRemoveFaces .* idxTrialSubsetClass1);
            NFaceLowConf = sum(~idx.binaryConf .* ~idx.toRemoveFaces .* idxTrialSubsetClass1);
            diffNTrialsBinConf = NFaceHighConf - NFaceLowConf;
            
            if diffNTrialsBinConf>0
                theseFaceHiConfTrials = find(idxTrialSubsetClass1 .* ~idx.toRemoveFaces .* idx.binaryConf);
                toRemoveFaces_h = randperm(length(theseFaceHiConfTrials), diffNTrialsBinConf);
                toRemoveFaces_h2 = theseFaceHiConfTrials(toRemoveFaces_h);
                idx.toRemoveFaces(toRemoveFaces_h2) = 1;
                
                theseHouseHiConfTrials = find(idxTrialSubsetClass2 .* ~idx.toRemoveHouses .* idx.binaryConf);
                toRemoveHouses_h = randperm(length(theseHouseHiConfTrials), diffNTrialsBinConf);
                toRemoveHouses_h2 = theseHouseHiConfTrials(toRemoveHouses_h);
                idx.toRemoveHouses(toRemoveHouses_h2) = 1;
                
            elseif diffNTrialsBinConf<0
                theseFaceLowConfTrials = find(idxTrialSubsetClass1 .* ~idx.toRemoveFaces .* idx.binaryConf);
                toRemoveFaces_h = randperm(length(theseFaceLowConfTrials), diffNTrialsBinConf);
                toRemoveFaces_h2 = theseFaceLowConfTrials(toRemoveFaces_h);
                idx.toRemoveFaces(toRemoveFaces_h2) = 1;
                
                theseHouseLowConfTrials = find(idxTrialSubsetClass2 .* ~idx.toRemoveHouses .* idx.binaryConf);
                toRemoveHouses_h = randperm(length(theseHouseLowConfTrials), diffNTrialsBinConf);
                toRemoveHouses_h2 = theseHouseLowConfTrials(toRemoveHouses_h);
                idx.toRemoveHouses(toRemoveHouses_h2) = 1;
            end
        end
        
        toRemoveFaces = ismember(find(idxTrialSubsetClass1), find(idx.toRemoveFaces));
        toRemoveHouses = ismember(find(idxTrialSubsetClass2), find(idx.toRemoveHouses));
        
        
        S.onsets_train_in_classifier{1}(toRemoveFaces) = [];
        S.onsets_train_in_classifier{2}(toRemoveHouses) = [];
        
        idxFinalFaceTrain = ismember(allTrainOns, S.onsets_train_in_classifier{1});
        idxFinalHouseTrain = ismember(allTrainOns, S.onsets_train_in_classifier{2});
        
        if sum(idxFinalFaceTrain)~=sum(idxFinalHouseTrain)
            error('subsampling error! classes are not balanced');
        end
        
        meanConfFaces = mean(idx.unsignedConf(logical(~idx.toRemoveHouses .* idxTrialSubsetClass2)));
        meanConfHouses = mean(idx.unsignedConf(logical(~idx.toRemoveFaces .* idxTrialSubsetClass1)));
        
        if S.balanceHiAndLowConf
            meanBinaryConfFaces = mean(idx.binaryConf(logical(~idx.toRemoveHouses .* idxTrialSubsetClass2)));
            meanBinaryConfHouses = mean(idx.binaryConf(logical(~idx.toRemoveFaces .* idxTrialSubsetClass1)));
            
            if meanBinaryConfFaces~=.5 || meanBinaryConfHouses~=.5
                error('subsampling error! different number of high vs. low confidence trials');
            end
        else
            meanConfFaces = mean(idx.unsignedConf(logical(~idx.toRemoveHouses .* idxTrialSubsetClass2)));
            meanConfHouses = mean(idx.unsignedConf(logical(~idx.toRemoveFaces .* idxTrialSubsetClass1)));
            
            if meanConfFaces~=meanConfHouses
                error('subsampling error! confidences are not equated across classes');
            end
        end
              
              
        sprintf('Subsampling performed.  %g per training bin', sum(idxFinalFaceTrain))

end

