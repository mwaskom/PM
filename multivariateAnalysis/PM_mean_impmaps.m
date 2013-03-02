function PM_mean_impmaps(S)


for it = 1:length(S.impType)
    for cds = 1:length(S.regNames)
        
        clear theMaps;
        cd ([S.importance_maps_dir '/' S.impType{it}]);
        
        for s = 1:length(S.subj_array)
            mapHelp = dir([S.subj_array{s} ['_' S.impType{it} '*' S.regNames{cds} '*img'] ]);
            %sHelp = dir([subjArray{s} '*scene*img']);
            
            theMaps_h{s} = mapHelp(1).name;
            %theMaps(s,:) = mapHelp(1).name;
            %FaceMaps(s,:) = fHelp(1).name;
            %SceneMaps(s,:) = sHelp(1).name;
            
        end
        theMaps = char(theMaps_h);
        
        spm_imcalc_ui(theMaps, ['mean' S.regNames{cds} 'ImpMap_smoothed.img'], 'sum(X)/size(X,1)', {'dmtx'});
        %spm_imcalc_ui(SceneMaps, 'meanSceneImpMap.img','sum(X)/size(X,1)', {'dmtx'});
        
        
    end
end