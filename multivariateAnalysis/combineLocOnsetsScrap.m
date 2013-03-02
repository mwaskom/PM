
c=0;
for C = [1:22]
    
    
    par1 = PM_Params(C, 'mnem');
    par2 = PM_Params(C, 'perc');
    if ~isempty(par2.substr)
        
        y1 = load(fullfile(par1.subdir, 'analysis_mvpa_Loc_3d/mvpa_ons.mat'));
        y2 = load(fullfile(par2.subdir, 'analysis_mvpa_Loc_3d/mvpa_ons.mat'));
        
        if ~ismember('noise', y1.names)
            y1.onsets{3} = [];
        end
        if ~ismember('noise', y2.names)
            y2.onsets{3} = [];
        end
        
        for i =1:3
            onsets{i} = [y1.onsets{i} y2.onsets{i} + sum(par1.numvols) * 2];
        end
        names = {'faceCor' 'houseCor' 'noise'};
        durations = {0 0 0};
    end
    
    
    analysisdir = fullfile(par1.subdir, 'analysis_mvpa_Loc_3d_2sess');
    
    if ~exist(analysisdir)
        mkdir(analysisdir)
    end
    cd(analysisdir);
    save mvpa_ons.mat onsets durations names;
    
end