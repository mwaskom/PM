function [] = PM_subUtil( par )
%unzip files that come from the CNI


% %% make a copy of the .nii files and prepend them with 'R'. 
% if exist(par.funcdir)
%     cd (par.funcdir);
% 
%     dr = dir('run*.nii');
%     dm = dir('mean.nii');
%     dv = dir('valid.nii');
%     
%     
% for i = 1:length(dr)
%     unix (['cp ' dr(i).name ' R' dr(i).name])
% end
%    
%     unix (['cp ' 'valid.nii' ' Rvalid.nii'])
%     unix (['cp ' 'mean.nii' ' Rmean.nii'])
% 
% else
%     warning(['no valid raw dir for subject ' par.substr]);
% end


% %%%% check for existence of preprocessed fucntional data. 
% if exist(par.funcdir)
%     cd (par.funcdir);
% 
%     dr = dir('Rrun*.nii');
%     
%     if isempty(dr)
%          warning(['no runs in functional dir for ' par.substr]);
%     end
% 
% else
%     warning(['no valid raw dir for subject ' par.substr]);
% end


% %%% move  
% if exist(par.rawdir)
%     cd (par.rawdir);
%     
%     dr = dir('run*.nii');
%     dm = dir('mean.nii');
%     dv = dir('valid.nii');
%     
%     if ~exist(par.funcdir)
%        mkdir(par.funcdir); 
%     end
%     
%     for i = 1:length(dr)
%         unix (['mv ' dr(i).name ' ' fullfile(par.funcdir, dr(i).name)]);
%     end
%     
%     unix (['mv ' dm.name ' ' fullfile(par.funcdir, dm.name)]);
%     unix (['mv ' dv.name ' ' fullfile(par.funcdir, dv.name)]);
%     
%     
% 
% else
%     warning(['no valid raw dir for subject ' par.substr]);
% end


%%% rename behavioral data in a consistent format. 
% if exist(par.behavdir)
%     cd (par.behavdir);
%     
%     %     dL = dir('*Loc*');
%     %
%     %     if isempty(dL)
%     %         fprintf('fixing loc filenames');
%     %         dl = dir('*loc*');
%     %         for i = 1:length(dl)
%     %            newName = strrep(dl(i).name, 'loc', 'Loc');
%     %            unix(['mv ' dl(i).name ' ' newName]);
%     %         end
%     %     end
%     %
%     %     dP = dir('*Perc*');
%     %
%     %     if isempty(dP)
%     %         fprintf('fixing perc filenames');
%     %         dp = dir('*perc*');
%     %         for i = 1:length(dp)
%     %            newName = strrep(dp(i).name, 'perc', 'Perc');
%     %            unix(['mv ' dp(i).name ' ' newName]);
%     %         end
%     %     end
%     %
%     dT = dir('AG4_test*');
%     if ~isempty(dT)
%         fprintf('fixing mnem filenames');
%         for i = 1:length(dT)
%             newName = strrep(dT(i).name, 'test', 'retrieve_');
%             newName = strrep(newName, 'out.mat', 'out\(1\).mat');
%             newName = strrep(newName, '_date', '');
%             unix(['mv ' dT(i).name ' ' newName]);
%         end
%     end
%     
% 
% else
%     warning(['no valid behav dir for subject ' par.substr]);
% end


% Determine the amount of perc and loc files in each subject
% if exist(par.behavdir)
%     cd (par.behavdir);
%     
%     dP = dir('*Perc*.mat');
%     dL = dir('*Loc*.mat');
%     hiddenFiles = dir('.*.mat');
%     
%     PercFileNames = setdiff({dP.name}, {hiddenFiles.name});
%     LocFileNames = setdiff({dL.name}, {hiddenFiles.name});
%     
%     w.P = 1:length(PercFileNames);
%     w.L = (1:length(LocFileNames)) + length(w.P);
% 
% else
%     warning(['no valid behav dir for subject ' par.substr]);
% end

% Remove poorly realigned data
% if exist(par.funcdir)
%     cd (par.funcdir);
%     
%     delete('sR*');
%     delete('R*');
% 
% else
%     warning(['no valid functional dir for subject ' par.substr]);
% end

% remove an analysis directory
% unix (['rm -r ' par.analysisdir]);


%% convert 4d to 3d niis

%store 4d files in their own directory.
% if ~exist(par.FourDNiiDir);
%     mkdir (par.FourDNiiDir)
% end
% 
% unix(['mv ' par.funcdir '/*run*.* ' par.FourDNiiDir]);

% cd (par.funcdir);
% for i=1:size(par.rascanfiles.all,1)
% 
%     im = par.rascanfiles.all(i,:);
%     
%     [pth,nm,ext] = fileparts(im);
% 
%     newfolder = pth;
% 
%     if newfolder(end) ~= filesep
%         newfolder(end+1) = filesep;
%     end
%     
%     thisRun = ['run' nm(end-1:end) '/'];
%     
%     if ~exist(fullfile(par.funcdir, thisRun))
%         mkdir(fullfile(par.funcdir, thisRun));
%     end
%     
%     V = spm_vol(im);
%     for i = 1:length(V)
%         Vnew = V(i);
%         Y = spm_read_vols(V(i));
% 
%         Vnew.fname = [par.funcdir,'/',thisRun,nm,num2str(i-1,'%04g'),'.nii'];
%         Vnew.n = [1 1];
%         
%         Vnew = spm_create_vol(Vnew);
%         spm_write_vol(Vnew,Y);
%         
%         %movefile(im, par.FourDNiiDir)
%     end    
% end



% remove analysisdir
%  unix(['rm -r ' par.analysisdir]);
% 
% rename behav files

%% rename behavioral data in a consistent format.
% if exist(par.behavdir)
%     cd (par.behavdir);
%     
%     dP = dir('*Perc*');
%     
%     if ~isempty(dP)
%         for i = 1:length(dP)
%             %             newName = strrep(dp(i).name, 'perc', 'Perc');
%             percIndex = strfind(dP(i).name, 'Perc');
%             FNLength = length(dP(i).name);
%             
%             if (FNLength - percIndex==4) || (FNLength - percIndex==8);
%                 thisNum = dP(i).name(percIndex+4);
%                 newNum = ['0' thisNum];
%                 
%                 newName = strrep(dP(i).name, ['Perc' thisNum], ['Perc' newNum]);
%                 unix(['mv ' dP(i).name ' ' newName]);
%             end
%         end
%     end
%     
% else
%     warning(['no valid behav dir for subject ' par.substr]);
% end

% d = dir([par.analysisdir '/pm*']);
% for i = 1:length(d)
%     cd (par.analysisdir);
%     newName = ['perc_coherence' d(i).name(end-3:end)] ;
%     unix(['mv ' d(i).name ' ' newName]);
% end


%%  normalize a group mask into subject-specific space
% groupMask = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/OTnoHipp.img';
% groupMaskDir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/OccTemp';
% 
% par.invNormalizationParams
% 
% % normalize the group mask into the subject's native space
% spm_write_sn(groupMask, par.invNormalizationParams, par.grwrflags);
% movefile(fullfile(groupMaskDir, 'wOccTemp7.nii'), fullfile(par.anatdir, 'nativeOccTemp.nii'));
% 
% 
% % threshold the native space group OT mask
% d{1} = fullfile(par.funcdir, 'Rmean.nii');
% d{2} = fullfile(par.anatdir, 'nativeOccTemp.nii');
% dirChar = char(d);
% %dirChar = fullfile(par.anatdir, 'nativeOccTemp.nii');
% imcalcText = 'i2>0.8 ';
% outputfile = fullfile(par.anatdir, 'tnativeOccTemp.nii');
% spm_imcalc_ui(dirChar,outputfile, imcalcText)
% 
% % mask group mask with individual subject's gray matter
% clear d dirChar
% d{1} = fullfile(par.anatdir, 'tnativeOccTemp.nii');
% d{2} = fullfile(par.anatdir, 'mask.nii');
%      
% dirChar = char(d);
% outputfile = fullfile(par.anatdir, 'tnativeOccTempGrey.nii');
% imcalcText = '(i1>0).*(i2>.2)';
% spm_imcalc_ui(dirChar,outputfile, imcalcText)


% %%% rename functionals  

% %% Select the top 500 voxels from each category
% %S.imgsDir = '/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/fmri_data2/mvpa_results/ImpMaps_11-Sep-2011/tsrmask23.hdr/pos';
% S.mask = '/Users/alangordon/Studies/AG1/Accumulator_fMRI/AG1/mvpa/OT_2012/OTnoHipp.img';
% 
% vMask = spm_vol(S.mask);
% volMask = spm_read_vols(vMask);
% 
% for s = 1:length(SA)
%     S.imgsDir = fullfile('/Users/alangordon/mounts/w5/alan/perceptMnemonic/fmri_data', SA{s}, '/analysis_retBinaryConf');
%     
%     vFace = spm_vol(fullfile(S.imgsDir, 'spmT_0001.img'));
%     %vScene = spm_vol(fullfile(S.imgsDir, [SA{s} '_pos_house.hdr']));
%     
%     volFace = spm_read_vols(vFace).* (volMask>.01);
%     %volScene = spm_read_vols(vScene).* (volMask>.01);
%     
%     vecFace = reshape(volFace, numel(volFace),1);
%     %vecScene = reshape(volScene, numel(volScene),1);
%     
% %     selectedVals = [];
% %     
% %         for i = 1:1000
% %            if(max(vecFace) >= max(vecScene))
% %                [voxVal, idxThisVox] = max(vecFace);
% %            else
% %                [voxVal, idxThisVox] = max(vecScene);
% %            end
% %            selectedVals(i) = voxVal;
% %            vecFace(idxThisVox) = [];
% %            vecScene(idxThisVox) = [];
% %        end
% %     
% %         voxThresh = selectedVals(500);
% 
%      vecFaceSorted = sort(vecFace);
%      %vecSceneSorted = sort(vecScene);
% %     
%      voxThresh(1) = vecFaceSorted(end - 499);
%      %voxThresh(2) = vecSceneSorted(end - 999);
% %     
%     d{1} = fullfile(S.imgsDir, 'spmT_0001.img');
%     %d{2} = fullfile(S.imgsDir, [SA{s} '_pos_house.hdr']);
% 
%     d{2} = S.mask;
%     
%     dirChar = char(d);
%     
%     outputfile = fullfile(S.imgsDir, ['/PConfVsSConf_OTNoHipp_500vox.img']);
%     
%     %imcalcText = ['((i1 > ' num2str(voxThresh(1)) ')+ (i2 > ' num2str(voxThresh(2)) ')) .* i3'];
%     imcalcText = ['(i1 >= ' num2str(voxThresh(1)) ') .* i2'];
%     %imcalcText = ['((i1 > ' num2str(voxThresh) ')+ (i2 > ' num2str(voxThresh) ')) .* i3'];
%     
%     spm_imcalc_ui(dirChar,outputfile, imcalcText)
% end    
    
%% normalize and then reslice

%       S.imgsDir = fullfile('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data', par.substr, '/anat');
% % 
% % %     imgToBeNormalized = fullfile(S.imgsDir, 'c2V001.nii');
% % %
% % %     graywrimg = par.graywrimg;
% % %     grwrflags = par.grwrflags;
% % %     [grpth,grnm] = fileparts(graywrimg);
% % %     grmatname        = fullfile(grpth,[grnm '_sn.mat']);
% % %     spm_write_sn(imgToBeNormalized, grmatname,grwrflags);
% % 
%      d{1} = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/pm_050812/functionalNoLPF/run01/Rrun010107.nii';
%      d{2} = fullfile(S.imgsDir, 'c2V001.nii');
% % 
%      flags.mask = 0;
%      flags.mean = 0;
%      flags.which = 1;
%      flags.interp = 1;
%      flags.prefix = 'rnative';
% 
%     spm_reslice(d,flags)

    % if exist(par.rawdir)
%     cd (par.rawdir);
%     
%     dr = dir('run*.nii');
%     dm = dir('mean.nii');
%     dv = dir('valid.nii');
%     
%     if ~exist(par.funcdir)
%        mkdir(par.funcdir); 
%     end
%     
%     for i = 1:length(dr)
%         unix (['mv ' dr(i).name ' ' fullfile(par.funcdir, dr(i).name)]);
%     end
%     
%     unix (['mv ' dm.name ' ' fullfile(par.funcdir, dm.name)]);
%     unix (['mv ' dv.name ' ' fullfile(par.funcdir, dv.name)]);
%     
%     
% 
% else
%     warning(['no valid raw dir for subject ' par.substr]);
% end


%     cd (par.analysisdir);
%     
%     dr = dir('trial*');
%       
%     for i = 1:length(dr)
%         unix (['mv ' dr(i).name ' ' fullfile(par.subdir, 'trialBetas/', ['spm_' dr(i).name])]);
%     end
%     
%     
%     

%% Mask a leave-one-out group contrast to isolate a single ROI
% %conImgs = {'conj_acc_p0.005.img' 'conj_conf_p0.005.img' 'conj_class_p0.005.img'	'conj_rt_p0.005.img'};
% rois = {'ts_leftAI.nii'	'ts_leftMFG.nii'	'ts_rightAI.nii'};
% %rois = {'ts_PCC.nii'	'ts_left_AnG.nii'	'ts_left_STG.nii'};
% 
% %for c=1:length(conImgs)
%     for r=1:length(rois)
%         d{1} = fullfile(par.mnem.analysisdir, 'percMnemConj_parModByConfAndRT_16Subs_1LeftOut', 'conj_RT_p.005_acrossPercAndMnem.img');
%         d{2} = fullfile('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/group_analyses/percMnemConj_parModByConfAndRT_16Subs/ROIs_RT', rois{r});
%         
%         dirChar = char(d);
%         outputfile = fullfile(par.mnem.analysisdir, 'percMnemConj_parModByConfAndRT_16Subs_1LeftOut', ['RT_p005_' rois{r}(4:end)]);
%     
%         imcalcText = '(i1>0).*(i2>0)';
%         spm_imcalc_ui(dirChar,outputfile, imcalcText)
%     end
% %end



%% Move old files into a special directory

% cd (par.subdir)
% 
% newDir = fullfile(par.subdir, 'oldAnalyses');
% if ~exist(newDir, 'dir');
%    mkdir(newDir); 
% end
% 
% dP = dir('patLin*');
% dA = dir('analysis*');
% dB = dir('percBy*');
% 
% dT = [dP; dA; dB];
% 
% for i=1:length(dT)
%     if isempty(strfind(dT(i).date, '2013'))
%         movefile(dT(i).name, newDir);
%     end    
% end

%% Move files back from special directory to initial directory

cd (fullfile(par.subdir, 'oldAnalyses'));

movefile('analysis_mvpa_mnemDM_3d', par.subdir);


