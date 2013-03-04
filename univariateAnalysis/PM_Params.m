function par = PM_Params(subNo, task, loadScans)
% sets up parameters for batching of PM experiment.
% <subNo>: a value in the range of 1:22, or the sub_id (e.g. 'pm_042912')
% <task>: 'mnem' or 'perc'
% <load scans>: 
% 1 - load all the scans (this is necessary for analyses of the fmri data).
% 0 - do not load scans (this saves time)

%% list of subjects in the present experiment
percSubs = {'pm_031711' 'pm_042211' 'pm_040711' 'pm_050111_2' '' 'pm_050111' 'pm_050811_2' 'pm_051111_2' '' 'pm_052611' 'pm_051211' ...
    'pm_051411'  '' 'pm_052211' 'pm_052311' 'pm_052511' '' 'pm_061011' 'pm_041512' 'pm_042912' 'pm_050612' 'pm_051312'};

mnemSubs = {'pm_052311_2' 'pm_032811' 'pm_040611' 'pm_042511' '' 'pm_042811' 'pm_050211' 'pm_050411' 'pm_050411_2' 'pm_050511' ...
'pm_051111' 'pm_051211_2'  'pm_051511' 'pm_051911' 'pm_052011' 'pm_052211_2' 'pm_060111' 'pm_060711' 'pm_040312' 'pm_041812' ...
'pm_042412' 'pm_050812'};


%% basic task information
par.subNo = subNo;
par.task = task;
par.subTask = 'DM';

%% convert input
if (nargin > 2)
    par.loadScans = loadScans; 
else
    par.loadScans = 1; % by default, load scans
end

% if subNo is given as a string, convert to subNo
if ischar(subNo) || iscellstr(par.subNo)
    if strcmp(par.task, 'perc')
        par.subNo = find(strcmp(percSubs, par.subNo));
    elseif strcmp(par.task, 'mnem')
        par.subNo = find(strcmp(mnemSubs, par.subNo));
    end
else
    par.subNo = subNo;
end



%% subject-specific information

% each subject has slightly different acquisition information.  Here is
% where we store that information.

%par.str.perc = name of perc dataset
%par.str.mnem = name of mnem dataset
%par.scansSelect.perc.DM = indices of perceptual decision scans
%par.scansSelect.perc.loc = indices of localizer scans
%par.perc.FM2Funcs = which fieldmap corresponds to which perc scan?
%par.scansSelect.mnem.DM = indices of mnemonic decision scans
%par.scansSelect.mnem.loc = indices of localizer decision scans
%par.mnem.FM2Funcs = which fieldmap crresponds to which mnemonic scan?
%par.mnem.numvols = # of volumes collected during mnem session
%par.perc.numvols = # of volumes collected during perc session
%par.goodSub = is this subject included in our analysis?
%par.flagIt = should we flag this subject for some special reason?

switch par.subNo
    case 1
        par.str.perc = 'pm_031711'; %name of perc dataset
        par.str.mnem = 'pm_052311_2'; %name of mnem dataset
        par.scansSelect.perc.DM = 1:6; %indices of perceptual decision scans
        par.scansSelect.perc.loc = 7:9; %indices of localizer scans
        par.perc.FM2Funcs = []; %which fieldmap corresponds to which perc scan?
        par.scansSelect.mnem.DM = 1:3; % indices of mnemonic decision scans
        par.scansSelect.mnem.loc = 4:6; % indices of localizer decision scans   
        par.mnem.FM2Funcs = []; % which fieldmap crresponds to which mnemonic scan?
        par.mnem.numvols = [250 250 250 176 176 176]; % # of volumes collected during mnem session
        par.perc.numvols = [176 176 176 176 176 176 176 176 176]; % # of volumes collected during perc session
        par.goodSub = 1; % is this subject included in our analysis?
        par.flagIt = 0; % should we flag this subject for some special reason?
        %notes: no fieldmaps collected for this subject.
        %RT not recorded for '9' button presses for this subject.
    case 2
        par.str.perc = 'pm_042211';
        par.str.mnem = 'pm_032811';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 10:12;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 1 1 1 2 2 2 2 2 3 3];
        par.mnem.numvols = [154 254 254 179 179 179];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
        par.extremeMotion = 1;
    case 3
        par.str.perc = 'pm_040711';
        par.str.mnem = 'pm_040611';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 10:13;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 2 2 3 3 4 4 5 5 6 6 6 6];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
    case 4
        par.str.perc = 'pm_050111_2';
        par.str.mnem = 'pm_042511';
        par.scansSelect.perc.DM = 2:11;
        par.scansSelect.perc.loc = 12:14;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [0 176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 0;
        par.flagIt = 1;
    case 5
        par.str.perc = '';%'pm_042811_2';
        par.str.mnem = '';%'pm_042711';
        par.scansSelect.perc.DM = 1:7;
        par.scansSelect.perc.loc = 8:9;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 1 1 1 1 1];
        par.perc.FM2Funcs = [1 1 2 2 2 2 2 2 2];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176];
        par.goodSub = 0;
        par.flagIt = 0;
        par.extremeMotion = 1;
        %Hi-res SPGR collected during perc but not mnem.  Thus, hi-res from
        %perc section will be used for mnem data
        
        %data will be thrown out, and thus not analyzed.  This is a bad
        %subject for multiple reasons.
    case 6
        par.str.perc = 'pm_050111';
        par.str.mnem = 'pm_042811';
        par.scansSelect.perc.DM = 1:10;
        par.scansSelect.perc.loc = 11:13;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 4 4 5 5 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 1;
        par.extremeMotion = 1;
    case 7
        par.str.perc = 'pm_050811_2';
        par.str.mnem = 'pm_050211';
        par.scansSelect.perc.DM = [1:4 6:8 10];
        par.scansSelect.perc.loc = 11:13;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [2 2 3 3 4 4 4 5 5 6 6];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 0 176 176 176 0 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 1;
        par.extremeMotion = 1;
        % perc runs 5 and 9 were aborted; will not be included in the data set.
    case 8
        par.str.perc = 'pm_051111_2';
        par.str.mnem = 'pm_050411';
        par.scansSelect.perc.DM = 1:10;
        par.scansSelect.perc.loc = 11:13;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 1;        

    case 9
        par.str.perc = '';
        par.str.mnem = 'pm_050411_2';
        par.scansSelect.perc.DM = [];
        par.scansSelect.perc.loc =[];
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [];
        par.goodSub = 0;
        par.flagIt = 0;
    case 10
        par.str.perc = 'pm_052611';
        par.str.mnem = 'pm_050511';
        par.scansSelect.perc.DM = 1:10;
        par.scansSelect.perc.loc = 11:13;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 1;
        % buttons coded incorrectly during perc
    case 11
        par.str.perc = 'pm_051211';
        par.str.mnem = 'pm_051111';
        par.scansSelect.perc.DM = 1:8;
        par.scansSelect.perc.loc = 9:11;
        par.perc.FM2Funcs = [];
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
    case 12
        par.str.perc = 'pm_051411';
        par.str.mnem = 'pm_051211_2';
        par.scansSelect.perc.DM = 1:10;
        par.scansSelect.perc.loc = 11:13;
        par.perc.FM2Funcs = [];
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 1;
    case 13
        par.str.perc = '';
        par.str.mnem = 'pm_051511';
        par.scansSelect.perc.DM = [];
        par.scansSelect.perc.loc = [];
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [];
        par.goodSub = 0;
        par.flagIt = 0;
    case 14
        par.str.perc = 'pm_052211';
        par.str.mnem = 'pm_051911';
        par.scansSelect.perc.DM = 1:10;
        par.scansSelect.perc.loc = 11:13;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 1;
    case 15
        par.str.perc = 'pm_052311';
        par.str.mnem = 'pm_052011';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 10:12;
        par.perc.FM2Funcs = [];
        par.scansSelect.mnem.DM = [1 2 6];
        par.scansSelect.mnem.loc = 3:5;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 5 5];
        par.mnem.numvols = [250 250 176 176 176 250];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
        par.extremeMotion = 1;
        %first mnemonic run was repeated at the end.
        %the original first mnemonic run is now treated as if it never
        %existed.  the fmri file for this run is buried in a subdirectory.
        %the first behavioral file used corresponds to retrieval run4, the
        %second to run2 and the third to run3.  This is taken care of in
        %Mnemonic_fMRIBehAnalysis_Retrieval.
    case 16
        par.str.perc = 'pm_052511';
        par.str.mnem = 'pm_052211_2';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 11:13;
        par.perc.FM2Funcs = [];
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 0 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
        
    case 17
        par.str.perc = '';
        par.str.mnem = 'pm_060111';
        par.scansSelect.perc.DM = [];
        par.scansSelect.perc.loc = [];
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [];
        par.goodSub = 0;
        par.flagIt = 0;
        %second mnemonic run was aborted and repeated.
    case 18
        par.str.perc = 'pm_061011';
        par.str.mnem = 'pm_060711';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 10:12;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
    case 19
        par.str.perc = 'pm_041512';
        par.str.mnem = 'pm_040312';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 10:12;
        par.scansSelect.mnem.DM = [1 2 3 5];
        par.scansSelect.mnem.loc = [4 6 7];
        par.mnem.FM2Funcs = [1 1 2 2 2 3 3];
        par.perc.FM2Funcs = [2 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [154 234 254 176 114 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
    case 20
        par.str.perc = 'pm_042912';
        par.str.mnem = 'pm_041812';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 10:12;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 4];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
    case 21
        par.str.perc = 'pm_050612';
        par.str.mnem = 'pm_042412';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 10:12;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
    case 22
        par.str.perc = 'pm_051312';
        par.str.mnem = 'pm_050812';
        par.scansSelect.perc.DM = 1:10;
        par.scansSelect.perc.loc = 11:13;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
    otherwise
        error('unrecognized subject');
end

par.substr = par.str.(par.task);  
subject = par.str.(task);

%% select which scans to include
par.scansSelect.perc.all = sort([par.scansSelect.perc.DM par.scansSelect.perc.loc]);
par.scansSelect.mnem.all = sort([par.scansSelect.mnem.DM par.scansSelect.mnem.loc]);
par.scans_to_include = par.scansSelect.(par.task).(par.subTask);
par.FM2Funcs = par.(task).FM2Funcs;
par.numvols = par.(task).numvols(par.scans_to_include);

%% directory information
par.exptdir = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/';
par.scriptsdir = fullfile(par.exptdir,'scripts');
par.fmridir = fullfile(par.exptdir, 'fmri_data');
par.subdir = fullfile(par.fmridir, subject);
par.anatdir = fullfile(par.subdir, 'anat');
par.funcdir = fullfile(par.subdir, 'functionalNoLPF');
par.behavdir = fullfile(par.subdir, 'behav');
par.rawdir = fullfile(par.subdir, 'raw');
par.rawNiiDir = fullfile(par.subdir, 'raw/niiArchives');
par.niiDir = fullfile(par.subdir, 'rawNiis');
par.logdir = fullfile(par.subdir, 'logfiles');
par.artrepdir = fullfile(par.subdir, 'artRepairDM');
par.classmat = fullfile(par.exptdir, '/fmri_data/mvpa_files/trainTestRet_balancedConf.mat');
par.meanfuncdir = fullfile(par.funcdir, 'meanFuncs');
%par.analysisdir = fullfile(par.subdir, 'analysis_retBinaryConf');
par.analysisdir = fullfile(par.subdir, 'AnalysisRetByRREv');
[~, par.thisAnalysis] = fileparts(par.analysisdir);
par.perc.analysisdir = fullfile(par.subdir, 'analysis_perc_parmodByConfAndRT_16Subs_leftOut');
par.mnem.analysisdir = fullfile(par.subdir, 'analysis_mnem_parmodByConfAndRT_16Subs_leftOut');
par.funcs3d = fullfile(par.subdir, 'funcs3d');
par.FourDNiiDir = fullfile(par.funcdir, '4dniis');

%% exceptions to directory information
if strcmp(par.substr, 'pm_042211');
     par.rawdir = fullfile(par.subdir, 'raw_copy');
end

%% anatomy info
dHiRes = dir(fullfile(par.rawdir, '0*FSPGR*'));
dInplane = dir(fullfile(par.rawdir, '0*Inplane*'));
par.hiresDir = fullfile(par.rawdir, dHiRes.name); 
par.inplaneDir = fullfile(par.rawdir, dInplane.name); 
par.inplaneimg = [par.anatdir '/In001.nii'];
par.hiresimg = [par.anatdir '/V001.nii'];

%% scan acquisition
par.numscans = sum(par.numvols); %number of scans
par.TR = 2; %TR in seconds
par.numslice = 36;
par.TA = par.TR-par.TR/par.numslice; %TA
par.sliceorder = 1:1:par.numslice; %slice order (assumes ascending)
par.refslice = floor(par.numslice/2); %reference slice (assumes middle)
par.dropvol = 6; %dropped volumes (Assumed at start of scan)
par.minvol = par.dropvol+1; %first non-dropped volume
par.slicetiming(1) = par.TA / (par.numslice -1);%timing var for slice timing
par.slicetiming(2) = par.TR - par.TA; %timing var for slice timing

%% fieldmap correction
par.FMCorrect = 1;

%% artifacts
par.art.motThresh = .5; % motion artefact threshold
par.art.sigThresh = 4; % signal intensity artefact threshold

%% realignment...
par.realflag.quality = 0.9;
par.realflag.fwhm = 5;
par.realflag.sep = 4;
par.realflag.rtm = 1;
par.realflag.interp = 2;

%% reslicing
par.reslflag.mask = 1;
par.reslflag.mean = 1;
par.reslflag.interp = 4;
par.reslflag.which = 0;

%% re-orienting
par.reorientmat = [105.0000 105.0000 59.4000 3.1416 0 -1.5708 1.0000 1.0000 1.0000 0 0 0]; 
%the proper rotation and translation to get the functional image preprocessed by kendrick's code to match the anatomical
dRN = dir([par.funcdir '/R*.nii']);
rawNiis_h = {};
for i = 1:length(dRN), 
    rawNiis_h{i}=fullfile(par.funcdir, dRN(i).name); %#ok<AGROW> %where the preprocessed .niis are stored
end 
par.rawNiis = char(rawNiis_h{:});

%% coregistration
par.cor_meanfunc = [par.funcdir '/Rmean.nii']; %can change this to make more robust!
par.cor_inimg = [par.anatdir '/In001.nii'];
par.cor_hiresimg = [par.anatdir '/V001.nii'];

%% segmentation
par.img2bSeg = par.cor_hiresimg;
par.segopts = struct('biascor',1,'GM',[0 0 1],'WM',[0 0 1],'CSF',[0 0 0],'cleanup',0);
par.segbiasfwhm = 60; % 60 is default in gui, 75 is default in command line for reasons unexplained

%% normalization:
par.graytemp = '/Applications/spm8/apriori/grey.nii';
par.grsegs(1,:) = [par.anatdir '/c1' 'V001.nii']; 
par.grsegs(2,:) = [par.anatdir '/c2' 'V001.nii']; 
par.graysrcimg = [par.anatdir '/c1' 'V001.nii']; 
par.graywrimg = [par.anatdir '/c1' 'V001.nii']; 
par.normalizationParams = [par.anatdir '/V001_seg_sn.mat'];
par.invNormalizationParams = [par.anatdir '/V001_seg_inv_sn.mat'];
par.grflags.smosrc = 8;
par.grflags.smoref = 0;
par.grflags.regtype = 'mni';
par.grflags.cutoff = 25;
par.grflags.nits = 16;
par.grflags.reg = 1;
par.grwrflags.preserve = 0;
par.grwrflags.bb = [[-78 -112 -50];[78 76 85]]; % where did this come from?
par.grwrflags.vox        = [3 3 3];
par.grwrflags.interp     = 1;
par.grwrflags.wrap       = [0 0 0];

%% spgm normalization:
par.spgrnrmflags.preserve = 0;
par.spgrnrmflags.bb = [[-78 -112 -50];[78 76 85]];
par.spgrnrmflags.vox = [2 2 2];
par.spgrnrmflags.interp     = 1;
par.spgrnrmflags.wrap       = [0 0 0];

%% smoothing kernels
par.s4mm_smoothkernel = [4 4 4];
par.s8mm_smoothkernel = [8 8 8];
par.s6mm_smoothkernel = [6 6 6];

%% specmaskvars
par.specwrflags.preserve = 0;
par.specwrflags.bb = [[-78 -112 -50];[78 76 85]];
par.specwrflags.vox        = [1 1 1];
par.specwrflags.interp     = 1;
par.specwrflags.wrap       = [0 0 0];
par.specsmooth = [3 3 3];%changed from [20 20 20] on 071409;

%% general mask vars
par.maskimg = [par.anatdir '/mask.nii']; 
par.smaskimg = [par.anatdir '/smask.nii'];
par.tsmaskimg = [par.anatdir '/tsmask.nii'];
par.wmaskimg = [par.anatdir '/wmask.nii'];
par.swmaskimg = [par.anatdir '/swmask.nii'];
par.tswmaskimg = [par.anatdir '/tswmask.nii'];
par.twmaskimg = [par.anatdir '/twmask.nii'];
par.addsegs = 'i1 + i2';
par.maskthresh = 'i1 > .6';

%% denoising
par.denoiseOpt.denoisespec = '10001';
par.denoiseDur = .5;
par.denoisingConds = {{'faceCor'} {'houseCor'}};

%% beta series estimation
par.betaDir = fullfile(par.subdir, 'betaSeries');
par.betaPrefix = 'retTrial';
par.betaAnalysisDir = fullfile(par.subdir, 'betaTmp');
par.betaEstimationAlg = 'glmEst';

%% model vars...
par.timing.fmri_t0 = par.refslice;  %micro-time resolution stuff (changed from 8)
par.timing.fmri_t = par.numslice; %used to be 16
par.timing.units = 'secs';
par.bases.hrf.derivs = [0 0]; 
par.volt = 1;%if you don't want to model volterra interactions (and you don't want to model volterra interactions) leave this at 1.
par.sess.multi = {fullfile(par.behavdir, 'ons.mat')};
par.sess.multi_reg = {fullfile(par.behavdir, 'regs.mat')};
par.sess.hpf = 128;  
par.cvi = 'AR(1)'; %note that this actually gets changed to AR(0.2) in spm_fmri_spm_ui.  
par.global = 'None';
par.mask = {'/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/groupMask/inclusive_mask.img'};% explicit mask

%% contrast vars
par.constat = 'T';



%% files

warning('off');
if ~isempty(par.substr) && par.loadScans % if we are supposed to load the scans, load them
    [par.betas par.betasByRun] = findScans(par, '*Trial*.nii', 3, par.betaDir);
    [par.rascanfiles par.rascanfilesByRun] = findScans(par, 'Rrun*.nii',3, par.funcdir);
    [par.wrascanfiles par.wrascanfilesByRun] = findScans(par, 'wRrun*.nii',3, par.funcdir);
    [par.swrascanfiles par.swrascanfilesByRun] = findScans(par, 'swRrun*.nii',3, par.funcdir);
else % otherwise, don't load the scans.
    par.rascanfiles = [];
    par.rascanfilesByRun = [];
    par.swrascanfiles = [];
    par.swrascanfilesByRun = [];
end
warning('on');
end



%%
function [scanfiles scansByRun] = findScans(par, matchstring, scanDim, scandir)
% load up relevant scans, using 'matchfiles'

% input:
% <par> - the relevant par structure
% <matchstring> - a string representing a class of image data of interest, e.g. 'Rrun*.nii'
% <scanDim> - how many dimensions is the data? (e.g. 3)
% <scanDir> - where the functional data are located

% returns: 
% <scanfiles> - lists all the images of interest, for the specified
% <scansByRun> - lists all scan files of interest by run


scanfiles = [];
scansByRun = [];

clear sf_h sf_h2
if ~isempty(par.scansSelect.(par.task).(par.subTask))
    
    %% scanfiles for the specified subtask
    if scanDim==3 % for 3d files   
        for i=par.scansSelect.(par.task).(par.subTask);
            ix = find(par.scansSelect.(par.task).all==i);
            scans_h{ix} = matchfiles (fullfile(scandir, ['run' prepend(ix) '/' matchstring ])); %#ok<AGROW>
        end
        scans_h2 = horzcat(scans_h{:});
    elseif scanDim==4 % for 4d files
        scans_h2 = matchfiles (fullfile(scandir, '4dniis', '/', matchstring ));
    end
    
    scanfiles.(par.subTask) = vertcat(scans_h2{:});
    
    %% scanfiles across subtasks
    clear sf_h sf_h2
    
    if scanDim==3 % for 3d files
        for i=1:length(par.scansSelect.(par.task).all);           
            scans_h{i} = matchfiles (fullfile(scandir, ['run' prepend(i) '/' matchstring])); %#ok<AGROW>
            ix = par.scansSelect.(par.task).all(i);
            scansByRun{ix} = char(scans_h{i}); %#ok<AGROW>
        end
        scans_h2 = horzcat(scans_h{:});
    elseif scanDim==4 % for 4d files
        scans_h2 = matchfiles (fullfile(scandir, '4dniis', '/', matchstring));
    end
    
    scanfiles.all = vertcat(scans_h2{:});
    clear sf_h sf_h2
end
end

