function par = PM_Params(subNo, task)
% function par = parParams(subject)
%sets up parameters for parmap batching.  Currently set-up to run once per
%subject...

%profile on

fprintf('\nEstablishing Parameters for Subject %g... \n', subNo);


par.subNo = subNo;
par.task = task;
par.FMCorrect = 1;
thisMachine = 'alan';
par.subTask = 'DM';
par.loadScans = 0;
par.usePremadePar = 1;
%% subject-specific stuff

switch par.subNo
    case 1
        par.str.perc = 'pm_031711';
        par.str.mnem = 'pm_052311_2';
        par.scansSelect.perc.DM = 1:6;
        par.scansSelect.perc.loc = 7:9;
        par.perc.FM2Funcs = [];
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;    
        par.mnem.FM2Funcs = [];
        par.fmCorrect = 0;
        par.mnem.numvols = [250 250 250 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176];
        par.goodSub = 1;
        par.flagIt = 0;
        %no fieldmaps collected.
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
        par.scansSelect.perc.DM = 1:10;
        par.scansSelect.perc.loc = 11:13;
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        %par.perc.numvols = [0 176 176 176 176 176 176 176 176 176 176 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176 176];
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
        
        %data will be thrown out, and thus not analyzed.  Too much extra
        %code would be necessary to write in order to analyze this subject,
        %who is a bad subject for multiple reasons anyway...
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
        par.str.perc = '';%'pm_052311'; This subject moves too much.
        par.str.mnem = '';%'pm_052011';
        par.scansSelect.perc.DM = 1:9;
        par.scansSelect.perc.loc = 10:12;
        par.perc.FM2Funcs = [];
        par.scansSelect.mnem.DM = [1 2 6];
        par.scansSelect.mnem.loc = [3:5];
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
        par.scansSelect.perc.loc = 10:12;
        par.perc.FM2Funcs = [];
        par.scansSelect.mnem.DM = 1:3;
        par.scansSelect.mnem.loc = 4:6;
        par.mnem.FM2Funcs = [1 1 2 2 3 3];
        par.perc.FM2Funcs = [1 1 2 2 2 3 3 3 4 4 4 5 5];
        par.mnem.numvols = [250 250 250 176 176 176];
        %par.perc.numvols = [176 176 176 176 176 176 176 176 176 0 176 176 176];
        par.perc.numvols = [176 176 176 176 176 176 176 176 176 176 176 176];
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
    otherwise
        error('unrecognized subject');
end

par.scansSelect.perc.all = [par.scansSelect.perc.DM par.scansSelect.perc.loc];
par.scansSelect.mnem.all = [par.scansSelect.mnem.DM par.scansSelect.mnem.loc];

par.scans_to_include = par.scansSelect.(par.task).(par.subTask);

par.substr = par.str.(par.task);  

subject = par.str.(task);
par.FM2Funcs = par.(task).FM2Funcs;
par.numvols = par.(task).numvols(par.scans_to_include);

%%
%----specify params-----:

% ---------Directory names----------
if strcmp(thisMachine, 'alan')
    par.exptdir = '/Users/alangordon/mounts/w5/alan/perceptMnemonic/';
elseif strcmp(thisMachine, 'jesse')
    par.exptdir = '/Users/Jesse/w5/biac3/wagner5/alan/perceptMnemonic/';
end

par.scriptsdir = fullfile(par.exptdir,'scripts');
par.fmridir = fullfile(par.exptdir, 'fmri_data');

par.subdir = fullfile(par.fmridir, subject);
par.anatdir = fullfile(par.subdir, 'anat');
par.funcdir = fullfile(par.subdir, 'rawNiis');
par.behavdir = fullfile(par.subdir, 'behav');
par.rawdir = fullfile(par.subdir, 'raw');
par.rawNiiDir = fullfile(par.subdir, 'raw/niiArchives');
par.niiDir = fullfile(par.subdir, 'rawNiis');
par.logdir = fullfile(par.subdir, 'logfiles');
par.artrepdir = fullfile(par.subdir, 'artRepairDM');
par.classmat = fullfile(par.exptdir, '/fmri_data/mvpa_files/trainLocTestPerc_3Subs.mat');
par.meanfuncdir = fullfile(par.funcdir, 'meanFuncs');
par.analysisdir = fullfile(par.subdir, ['analysis_percByEv_train1Sess']);
par.funcs3d = fullfile(par.subdir, 'funcs3d');
par.FourDNiiDir = fullfile(par.funcdir, '4dniis');

%% exceptions
if strcmp(par.substr, 'pm_042211');
     par.rawdir = fullfile(par.subdir, 'raw_copy');
end

%%
dHiRes = dir(fullfile(par.rawdir, '*FSPGR*'));
dInplane = dir(fullfile(par.rawdir, '*Inplane*'));
par.hiresDir = fullfile(par.rawdir, dHiRes.name); 
par.inplaneDir = fullfile(par.rawdir, dInplane.name); 
%can add more later



%----------Scan Params----------------
%par.scans_to_include = [1 2 11];


par.numscans = sum(par.numvols); %number of scans
par.art.motThresh = .3; % motion artefact threshold
par.art.sigThresh = 3; % signal intensity artefact threshold
par.TR = 2; %TR in seconds

par.numslice = 36;
par.TA = par.TR-par.TR/par.numslice; %TA
par.sliceorder = [36:-2:2 35:-2:1]; %slice order (interleaved, dorsal to ventral);
par.refslice = floor(par.numslice/2); %reference slice (assumes middle)

%par.maxvol = [126 600 600 600 600]; %highest volume number FOR EACH RUN
par.dropvol = 6; %dropped volumes (Assumed at start of scan)
par.minvol = par.dropvol+1; %first non-dropped volume
%par.numvols = par.maxvol-par.dropvol; %number of volumes per scan
par.slicetiming(1) = par.TA / (par.numslice -1);%timing var for slice timing
par.slicetiming(2) = par.TR - par.TA; %timing var for slice timing

% anat info
par.inplaneimg = [par.anatdir '/In001.nii'];
par.hiresimg = [par.anatdir '/V001.nii'];



% variables for realigning/reslicing

% realign flags...
par.realflag.quality = 0.9;
par.realflag.fwhm = 5;
par.realflag.sep = 4;
par.realflag.rtm = 1;
% par.realflag.PW;  %if field exists, then weighting is done...
par.realflag.interp = 2;

% reslice flags
par.reslflag.mask = 1;
par.reslflag.mean = 1;
par.reslflag.interp = 4;
par.reslflag.which = 0;

% re-orienting info
par.reorientmat = [105.0000 105.0000 59.4000 3.1416 0 -1.5708 1.0000 1.0000 1.0000 0 0 0]; %the proper rotation and translation to get the functional image preprocessed by kendrick's code to match the anatomical
dRN = dir([par.funcdir '/R*.nii']);
rawNiis_h = {};
for i = 1:length(dRN), rawNiis_h{i}=fullfile(par.funcdir, dRN(i).name); end %where the preprocessed .niis are stored
par.rawNiis = char(rawNiis_h{:});

% coregistration info
par.cor_meanfunc = [par.funcdir '/scan01/meanascan01_0007.nii']; %can change this to make more robust!
par.cor_inimg = [par.anatdir '/In001.nii'];
par.cor_hiresimg = [par.anatdir '/V001.nii'];



% segment info
par.img2bSeg = par.cor_hiresimg;
par.segopts = struct('biascor',1,'GM',[0 0 1],'WM',[0 0 1],'CSF',[0 0 0],'cleanup',0);
par.segbiasfwhm = 60; % 60 is default in gui, 75 is default in command line for reasons unexplained
% see spm_config_preproc and spm_preproc(_write) for details


% normalization:
% gray matter:
par.graytemp = '/Applications/MATLAB_R2007b/spm5/apriori/grey.nii';
par.grsegs(1,:) = [par.anatdir '/c1' 'V001.nii']; 
par.grsegs(2,:) = [par.anatdir '/c2' 'V001.nii']; 
par.graysrcimg = [par.anatdir '/c1' 'V001.nii']; 
par.graywrimg = [par.anatdir '/c1' 'V001.nii']; 
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

% spgm normalization:
par.spgrnrmflags.preserve = 0;
par.spgrnrmflags.bb = [[-78 -112 -50];[78 76 85]];
par.spgrnrmflags.vox = [2 2 2];
par.spgrnrmflags.interp     = 1;
par.spgrnrmflags.wrap       = [0 0 0];

% smoothing funcs
par.s4mm_smoothkernel = [4 4 4];
par.s8mm_smoothkernel = [8 8 8];
par.s6mm_smoothkernel = [6 6 6];

% specmaskvars

% specmaskvars
par.specwrflags.preserve = 0;
par.specwrflags.bb = [[-78 -112 -50];[78 76 85]];
par.specwrflags.vox        = [1 1 1];
par.specwrflags.interp     = 1;
par.specwrflags.wrap       = [0 0 0];

par.specsmooth = [8 8 8];%changed from [20 20 20] on 071409;
par.maskimg = [par.anatdir '/mask.nii']; 
par.smaskimg = [par.anatdir '/smask.nii'];
par.tsmaskimg = [par.anatdir '/tsmask.nii'];
par.wmaskimg = [par.anatdir '/wmask.nii'];
par.swmaskimg = [par.anatdir '/swmask.nii'];
par.tswmaskimg = [par.anatdir '/tswmask.nii'];
par.twmaskimg = [par.anatdir '/twmask.nii'];
par.addsegs = 'i1 + i2';
par.maskthresh = 'i1 > .2';


% model vars...
par.timing.fmri_t0 = par.refslice;  %micro-time resolution stuff (changed from 8)
par.timing.fmri_t = par.numslice; %used to be 16-- changed based on conversation with melina
par.timing.units = 'secs';
par.bases.hrf.derivs = [0 0]; % Melina says no cost to doing time and dispersion derivs ([0 0] = no derivs) NOTE that this will impact how you make contrasts!!
par.volt = 1;%if you don't want to model volterra interactions (and you don't want to model volterra interactions) leave this at 1.
% par.sess.scans is specified after populating list names...
par.sess.multi = {fullfile(par.behavdir, 'ons.mat')};
par.sess.multi_reg = {fullfile(par.behavdir, 'regs.mat')};
par.sess.hpf = 128;  % has anyone played around with this AND linear regs?
par.cvi = 'AR(1)'; %note that this actually gets changed to AR(0.2) in spm_fmri_spm_ui.  
% It looks like you might be able to give it a custom value
% by simply putting a number in here, but I haven't played around with it
par.global = 'None';
% explicit mask
par.mask = {};%{par.tswmaskimg};

% contrast vars
par.constat = 'T';


%par.srascanfiles = findScans(par, 'sRrun*.mat');
%par.rascanfiles = findScans(par, 'Rrun*.nii');
%par.scanfiles = findScans(par, 'run*.nii');
%par.valid = fullfile(par.funcdir, 'valid.nii');
%par.mean = fullfile(par.funcdir, 'mean.nii');
%par.rascanfilesPlusValidAndMean = findScans(par, 'R*.nii');
%par.wrascanfilesPlusValidAndMean = findScans(par, 'wR*.nii');
% par.scanfiles = findScans(par, 'run*.nii');
% par.rascanfiles = findScans(par, 'Rrun*.nii');
% par.wrascanfiles = findScans(par, 'wRrun*.nii');

if ~isempty(par.substr) && par.loadScans
    %[par.rascanfiles par.rascanfilesByRun] = findScans(par, 'Rrun*.nii','preprocessed');
    %[par.swrascanfiles par.swrascanfilesByRun] = findScans(par, 'swRrun*.nii','preprocessed');
    [par.scanfiles par.scanfilesByRun]= findScans(par, 'scan*.nii','nonpreprocessed');
    [par.ascanfiles par.ascanfilesByRun]= findScans(par, 'ascan*.nii','nonpreprocessed');
    [par.wascanfiles par.wascanfilesByRun]= findScans(par, 'wascan*.nii','nonpreprocessed');
    [par.swascanfiles par.swascanfilesByRun]= findScans(par, 'swascan*.nii','nonpreprocessed');
else
    par.rascanfiles = [];
    par.rascanfilesByRun = [];
    par.swrascanfiles = [];
    par.swrascanfilesByRun = [];
end

end

% function scanfiles = findScans(par, matchstring)
% 
% %include something that gets rid of files that start with '.'
% %also, change names from srascanfilenames to filenames.
% scanfiles.loc = [];
% scanfiles.DM = [];
% scanfiles.all = [];
% 
% if ~isempty(par.str.(par.task));
%     allsrascanfiles_h = dir(fullfile(par.funcdir, matchstring));
%     allsrascanfilenames = {allsrascanfiles_h.name};
%     
%     if ~isempty(allsrascanfilenames)
%         allsrascanfiles.loc_h = allsrascanfilenames(par.scansSelect.(par.task).loc);
%         allsrascanfiles.DM_h = allsrascanfilenames(par.scansSelect.(par.task).DM);
%         
%         for i = 1:length(allsrascanfiles.loc_h)
%             allsrascanfiles.loc_h2{i} = fullfile(par.funcdir, allsrascanfiles.loc_h{i});
%         end
%         
%         for i = 1:length(allsrascanfiles.DM_h)
%             allsrascanfiles.DM_h2{i} = fullfile(par.funcdir, allsrascanfiles.DM_h{i});
%         end
%         
%         for i = 1:length(allsrascanfilenames)
%             allsrascanfiles.all{i} = fullfile(par.funcdir, allsrascanfilenames{i});
%         end
%         
%         scanfiles.loc = char(allsrascanfiles.loc_h2{:});
%         scanfiles.DM = char(allsrascanfiles.DM_h2{:});
%         scanfiles.all = char(allsrascanfiles.all{:});
%     end
% else
%     
% end
% 
% end

%%
function [scanfiles scansByVec] = findScans(par, matchstring, prepflag)
scanfiles = [];
scansByVec = [];

if ~isempty(par.scansSelect.(par.task).(par.subTask))
    
    for i=par.scansSelect.(par.task).(par.subTask);
        if strcmp(prepflag, 'preprocessed')
            sf_h{i} = matchfiles (fullfile(par.funcdir, ['run' prepend(num2str(i)) '/' matchstring ]));
        elseif strcmp(prepflag, 'nonpreprocessed')
            sf_h{i} = matchfiles (fullfile(par.niiDir, ['scan' prepend(num2str(i)) '/' matchstring ]));
        end
    end
    sf_h2 = horzcat(sf_h{:});
    
    
end

scanfiles.(par.subTask) = vertcat(sf_h2{:});

clear sf_h sf_h2

for i=par.scansSelect.(par.task).all;
    if strcmp(prepflag, 'preprocessed')
        sf_h{i} = matchfiles (fullfile(par.funcdir, ['run' prepend(num2str(i)) '/' matchstring]));
    elseif strcmp(prepflag, 'nonpreprocessed')
        sf_h{i} = matchfiles (fullfile(par.niiDir, ['scan' prepend(num2str(i)) '/' matchstring]));
    end
    scansByVec{i} = char(sf_h{i});
end
sf_h2 = horzcat(sf_h{:});


scanfiles.all = vertcat(sf_h2{:});

end


