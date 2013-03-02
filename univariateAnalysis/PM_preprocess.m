% This is a script that calls preprocessfmri_standard.m.
% Edit the top section of this script to fit your needs and then run it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT THIS:
function PM_preprocess(par)

% start parallel MATLAB to speed up execution.
if matlabpool('size')==0
  matlabpool open;
end

subdir = par.subdir;
% what directory do the data live in?
datadir = par.rawdir;

% where should i save figures to?
figuredir = fullfile(subdir, 'figures2');

if par.FMCorrect
    % where are the raw fieldmap files to use?
    % (to omit the fieldmap-based undistortion process, just set fieldmapfiles to [].)
    fieldmapfiles  = matchfiles([datadir '/*spiral*/P*.7']); % oldest first
    %matchfiles([datadir '/*spiral*/P*.7']); % oldest first
    %matchfiles([datadir '/raw_spirals/P*.7']); % oldest first
else
    fieldmapfiles = [];
end

% what are the time values to associate with the fieldmaps?
% if [], default to 1:N where N is the number of fieldmaps.
fieldmaptimes = [];

% what rotation and flipping should we apply to the fieldmaps?
% [0 1 0] should be a set of values that should work in all circumstances.
% however, we may find as we collect new data that this is not true.  
% it is important to always check the result figures!
% [] means to use an interactive figure window to obtain the values to use.
% if you do pass in [], there must be at least one in-plane run supplied, 
% since the interactive figure window uses the in-plane run to compare against.
fieldmaporient = [0 1 0];

% what is the voxel size for the fieldmaps?
fieldmapsizes = [210/256 210/256 3.3];

% what is the difference in TE (in milliseconds) for the two volumes in the fieldmaps?
% (hint: after entering in the value of map_deltaf, check the value of map_delta in 
% the CV vars of the spiral fieldmap sequence.)
fieldmapdeltate = 2.272;

% should we attempt to unwrap the fieldmaps? (note that 1 defaults to a fast, 2D-based strategy; 
% see preprocessfmri.m for details.)  if accuracy is really important to you and the 2D strategy 
% does not produce good results, consider switching to a full 3D strategy like 
% fieldmapunwrap = '-f -t 0' (however, execution time may be very long).
fieldmapunwrap = 1;

% how much smoothing (in millimeters) along each dimension should we use for the fieldmaps?
% the optimal amount will depend on what part of the brain you care about.
% I have found that 7.5 mm may be a good general setting.
fieldmapsmoothing = [7.5 7.5 7.5];

% what DICOM directories should we interpret as in-plane runs?
% it is okay if you also match things that are not DICOM directories; we'll just ignore them.
inplanefilenames = matchfiles([datadir '/*Inplane*']);

% what DICOM directories should we interpret as EPI runs?
% it is okay if you also match things that are not DICOM directories; we'll just ignore them.
epifilenames = matchfiles([datadir '/*EPItrig*']);
%matchfiles([datadir '/*EPItrig*']);
%matchfiles([datadir '/EPI*']);

% this input is important only if you acquired something other than just the magnitude images
% for the EPI data.  the format is documented in dicomloaddir.m (the input <phasemode>).
% for the purposes of this script, the second element of <epiphasemode> must have exactly one
% element (and that element should probably be either 1 or 5).  if you acquired just the 
% magnitude images, just leave epiphasemode set to [].
epiphasemode = [];

% what is the desired in-plane matrix size for the EPI data?
% this is useful for downsampling your data (in order to save memory) 
% in the case that the data were reconstructed at too high a resolution.  
% for example, if your original in-plane matrix size was 70 x 70, the 
% images might be reconstructed at 128 x 128, in which case you could 
% pass in [70 70].  what we do is to immediately downsample each slice
% using lanczos3 interpolation.  if [] or not supplied, we do nothing special.
epidesiredinplanesize = [];

% what is the TR in seconds for the EPI runs?
epitr = 2;

% what is the slice order for the EPI runs?
episliceorder = 'interleaved';

% what is the phase-encode direction for the EPI runs? (see preprocessfmri.m for details.)
% (note that 'InPlanePhaseEncodingDirection' in the dicominfo of the EPI files gives either COL 
% which means the phase-encode direction is up-down in the images (which is 1 or -1 in our
% convention) or ROW which means the direction is left-right in the images (which is 2 or -2 in
% our convention).  the problem is that we don't know if the phase direction has been flipped,
% which you can do manually via CV vars.  it appears that if you don't explicitly flip, you should
% use the positive version (i.e. 1 or 2), meaning that the direction is towards the bottom
% (in the case of 1) or to the right (in the case of 2).  but in any case you should always
% check the sanity of the results!
epiphasedir = 1;

% what is the total readout time in milliseconds for an EPI slice?
% (note that 'Private_0043_102c' in the dicominfo of the EPI files gives the time per phase-encode line in microseconds.
% I confirmed that this time is correct by checking against the waveforms displayed by plotter.)
epireadouttime = 516/1000*64/2;  % divide by 2 if you are using 2x acceleration

% what fieldmap should be used for each EPI run? ([] indicates default behavior, which is to attempt
% to match fieldmaps to EPI runs 1-to-1, or if there is only one fieldmap, apply that fieldmap
% to all EPI runs, or if there is one more fieldmap than EPI runs, interpolate each successive
% pair of fieldmaps; see preprocessfmri.m for details.)
epifieldmapasst = par.FM2Funcs;

% how many volumes should we ignore at the beginning of each EPI run?
numepiignore = 6;

% by default, we tend to use double format for computation.  but if memory is an issue,
% you can try setting <dformat> to 'single', and this may reduce memory usage.
dformat = 'single';

% what volume should we use as reference in motion correction? ([] indicates default behavior which is
% to use the first volume of the first run; see preprocessfmri.m for details.  set to NaN if you
% want to omit motion correction.)
motionreference = [];

% what cut-off frequency should we use for filtering motion parameter estimates? ([] indicates default behavior
% which is to low-pass filter at 1/90 Hz; see preprocessfmri.m for details.)
motioncutoff = [Inf];

% what extra transformation should we use in the final resampling step? ([] indicates do not perform an extra transformation.)
extratrans = [];

% what is the desired resolution for the resampled volumes? ([] indicates to just use the original EPI resolution.)
targetres = [];

% should we perform slice shifting?  if so, specify band-pass filtering cutoffs in Hz, like [1/360 1/20].
% probably should be left as [] which means to do nothing special.
sliceshiftband = [];

% these are constants that are used in fmriquality.m.  it is probably 
% fine to leave this as [], which means to use default values.
fmriqualityparams = [];

% what kind of time interpolation should we use on the fieldmaps (if applicable)?
% ([] indicates to use the default, which is cubic interpolation.)
fieldmaptimeinterp = [];

% should we use a binary 3D ellipse mask in the motion parameter estimation?
% if [], do nothing special (i.e. do not use a mask).
% if {}, then we will prompt the user to interactively determine the
%   3D ellipse mask (see defineellipse3d.m for details).  upon completion,
%   the parameters will be reported to the command window so that you can
%   simply supply those parameters if you run again (so as to avoid user interaction).
% if {MN SD}, then these will be the parameters that determine the mask to be used.
mcmask = [];

% how should we handle voxels that have NaN values after preprocessing?
% if [], we use the default behavior which is to zero out all voxels that have a NaN
% value at any point in the EPI data.  see preprocessfmri.m for other options.
maskoutnans = [];

% savefile:  what .nii files (accepting a 1-indexed integer) should we save the final EPI data to?
%            (in the special EPI flattening case, we save the data to raw binary files (time x voxels) instead of .nii files.)
% savefileB: what .nii file should we save the valid voxels (binary mask) to?  ([] means do not save.)
% savefileC: what .nii file should we save the mean volume to?  ([] means do not save.)
% (we automatically make parent directories if necessary.)
savefile = [datadir '/run%02d.nii'];
savefileB = [datadir '/valid.nii'];
savefileC = [datadir '/mean.nii'];

% what .txt file should we keep a diary in?
diaryfile = [datadir '/diary.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW:

  mkdirquiet(stripfile(diaryfile));
  diary(diaryfile);
preprocessfmri_standard;
  diary off;
