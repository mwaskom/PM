function PM_anatDicoms2Analyze(subpar)
%
% convert dicom functional images into analyze images, move them into
% appropriate directories, and name them something sensible
% 
% alan gordon, 3/21/11


origdir = pwd;
% ---load par params if need be---
if isstruct(subpar) % if it is par_params struct
    par = subpar;
else % assume subject string
    par = par_params(subpar);
end

% which directories contain functional dicoms of interest?
InPlane = dir(fullfile(par.rawdir, '*Inplane*'));
SPGR = dir(fullfile(par.rawdir, '*SPGR*'));


%% Inplane
cd (fullfile(par.rawdir, InPlane.name));

pdcms = dir('*.dcm');
P = vertcat(pdcms.name);

%read dicom headers
display(sprintf ('reading dicom headers for Inplane'));
dcmhdrs = spm_dicom_headers(P);

%convert dicoms to analyze format, write into current dir
display(sprintf ('converting dicoms for Inplane'));
spm_dicom_convert(dcmhdrs);

par.thisanatdir = fullfile(par.anatdir);

if ~exist(par.thisanatdir)
    mkdir(par.thisanatdir);
end
    
system(['mv *.hdr ' par.thisanatdir]);
system(['mv *.img ' par.thisanatdir]);    

cd(par.thisanatdir);

dIn = dir('sA*.img');

v = spm_vol(dIn.name);
vol = spm_read_vols(v);
v.fname = 'In001.img';

spm_write_vol(v, vol);
system('rm sA*');


%% SPGR
cd (fullfile(par.rawdir, SPGR.name));

pdcms = dir('*.dcm');
P = vertcat(pdcms.name);

%read dicom headers
display(sprintf ('reading dicom headers for SPGR'));
dcmhdrs = spm_dicom_headers(P);

%convert dicoms to analyze format, write into current dir
display(sprintf ('converting dicoms for SPGR'));
spm_dicom_convert(dcmhdrs);

par.thisanatdir = fullfile(par.anatdir);

if ~exist(par.thisanatdir)
    mkdir(par.thisanatdir);
end
    
system(['mv *.hdr ' par.thisanatdir]);
system(['mv *.img ' par.thisanatdir]);    



cd(par.thisanatdir);

dV = dir('sA*.img');

v = spm_vol(dV.name);
vol = spm_read_vols(v);
v.fname = 'V001.img';

spm_write_vol(v, vol);

system('rm sA*');

