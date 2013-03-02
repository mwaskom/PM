function PM_anatDicomsToNii(subpar)

%construct inplane images
par = subpar;
cd(par.inplaneDir);

dcmFiles = dir('M*.dcm');
P = char(dcmFiles.name);
dcmHdrs = spm_dicom_headers(P);
spm_dicom_convert(dcmHdrs, 'all', 'flat', 'nii');


if ~exist(par.anatdir);
    mkdir (par.anatdir);
end

syserror = system (['mv *.nii ' par.anatdir '/In001.nii']);
if syserror
    sprintf('error constructing inplanes');
end

%construct hires images
cd(par.hiresDir);

dcmFiles = dir('M*.dcm');
P = char(dcmFiles.name);
dcmHdrs = spm_dicom_headers(P);
spm_dicom_convert(dcmHdrs, 'all', 'flat', 'nii');


if ~exist(par.anatdir);
    mkdir (par.anatdir);
end

syserror = system (['mv *.nii ' par.anatdir '/V001.nii']);
if syserror
    sprintf('error constructing hires');
end