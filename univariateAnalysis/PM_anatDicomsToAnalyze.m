function PM_anatDicomsToAnalyze(subpar)
%transform dicom images into analyze.  This is pretty slow...

par = subpar;

cd(par.inplaneDir);

dcmFiles = dir('*.dcm');

P = char(dcmFiles.name);

dcmHdrs = spm_dicom_headers(P);

spm_dicom_convert(dcmHdrs, 'all', 'flat', 'nii');


if ~exist(par.anatdir);
    mkdir (par.anatdir);
end

system (['mv *.nii ' par.anatdir '/In001.nii']);



cd(par.hiresDir);

dcmFiles = dir('*.dcm');

P = char(dcmFiles.name);

dcmHdrs = spm_dicom_headers(P);

spm_dicom_convert(dcmHdrs, 'all', 'flat', 'nii');


if ~exist(par.anatdir);
    mkdir (par.anatdir);
end

system (['mv *.nii ' par.anatdir '/V001.nii']);