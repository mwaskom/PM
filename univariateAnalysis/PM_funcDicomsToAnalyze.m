function PM_funcDicomsToAnalyze(subpar)
%transform dicom images into analyze.  This is pretty slow...

par = subpar;

cd(par.subdir);

% make required dirs...
if ~exist(par.funcdir,'dir')
    mkdir(par.funcdir)
end
    

EPIDirs = dir('*EPItrig*');

for i = 1:length(EPIDirs)
    
    cd (fullfile(par.subdir, EPIDirs(i).name))
    dcmFiles = dir(fullfile(par.subdir, EPIDirs(i).name, '*.dcm'));
    
    
    P = char(dcmFiles.name);
    
    dcmHdrs = spm_dicom_headers(P);
    
    
    spm_dicom_convert(dcmHdrs);
    
    thisScanDir = fullfile(par.funcdir, ['scan' prepend(num2str(i)) ]);
    if ~exist(thisScanDir );
        mkdir (thisScanDir);
    end
    
    movefile('*.img', thisScanDir);
    movefile('*.hdr', thisScanDir);
    
    cd (thisScanDir)
    
    d = dir(fullfile(thisScanDir, '*.hdr'));
    
    for s = 1:length(d)
        v = spm_vol(d(s).name);
        vol = spm_read_vols(v);
        
        %v.mat = newVolMat.mat;
        JFormatted = prepend(num2str(s),3);
        v.fname = ['scan' prepend(num2str(i)) '.V'  JFormatted '.img'];
        
        spm_write_vol(v, vol);
        
    end
end

