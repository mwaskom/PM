function [res] = PM_concatenateClassifierResults(filePattern)
cwd = pwd;
cd /biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/mvpa_files
d = dir([filePattern '*']);
dF = {d.name};

if length(dF)>17
   error('too many files match the filePattern'); 
end

for i = 1:length(dF)
    thisRes = load(dF{i});
    res.subj{i} = thisRes.res.subj{1};
    res.subjArray(i) = thisRes.res.subjArray;
end

save (filePattern, 'res');

for i = 1:length(dF)
    delete(dF{i})
end

cd (cwd)
end

