function [] = PM_make3dNiftis(par)


d = dir(fullfile(par.rawNiiDir, '*EPItrig*.nii.gz'));
fileNames = {d.name};

dotFiles = cell2mat(cellfun(@(x) x(1)=='.', fileNames, 'UniformOutput', false));
fileNames(find(dotFiles)) = []; %remove hidden files that are prepended with dots.

if (length(d) ~= length(par.numvols))
   warning('scan length differs from par'); 
end

for i = 1:length(fileNames)
    
    %rename EPITrig to 'scan##'
    cd (par.rawNiiDir);
    thisRawScan = ['scan' prepend(num2str(i)) '.nii.gz'];
    movefile(fileNames{i}, thisRawScan);
    
    %unzip it.
    unix(['gunzip ' thisRawScan]);
    thisUnzippedScan =  thisRawScan(1:end-3);
    
    % make a run-specific diretory
     scanDir = ['scan' prepend(num2str(i))];
     newDir = fullfile(par.niiDir, scanDir);
     mkdir(newDir);
%     
%     % move the 4d nii to the run-specific directory
     movefile(fullfile(par.rawNiiDir, thisUnzippedScan), newDir);
     cd (newDir);
%     
%     % 4d nii to 3d nii
     unix(['fslsplit ' thisUnzippedScan ' ' thisUnzippedScan(1:end-4) '_' ' -t']);
%     
%     % move raw scan back to its directory
     movefile(fullfile(par.niiDir, scanDir, thisUnzippedScan), par.rawNiiDir);
%     
%     % unzip the newly made niis
     unix('gunzip scan*.gz');

     
    % hide discarded volumes
    discardDir = fullfile(par.niiDir, scanDir, 'discardedVols');
    mkdir(discardDir);
    dS = dir('scan*.nii');

     dropVols = {dS(1:par.dropvol).name}';
    
    for j = 1:length(dropVols)
        movefile(dropVols{j}, discardDir);
    end
end