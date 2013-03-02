function [ ] = PM_unzip( par )
%unzip files that come from the CNI

if exist(par.rawdir)
    cd (par.rawdir);
    
    d = dir('*.zip');
    
    if ~isempty(d);
        unix ('unzip -qo \*.zip');
        
        if ~exist('dcmArchives')
            unix ('mkdir dcmArchives');
        end
        
        if ~exist('niiArchives')
            unix ('mkdir niiArchives');
        end
        
        
        
        unix ('mv *.zip dcmArchives');
        unix ('mv *.nii.gz niiArchives');
        unix ('mv *.nii niiArchives');
    end
else
    warning(['no valid raw dir for subject ' par.substr]);
end

