

for i = 1:length(sa)
    par = PM_Params(sa(i), 'mnem');
    
    cd(par.analysisdir);
    
    
    delete('spmT*')
    delete('con*')
    
    load('SPM');
    SPM.xCon = [];
    
    save SPM SPM
    
end


