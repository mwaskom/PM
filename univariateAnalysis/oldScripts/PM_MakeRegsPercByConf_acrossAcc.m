function PM_MakeRegsPercByConf_acrossAcc(par)


cd (par.behavdir);

dFN = dir('*Perc*');
fileNames = {dFN.name}';

countdown = 12; %length of countdown period (in secs) at the beginning of each run 
         
         


[res psy idx] = Perceptual_fMRIBehAnalysis(par);

%perf_set = {'corWithZeros' 'inc'};
class_set = {'face' 'house'};
resp_set = {'resp_face' 'resp_house'};


i = 0;
for c=1:length(idx.conf_unsigned)
    for l=1:length(class_set)
        %for p=1:length(perf_set)
            
            %idx.thisOns = (idx.conf_unsigned==c) .* idx.(perf_set{p}) .* idx.(resp_set{l});
            idx.thisOns = (idx.conf_unsigned==c) .* idx.(resp_set{l});
            
            if ~isempty(idx.alltrials(find(idx.thisOns)))
                
                i = i+1;
                %fName = sprintf('%s_conf%s_%s', class_set{l}, num2str(c), perf_set{p});
                fName = sprintf('%s_conf_%s', class_set{l}, num2str(c));
                
                stimOnsets{i}= idx.alltrials(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = 0;
            end
        %end
    end
end


onsets = stimOnsets;
names = stimNames;
durations = stimDurations;

if ~exist(par.analysisdir)
    mkdir (par.analysisdir);
end



sessReg = zeros(sum(par.numvols),length(par.numvols) - 1);
for i = 1:(length(par.numvols) - 1)
    sessReg((i-1)*par.numvols(i)+1:i*par.numvols(i),i) = ones(par.numvols(i),1);
end

% for i = 1:par.numscans
%     cd(fullfile(par.subdir, 'functional', ['scan' (prepend(num2str(par.scans_to_include(i)),2))]   ));
%     motTxt = dir('rp*');
%     motRegs_h{i} = textread(motTxt.name);
% end
% 
% motRegs = vertcat(motRegs_h{:});

R = horzcat(sessReg);




cd (par.analysisdir);

save ons.mat onsets durations names;
save regs.mat R