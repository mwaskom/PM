function PM_MakeRegsPerc_parModByConfAndRT(par)


cd (par.behavdir);

dFN = dir('*Perc*');
fileNames = {dFN.name}';

countdown = 12; %length of countdown period (in secs) at the beginning of each run 
               
[res psy idx] = Perceptual_fMRIBehAnalysis(par);

%perf_set = {'cor' 'inc'};
class_set = {'face' 'house'};
resp_set = {'resp_face' 'resp_house'};

i = 0;

for l=1:length(class_set)
    %for p=1:length(perf_set)
        
        %idx.thisOns = idx.(perf_set{p}) .* idx.(resp_set{l}) .* ~isnan(idx.rt);
        idx.thisOns = (idx.(class_set{l}) | idx.(resp_set{l}) .* idx.corWithZeros ) .* ~isnan(idx.rt);
        
        if ~isempty(idx.alltrials(find(idx.thisOns)))
            
            i = i+1;
            %fName = sprintf('%s_%s', class_set{l}, perf_set{p});
            fName = sprintf('%s_%s', class_set{l});
            
            stimOnsets{i}= idx.alltrials(find(idx.thisOns));
            stimNames{i} = fName;
            stimDurations{i} = 0;
            
            thisConf = idx.conf_unsigned(find(idx.thisOns));
            thisRT = idx.rt(find(idx.thisOns));  
            
            pmod(i).param{1} = thisRT;
            pmod(i).param{2} = thisConf;
            
            pmod(i).name{1} = 'RT';
            pmod(i).name{2} = 'Conf';
            
            pmod(i).poly{1} = 1;
            pmod(i).poly{2} = 1;
            
        end
    %end
end

% if sum(idx.coh_==0 & idx.resp_face)>0
%     i = i+1;
%     
%     idx.thisOns = idx.coh_==0 & idx.resp_face & ~isnan(idx.rt);
%     
%     stimOnsets{i} = idx.alltrials(find(idx.thisOns));
%     stimNames{i} = 'zeroCohRespFace';
%     stimDurations{i} = 0;
%     
%     thisConf = idx.conf_unsigned(find(idx.thisOns));
%     thisRT = idx.rt(find(idx.thisOns));
%     
%     pmod(i).param{1} = thisRT;
%     pmod(i).name{1} = 'RT';
%     pmod(i).poly{1} = 1;
%     
%     if (sum(idx.thisOns)>2 && var(thisConf)>0)
%         pmod(i).param{2} = thisConf;
%         pmod(i).name{2} = 'Conf';
%         pmod(i).poly{2} = 1;
%         
%     end
% end
% 
% if sum(idx.coh_==0 & idx.resp_house)>0
%     i = i+1;
%     
%     idx.thisOns = idx.coh_==0 & idx.resp_house & ~isnan(idx.rt);
%     
%     stimOnsets{i} = idx.alltrials(find(idx.thisOns));
%     stimNames{i} = 'zeroCohRespHouse';
%     stimDurations{i} = 0;
%     
%     
%     thisConf = idx.conf_unsigned(find(idx.thisOns));
%     thisRT = idx.rt(find(idx.thisOns));
%     
%     pmod(i).param{1} = thisRT;
%     pmod(i).name{1} = 'RT';
%     pmod(i).poly{1} = 1;
%     
%     if (sum(idx.thisOns)>2 && var(thisConf)>0)
%         pmod(i).param{2} = thisConf;
%         pmod(i).name{2} = 'Conf';
%         pmod(i).poly{2} = 1;
%     end
% end

if sum(idx.junk>0)
    i = i+1;
    stimOnsets{i} = idx.alltrials(find(idx.junk + isnan(idx.rt)));
    stimNames{i} = 'junk';
    stimDurations{i} = 0;
end

onsets = stimOnsets;
names = stimNames;
durations = stimDurations;


if ~exist(par.analysisdir)
    mkdir (par.analysisdir);
end



sessReg = zeros(sum(par.numvols),length(par.numvols)-1);
for i = 1:(length(par.numvols)-1)
    if (i==1)
        sessReg(1:par.numvols(i),i) = ones(par.numvols(i),1);
    else
        sessIdx = (1+sum(par.numvols(1:(i-1)))):sum(par.numvols(1:i));
        sessReg(sessIdx,i) = ones(par.numvols(i),1);
    end
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

save ons.mat onsets durations names pmod;
save regs.mat R