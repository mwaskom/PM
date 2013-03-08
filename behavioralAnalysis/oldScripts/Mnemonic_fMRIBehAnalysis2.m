function [res M] = Mnemonic_fMRIBehAnalysis()

%Fix RTs to include early responses
% in general, clean this up!

exptDir = '/Users/alangordon/mounts/wagner5/alan/perceptMnemonic/fmri_data';
subs = {'pm_012611'};


hands = {{'1' '2' '3' '4'} {'5' '6' '7' '8'}};

analyzeCats = {'face' 'house'};


for s=1:length(subs)
    
    fName  = dir(fullfile(exptDir, subs{s}, 'behav', 'AG4*retrieve*'));
    
    %dat = load(fullfile(exptDir, subs{s}, 'behav', fName.name));
    
    
    
    reh.wordList = {};
    reh.stimresp = {};
    reh.judgeResp = {};
    reh.judgeRT = {};
    reh.stimRT = {};
    reh.cond = [];
    for r = 1:3
        
        
        thisReh = dir(['*retrieve' '*(' num2str(r) ')*']);
        
        
        RL(r).dat = load(thisReh.name);
        
        
        
        reh.wordList = vertcat(reh.wordList, RL(r).dat.list.words);
        reh.stimResp = horzcat(reh.stimresp, RL(r).dat.theData.stimresp);
        reh.judgeResp = horzcat(reh.judgeResp, RL(r).dat.theData.judgeResp);
        reh.judgeRT = horzcat(reh.judgeRT, RL(r).dat.theData.judgeRT);
        reh.stimRT = horzcat(reh.judgeRT, RL(r).dat.theData.stimRT);
        reh.cond = horzcat(reh.cond, RL(r).dat.theData.cond);
        
    end
    
    respVec = reh.judgeResp;
    RTVec = reh.judgeRT;
    idx.earlyResp = find(~strcmp(reh.stimresp, 'noanswer'));
    
    respVec(idx.earlyResp) = reh.stimresp(idx.earlyResp);
    RTVec(idx.earlyResp) = reh.stimRT(idx.earlyResp);
    
    %takes the first recorded response
    %also, changes '1!' into '1'
    reh.resp = arrayfun(@(x) x{1}(1), respVec, 'UniformOutput',false);
    reh.RT = arrayfun(@(x) x{1}, RTVec);
    
    idx.face = reh.cond==1;
    idx.house = reh.cond==2;
    
    HN = RL(1).dat.S.retHandNum;
    

    
    
    idx.respFace = ismember(reh.resp, hands{HN});
    idx.respHouse =  ismember(reh.resp, hands{3-HN});
    
    h.idx.cor.Face = idx.respFace.*idx.face;
    h.idx.cor.House = idx.respHouse.*idx.house;
    
    idx.cor = idx.respFace.*idx.face + idx.respHouse.*idx.house;
    idx.inc = idx.respFace.*idx.house + idx.respHouse.*idx.face;
    
    idx.cleanResp = (idx.respFace + idx.respHouse) ;
    
    
    
    for a=1:length(analyzeCats)
            
            res.pctCor.(analyzeCats{a}) = sum(idx.(analyzeCats{a}) .* idx.cor .*idx.cleanResp)/sum(idx.(analyzeCats{a}) .*idx.cleanResp);
            
            res.pctHouseResps = sum(idx.respHouse .* idx.cleanResp)/sum(idx.cleanResp);
                      
            res.pctCleanAnswers.(analyzeCats{a}) = mean(idx.cleanResp);
            
            
%             M.pctCor(s,a) = res.sub(s).pctCor.(analyzeCats{a}).(thisCoh);
%             M.pctCorOnlyRehCor(s,a) = res.sub(s).pctCorOnlyRehCorTrials.(analyzeCats{a}).(thisCoh);
%             M.pctHouseResps(s) = res.sub(s).pctHouseResps.(thisCoh);
%             M.medRT(s,a) = nanmedian(RT(find(idx.(analyzeCats{a}) .* idx.Coh .* idx.inc .*idx.cleanNonLures)));
%             
%             M.rehHits(s,a) = res.sub(s).reh.Hits.(analyzeCats{a});
%             M.rehCRs(s,a) = res.sub(s).reh.CRs.(analyzeCats{a});
%             M.HitsMinusFAs(s,a) = res.sub(s).reh.Hits.(analyzeCats{a}) - (1 - res.sub(s).reh.CRs.(analyzeCats{a})) ;
    end 
end