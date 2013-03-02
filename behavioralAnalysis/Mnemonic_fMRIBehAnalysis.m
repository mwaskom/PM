function [res M] = Mnemonic_fMRIBehAnalysis()

%Fix RTs to include early responses
% in general, clean this up!

exptDir = '/Users/alangordon/mounts/wagner5/alan/perceptMnemonic/fmri_data';
subs = {'pm_012211b'};

%rawDatDir  = '/Users/alangordon/mounts/wagner5/alan/perceptMnemonic/mnemonicPilotDat/ver4/data';

rehNums = [22 24 26 28 30 32 34 36 38 ];

hands = {'q' 'p'};

analyzeCats = {'face' 'house'};


for s=1:length(subs)
    
    fName  = dir(fullfile(exptDir, subs{s}, 'behav', 'AG4*retrieve*'));
    
    dat = load(fullfile(exptDir, subs{s}, 'behav', fName.name));
    
    if ~isfield(dat, 'retData')
        dat.retData = dat.theData;
    end
    
    idx.face = dat.retData.cond==1;
    idx.house = dat.retData.cond==2;
    
    HN = dat.S.retHandNum;
    cohs = unique(dat.retData.coh);
    
    respVec = dat.retData.judgeResp;
    idx.earlyResp = find(~strcmp(dat.retData.stimresp, 'noanswer'));
    respVec(idx.earlyResp) = dat.retData.stimresp(idx.earlyResp);
    
    
    idx.respFace = cellfun(@(x) strcmp(x(1), hands{HN}), respVec);
    idx.respHouse = cellfun(@(x) strcmp(x(1),hands{3-HN}), respVec);
    
    RT = .5 + cellfun(@(x) x(1), dat.retData.judgeRT);

    RTEarly = cellfun(@(x) x(1), dat.retData.stimRT);
    
    RT(find(idx.earlyResp)) = RTEarly(find(idx.earlyResp));
    
    h.idx.cor.Face = idx.respFace.*idx.face;
    h.idx.cor.House = idx.respHouse.*idx.house;
    
    idx.cor = idx.respFace.*idx.face + idx.respHouse.*idx.house;
    idx.inc = idx.respFace.*idx.house + idx.respHouse.*idx.face;
    
    cohs = unique(dat.retData.coh);
    
    
    %rehLists = dir ([listDir '/' '*Reh*' prepend subs{s}(1) '*']);
    
    reh.wordList = {};
    reh.rehLures = [];
    reh.judgeResp = {};
    for r = 1:length(rehNums)
        
        if s< 5
        thisReh = dir(['*retrieve_' num2str(s) '*(' num2str(rehNums(r)) ')*']);
        else
        thisReh = dir(['*rehearse*sub' num2str(s) '*(' num2str(r + length(rehNums)) ')*']);    
        end
        
        RL(r) = load(thisReh.name);
        
        
        
        reh.wordList = vertcat(reh.wordList, RL(r).list.words);
        reh.rehLures = vertcat(reh.rehLures, RL(r).list.rehLures);
        reh.judgeResp = horzcat(reh.judgeResp, RL(r).theData.judgeResp);
        

    end
    
    idx.reh.respOld = strcmp(hands{dat.S.encHandNum}, reh.judgeResp);
    idx.reh.respNew = strcmp(hands{3-dat.S.encHandNum}, reh.judgeResp);
    idx.reh.cleanResp = idx.reh.respOld + idx.reh.respNew;
    
    for i = 1:length(reh.wordList)
        idx.reh2Test(i) = find(strcmp(dat.retData.item(i), reh.wordList));
    end
    
    %rehearsal through learning was correct both times
    idx.rehCor = strcmp(hands{dat.S.encHandNum}, reh.judgeResp);
    
    %the words were not paired with a lure img during rehearsal
    idx.rehNonLure = reh.rehLures(idx.reh2Test);
    
    
    idx.cleanResp = (idx.respFace + idx.respHouse) ;
    
    idx.cleanNonLures = idx.cleanResp.* idx.rehNonLure';
    
    
    
    for a=1:length(analyzeCats)
        for c=1:length(cohs)
            thisCoh = ['coh_' num2str(cohs(c))];
            
            idx.Coh = dat.retData.coh==cohs(c);
            
            
            res.sub(s).pctCor.(analyzeCats{a}).(thisCoh) = sum(idx.(analyzeCats{a}) .* idx.Coh .* idx.cor .*idx.cleanNonLures)/sum(idx.(analyzeCats{a}) .* idx.Coh .*idx.cleanNonLures);
            res.sub(s).pctCorOnlyRehCorTrials.(analyzeCats{a}).(thisCoh) = sum(idx.(analyzeCats{a}) .* idx.Coh .* idx.cor .* idx.rehCor .* idx.cleanNonLures)/sum(idx.(analyzeCats{a}) .* idx.Coh .* idx.rehCor .* idx.cleanNonLures);
            
            
            %
            res.sub(s).pctHouseResps.(thisCoh) = sum(idx.respHouse .* idx.Coh .* idx.cleanNonLures)/sum(idx.Coh .* idx.cleanNonLures);
            %
            res.sub(s).reh.Hits.(analyzeCats{a}) = sum(idx.(analyzeCats{a}) .* idx.reh.respOld .* idx.reh.cleanResp .* ~reh.rehLures') / sum(idx.(analyzeCats{a}) .* ~reh.rehLures' .* idx.reh.cleanResp);
            
            res.sub(s).reh.CRs.(analyzeCats{a}) = sum(idx.(analyzeCats{a}) .* idx.reh.respNew .* idx.reh.cleanResp .* reh.rehLures') / sum(idx.(analyzeCats{a}) .* reh.rehLures' .* idx.reh.cleanResp);
            
            res.sub(s).test.pctCleanAnswers.(analyzeCats{a}) = sum(idx.cleanNonLures)/ sum(idx.rehNonLure);
            res.sub(s).reh.pctCleanAnswers.(analyzeCats{a}) = mean(idx.reh.cleanResp);
            
            
            M.pctCor(s,a) = res.sub(s).pctCor.(analyzeCats{a}).(thisCoh);
            M.pctCorOnlyRehCor(s,a) = res.sub(s).pctCorOnlyRehCorTrials.(analyzeCats{a}).(thisCoh);
            M.pctHouseResps(s) = res.sub(s).pctHouseResps.(thisCoh);
            M.medRT(s,a) = nanmedian(RT(find(idx.(analyzeCats{a}) .* idx.Coh .* idx.inc .*idx.cleanNonLures)));
            
            M.rehHits(s,a) = res.sub(s).reh.Hits.(analyzeCats{a});
            M.rehCRs(s,a) = res.sub(s).reh.CRs.(analyzeCats{a});
            M.HitsMinusFAs(s,a) = res.sub(s).reh.Hits.(analyzeCats{a}) - (1 - res.sub(s).reh.CRs.(analyzeCats{a})) ;
        end
        
        res.sub(s).pctCor.(analyzeCats{a}).allCohs = sum(idx.(analyzeCats{a}) .* idx.cor .*idx.cleanResp )/sum(idx.(analyzeCats{a}) .* idx.cleanResp);
    end
    
    
end