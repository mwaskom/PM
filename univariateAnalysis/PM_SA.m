function [sa SA] = PM_SA()

sa.sa_all = [1 2 3 6 7 8 10 11 12 14 15 16 18 19 20 21 22];
sa.sa16_all = [1 2 3 6 8 10 11 12 14 15 16 18 19 20 21 22];
sa.sa14_all = [1 3 6 8 10 11 12 14 16 18 19 20 21 22];

sa.sa16_CorVsInc = [1 2 3 6 8 10 11 12 14 15 16 18 20 21 22];
sa.sa16_Conf = [1 2 3 6 8 10 11 12 14 15 16 18 20 21];

sa.sa14_CorVsInc = [1 3 6 8 10 11 12 14 16 18 20 21 22];
sa.sa14_Conf = [1 3 6 8 10 11 12 14 16 18 20 21];

sa.sa16_mnemFaceConf1 = [1 3 6 8 10 11 12 14 15 16 18 21];
fn = fieldnames(sa);
tsk =  {'perc', 'mnem'};

for f=1:length(fn)
    for t = 1:length(tsk)
        thisStruct = sa.(fn{f});
        for i=1:length(thisStruct)
            thisPar = PM_Params(thisStruct(i), tsk{t}, 0);
            SA.(tsk{t}).(fn{f}){i} = thisPar.substr;
        end
    end
end

end