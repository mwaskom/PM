function SA = PM_batchshebang(sa)

if nargin<1
    sa = 1:22;
end

tasks = {'perc'};

i=0;
for C = sa
    for t =1:length(tasks)
        par = PM_Params(C, tasks{t}, 1);
        
        if par.goodSub
            %PM_subUtil(par);
            PM_wholeshebang(par, 't')
        end
    end
end


