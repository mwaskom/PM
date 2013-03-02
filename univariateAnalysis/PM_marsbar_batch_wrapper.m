function SA = PM_marsbar_batch_wrapper(sa)

if nargin<1
    sa = 1:22;
end

tasks = {'mnem'};

i=0;
for C = sa
    for t =1:length(tasks)
        par = PM_Params(C, tasks{t}, 1);
        
        if par.goodSub
            i = i+1;
            %PM_PatternwiseRegression(C, tasks{t})
            PM_wholeshebang(par, 't')
            SA{i} = par.substr;

        end

    end
end
