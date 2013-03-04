function SA = PM_batchshebang(sa)

if nargin<1
    sa = 1:22;
end

tasks = {'mnem'};

i=0;
for C = sa
    for t =1:length(tasks)
        par = PM_Params(C, tasks{t}, 1);
        
        if par.goodSub
            PM_wholeshebang(par, 't')
        end
    end
end


% %                 PM_subUtil(par);
%         %PM_conjoin(C);
%         %         if par.goodSub
%         %             this_sa = setdiff(sa, C);
%         %             gpar = PM_GroupParams(this_sa,par);
%         %             % %
%         %             PM_GroupWholeshebang(gpar, 'mec');
%         %         end
%         
%         
%                     %             SA{i} = par.substr;
%                     
%                                 %             i = i+1;
%             %             %PM_PatternwiseRegression(C, tasks{t})