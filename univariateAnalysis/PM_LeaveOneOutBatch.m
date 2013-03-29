function SA = PM_LeaveOneOutBatch(sub)

[sa SA] = PM_SA;
% thisSA = SA.perc.sa16_all;
% task = 'perc';
% 
% i=0;
% 
% par = PM_Params(sub, task, 1);
%         
% theseSubs = setdiff(thisSA, par.substr);
% gpar = PM_GroupParams(theseSubs, par);
% PM_GroupWholeshebang(gpar, 'mec');

PM_conjoinWithinSubject(sub)