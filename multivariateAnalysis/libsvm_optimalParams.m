
function [scratchpad] = train_liblinear(trainpats,traintargs,in_args,cv_args)
% 
%   
% Function written by Alan Gordon to fit mvpa toolbox conventions
 
v = sparse(double(trainpats'));

if (~isempty(strfind(in_args.libsvm, '-s 3')) || ~isempty(strfind(in_args.libsvm, '-s 4')))
    choice = traintargs';
else
%     if size(traintargs,1)==2
%         choice = 2*(traintargs(1,:)') - 1; %must use -1 and 1 labels
%     else
        for i=1:size(traintargs,2)
            choice(i,1) = find(traintargs(:,i))';
        end
%     end
end

voxel_num = size(v,2);
lambda_beta = 0;

trainOpts_orig = in_args.libsvm ;
%%  pick the optimal cost parameter
%MUST MAKE THIS GENERAL FOR >2 CLASSES
p1Set = [1000];
p2Set = [0 .1 1 10 100 1000 10000 10^5];
p3Set = [ 1e-5 1e-4 1e-3 1e-2 1e-1 1 10];

for p1 = 1:length(p1Set)
    for p2 = 1:length(p2Set)
        for p3 = 1:length(p3Set)
            
            
            thisChoice = choice;
            
            thisV = v;
            
            
            trainOptsOptimize = sprintf('%s -c %g -v 10 -r %g -g %g', trainOpts_orig, p1Set(p1), p2Set(p2), p3Set(p3));
              
            m = svmtrain(thisChoice, thisV, trainOptsOptimize);
            
            scratchpad.paramSearch.p1(p1).p2(p2).p3(p3) = m;
        end
    end
end





% train it, with leave one out cross validation
model = svmtrain(choice, v, trainOptsOptimize);

% liblinear naturally picks the first choice direction as the first
% choice.  This differs from mvpa toolbox, which keeps the '1' labels as
% the first choice regardless of what the first presented trial label is.
% This variable ensures consistency across these toolboxes.

choice_set = unique(choice);
for i = 1:length(choice_set)
    thisChoice = choice_set(i);
    findChoice = find(choice==thisChoice);
    earliestInstanceOfChoice(i) = findChoice(1);
end



%scratchpad.logreg.betas = [model.w(end) model.w(1:(end-1))];

% if size(unique(choice),1)==2
%     scratchpad.liblinOrientation = choice(1);
%     
%     if scratchpad.liblinOrientation==1
%         scratchpad.logreg.betas(:,1) = [model.w(end) model.w(1:(end-1))];
%     else
%         scratchpad.logreg.betas(:,1) = [-model.w(end) -model.w(1:(end-1))];
%     end
% else
%     scratchpad.logreg.betas = model.w;
% end

scratchpad.constant = in_args.constant;
scratchpad.choice = choice;
scratchpad.model = model;
scratchpad.args = in_args;
