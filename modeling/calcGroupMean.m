
function [m mse] = calcGroupMean(v, g, groups, data_type)
% 
%   function [m mse] = calcGroupMean(v, g, groups, data_type)
% 
% calculates mean and se of 'v' after dividing it into separate groups    
% g is the grouping variable and has the same length as v
% groups defines divisions of g. if unspecified groups is set to unique(g) 
% data_type defines the type of variable: 'binary' or 'continuous' (default) 
% 
% example
% calculate probability correct
%   [p, p_se] = calcGroupMean(cor==1, coh, coh_set, 'binary')
% calculate rt
%   [rt, rt_se] = calcGroupMean(rt, coh, coh_set)

if nargin<3 || isempty(g),
    groups = unique(g);
end;

if nargin<4 || isempty(data_type),
    data_type = 'continuous';
end;

if ~iscell(groups),
    m = arrayfun(@(s) mean(v(g==s)), groups);
    switch data_type,
        case 'binary',
            mse = arrayfun(@(s) sqrt(m(groups==s).*(1-m(groups==s))./sum(g==s)), groups);
        case 'continuous',
            mse = arrayfun(@(s) std(v(g==s))/sqrt(sum(g==s)), groups);
    end;        
else
    m = cellfun(@(s) mean(v(ismember(g,s))), groups);
    switch data_type,
        case 'binary',
            findcell = @(c,pat) cellfun(@(x,s) isequal(x,s), c, repmat({pat},size(c)));
            mse = cellfun(@(s) sqrt(m(findcell(groups,s)).*(1-m(findcell(groups,s)))./sum(ismember(g,s))), groups);
        case 'continuous',
            mse = cellfun(@(s) std(v(ismember(g,s)))/sqrt(sum(ismember(g,s))), groups);
    end;        
end;
        
