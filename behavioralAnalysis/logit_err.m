function err = logit_err(q)
% logit_err(q) returns the negative of the log likelihood of
%    obtaining Data in the global arrays given a logistic model with
%    parameters in vector Q.
%	The data is in columns such that Data(:,1) is ones (you must specify
%	it though.  Data(:,[2:m-1]) are the values for the other observed
%	variables. and Data(:,m) is a 1 or a 0.
%
global Data Plothandle
% fake data for debugging
%$$$ Data = [1 0 .5 0;
%$$$     1 2 .5 0
%$$$     1 2 .5 0
%$$$     1 2 .5 0
%$$$     1 2 .5 1
%$$$     1 2 1 1
%$$$     1 2 1 1
%$$$     1 2 1 0
%$$$     1 4 1 1]
%$$$ q = [0 1 2]

[m n] = size(Data);
qcol = q(:);				% we need a column vector of params
z = exp(Data(:,1:n-1) * qcol);		% mxn * nx1 = mx1 column vector
p1 = z ./ (1 + z);			% calculate the prob of getting a 1
% for each data entry. 
% Now calculate the joint probability of obtaining the data set conceived as
% a list of Bernoulli trials.  This is just p1 for trials = 1 and 1-p1 for
% trials of 0.

p = [p1(Data(:,n)==1); 1-p1(Data(:,n)==0)];
err = -sum(log(p));			% minus log likelihood
