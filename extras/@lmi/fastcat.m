function X = fastcat(varargin)

X = varargin{1};
X = flatten(X);
nTOT = length(X.clauses);
for i = 2:nargin
    X.clauses{i} = varargin{i}.clauses{1};   
    nTOT = nTOT + length(varargin{i}.clauses);
    X.LMIid = [X.LMIid varargin{i}.LMIid];
end

% VERY FAST UNIQUE BECAUSE THIS IS CALLED A LOT OF TIMES....
i = sort(X.LMIid);
i = i(diff([i NaN])~=0);
if length(i)<nTOT
    [i,j] = unique(X.LMIid);
    X = subsref(X,struct('type','()','subs',{{j}}));
end
