function X = plus(X,Y)
%PLUS Merges two LMI objects to one LMI

% Author Johan Löfberg
% $Id: plus.m,v 1.10 2006-08-10 13:15:45 joloef Exp $

if isa(X,'constraint')
    X = set(X);
elseif isa(X,'sdpvar')
    X = set(X);
end

if isa(Y,'constraint')
    Y = set(Y);
elseif isa(Y,'sdpvar')
    Y = set(Y);
end

% Support set+ []
if isempty(X)
    X = Y;
    return
elseif isempty(Y)
    X = X;
    return
end

if ~((isa(X,'lmi')) & (isa(Y,'lmi')))
    error('Both arguments must be SET objects')
end

nX = length(X.LMIid);
nY = length(Y.LMIid);
if nX==0
    X = Y;
    return
end
if nY == 0
    return;
end

for i = 1:length(Y.clauses)
    X.clauses{end+1} = Y.clauses{i};
end

X.LMIid = [X.LMIid Y.LMIid];

% VERY FAST UNIQUE BECAUSE THIS IS CALLED A LOT OF TIMES....
i = sort(X.LMIid);
i = i(diff([i NaN])~=0); 
if length(i)<nX+nY    
    [i,j] = unique(X.LMIid);
    X = subsref(X,struct('type','()','subs',{{j}}));
end
