function S = spy(X)
%SPY (overloaded)

% Author Johan Löfberg 
% $Id: spy.m,v 1.3 2005-11-25 14:15:14 joloef Exp $   

if length(X)>1
    S = [];
    for i = 1:length(X)
        S = blkdiag(S,spy(X.clauses{1}.data));
    end    
    if nargout == 0
        spy(S);
    end
else
    X = spy(sdpvar(X));
    if nargout==0
        spy(X);
    else
        S = X;
    end
end