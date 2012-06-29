function sys = and(X,Y)
%AND Overloaded
%   
%   See also   LMI

% Author Johan Löfberg
% $Id: and.m,v 1.3 2005-02-04 10:10:26 johanl Exp $

% TODO : Check if binaries etc
if isa(X,'sdpvar')
    X = true(X);
end
if isa(Y,'sdpvar')
    Y = true(Y);
end

if isa(X,'constraint')
    X = set(X);
end
if isa(Y,'constraint')
    Y = set(Y);
end

nX = size(X.clauses,2); % Number of clauses in X
nY = size(Y.clauses,2); % Number of clauses in Y
sys = X;
% Add the two objects (append Y to X)
for i = 1:nY;    
    sys.clauses{i+nX} = Y.clauses{i};    
end
sys.LMIid = [X.LMIid Y.LMIid];