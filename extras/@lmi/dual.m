function sys = dual(X)
%DUAL Extract dual variable
%   
%   Z = DUAL(F)     Returns the dual variable for the constraint F
% 
%   See also SOLVESDP
  
% Author Johan Löfberg
% $Id: dual.m,v 1.4 2009-04-29 12:44:40 joloef Exp $

nlmi = size(X.clauses,2);

% Is it an empty SET
if (nlmi == 0) 
    sys = [];
    return
end

if nlmi>1
    if ~all(is(X,'elementwise'))
        error('Dual not applicable on list of constraints. Use dual(F(index)) or dual(F(''tag''))')
    end
end

% Get dual from repospitory
sys = [];
for i = 1:length(X.LMIid)
    sys = [sys;yalmip('dual',X.LMIid(i))];
end

% If no dual available, returns NaNs with correct dimension
if isempty(sys)
    sys = real(double(X))+NaN;
end