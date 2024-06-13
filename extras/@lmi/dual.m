function sys = dual(X)
%DUAL Extract dual variable
%   
%   Z = DUAL(F)     Returns the dual variable for the constraint F
% 
%   See also SOLVESDP
  
X = flatten(X);
nlmi = length(X.LMIid);

% Is it an empty SET
if (nlmi == 0) 
    sys = [];
    return
end

if nlmi>1
    if ~(all(is(X,'elementwise')) || all(is(X,'equality')))
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