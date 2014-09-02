function [p,d] = check(X)
% CHECK(F)  Displays/calculates constraint residuals on constraint F
%
% [pres,dres] = CHECK(F)
%
% pres : Primal constraint residuals
% dres : Dual constraint residuals
%
% If no output argument is supplied, tabulated results are displayed
%
% Primal constraint residuals are calculated as:
%
%  Semidefinite constraint F(x)>0 : min(eig(F))
%  Element-wise constraint F(x)>0 : min(min(F))
%  Equality constraint F==0       : -max(max(abs(F)))
%  Second order cone t>||x||      : t-||x||
%  Integrality constraint on x    : max(abs(x-round(x)))
%  Sum-of-square constraint       : Minus value of largest (absolute value) coefficient 
%                                   in the polynomial p-v'*v
%
% Dual constraints are evaluated similarily.
%
%  See also  SOLVESDP, SOLVESOS, SOSD, DUAL

switch nargout
    case 0
        check(lmi(X));
    case 1
        p = check(lmi(X));
    case 2
        [p,d] = check(lmi(X));
end
