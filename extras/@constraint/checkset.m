function [p,d] = checkset(X)
% CHECKSET(F)  Displays/calculates constraint residuals on constraint F
%
% [pres,dres] = CHECKSET(F)
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
%  Rank constraint rank(X) < r     : r-rank(X)
%  Sum-of-square constraint       : Minus value of largest (absolute value) coefficient 
%                                   in the polynomial p-v'*v
%
% Dual constraints are evaluated similarily.
%
%    See also   SET, SOLVESDP, SOLVESOS, SOSD, DUAL

% Author Johan Löfberg
% $Id: checkset.m,v 1.19 2008-02-21 15:32:09 joloef Exp $

switch nargout
    case 0
        checkset(lmi(X));
    case 1
        p = checkset(lmi(X));
    case 2
        [p,d] = checkset(lmi(X));
end
