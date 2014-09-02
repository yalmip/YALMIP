function varargout = solvemp(P,x,y,options)
%SOLVEMP Computes solution to multi-parametric optimization problem
%
% min_z(x)   h(x,z)
% subject to
%            F(x,z) > 0
%
%
% [SOL, DIAGNOSTIC,Z,HPWF,ZPWF] = SOLVEMP(P,options,x,y)
%
% SOL        : Multi-parametric solution (see MPT toolbox)
%
% DIAGNOSTIC : struct with diagnostic information
%
% Z          : SDPVAR object with the detected decision variable z
%
% HPWF       : The value function as a pwf function
%
% ZPWF       : The optimal decision variable as a pfw function
%
% Input
%    P        : Optimization model
%    options  : solver options. See SDPSETTINGS. Can be [].
%    x        : Parametric variables
%    y        : Requested decision variables (subset of z)
%
% NOTE : If you are solving a problem leading to an mpMILP, the
% output SOL will be a set-valued map. To obtain the minimal
% solution (without so called overlaps), run removeOverlaps(SOL). If you
% have requested the 5th output ZPWF, overlaps are automatically removed.
% If your problem leads to an mpMIQP,  the output SOL will also be a
% set-valued map, but there is currently no way in MPT to obtain a
% non-overlapping solution. To use the solution in MPT, the command
% mpt_mergeCS(SOL) can be useful. Notice that the fifth output argument
% not will be available for mpMIQP problems.
%
% See also PARAMETRIC, SDPSETTINGS, YALMIPERROR

if nargin < 4
    options = sdpsettings;
else
    options = P.Options;
end

[sol,info] = solvemp(P.Constraints,P.Objective,options,x,y)
