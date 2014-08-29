function [SOL, DIAGNOSTIC,Z,HPWF,ZPWF] = solvemp(P,Domain)
%SOLVEMP Computes solution to multi-parametric optimization problem
%
% min_z(x)   h(x,z)
% subject to
%            F(x,z) >= 0
%
%
% [SOL, DIAGNOSTIC,Z,HPWF,ZPWF] = SOLVEMP(P)
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
%    P        : Optimizer object
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
% See also PARAMETRIC, SET, SDPSETTINGS, YALMIPERROR

%P.F = [P.F,0 <= P.input.expression <= 100];
[SOL, DIAGNOSTIC,Z,HPWF,ZPWF] = solvemp(P.F,P.h,[],P.input.expression,P.output.expression);
