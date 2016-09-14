function [P, SOL, DIAGNOSTIC] = solvemp(P)
%SOLVEMP Computes solution to parametric optimization problem
%
% min_z(x)   h(x,z)
% subject to
%            F(x,z) >= 0
%
%
% [P, SOL, DIAGNOSTIC] = SOLVEMP(P)
%
% P          : Optimizer object with parametric solution attached
%
% SOL        : Multi-parametric solution (see MPT toolbox)
%
% DIAGNOSTIC : struct with diagnostic information
%
% Input
%    P        : Optimizer object
% Examle
%  sdpvar x y u1 u2
%  Model = [-1 <= [u1 u2] <= 1, -5 <= [x,y] <= 5];
%  Objective = (x-u1)^2+(y-u2)^2;
%  P = optimizer(Model,Objective,[],[x y],[u1 u2])
%  P = solvemp(P)
%  P{[1 .3]}
%
% See also OPTIMIZER

if ~isempty(P.output.z)
    error('Parametric solution can only be computed on optimizer object with simple output');
end
if length(P.dimoutOrig)>1
    error('Parametric solution can only be computed on optimizer object with as single matrix output');
end

[SOL, DIAGNOSTIC] = solvemp(P.F,P.h,P.model.options,P.input.expression,P.output.expression);
if ~isempty(SOL) && DIAGNOSTIC.problem == 0
    P.ParametricSolution = SOL;
end