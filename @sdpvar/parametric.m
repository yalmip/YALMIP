function x = parametric(x)
%PARAMETRIC Defines a variable as parametric
%
%   F = PARAMETRIC(x) is used to describe the set of parametric variables
%   in a multi-parametric program, as an alternative to using the 4th input
%   in SOLVEMP
%
%   It can also be used to define the parametric variables in a SOS program
%
%   INPUT
%    x : SDPVAR object
%
%   OUTPUT
%    F : Constraint
%
%   EXAMPLE
%    F = parametric(x);               % Short notation
%
%   See also SOLVEMP, SOLVESOS, SDPVAR

x.typeflag = 13;
x = lmi(x);