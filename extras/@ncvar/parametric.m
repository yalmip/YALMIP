function x = parametric(x)
%PARAMETRIC Defines a variable as parametric
%
%   F = PARAMETRIC(x) is used to describe the set of parametric variables
%   in a multi-parametric program, as an alternative to using the 4th input
%   in SOLVEMP
%
%
%   INPUT
%    x : SDPVAR object
%
%   OUTPUT
%    F : Constraint object
%
%   EXAMPLE
%    F = parametric(x);             
%
%   See also SOLVEMP, SDPVAR

x.typeflag = 13;
x = lmi(x);