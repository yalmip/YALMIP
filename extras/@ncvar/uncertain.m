function x = uncertain(x)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(x) is used to describe the set of parametric variables
%   in an uncertain program, as an alternative to using the 4th input in
%   SOLVEROBUST.
%
%   If an uncertain multi-parametric problem is solved, UNCERTAIN has to be
%   used to declare the set of uncertain variables
%
%   INPUT
%    x : SDPVAR object
%
%   OUTPUT
%    F : SET object
%
%   See also SOLVEROBUST, ROBUSTIFY

% Author Johan Löfberg
% $Id: uncertain.m,v 1.1 2006-08-10 18:00:23 joloef Exp $

x.typeflag = 15;
x = lmi(x);