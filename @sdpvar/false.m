function x = false(x)
% FALSE Constrains a variable to be negative
%
% FALSE(x) returns the constraint set(x<=0).
%
% Example
%
% SET(FALSE(a))    constrains binary variable a to be 0 (i.e. false)
% SET(FALSE(a|b))  constrains binary variables a and b to be 0 (i.e. false)
%
%   See also SDPVAR/TRUE, SDPVAR/AND, SDPVAR/OR, SDPVAR/NOT, BINVAR, BINARY

x = set(x<=0);