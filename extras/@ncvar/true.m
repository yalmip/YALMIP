function x = true(x)
% TRUE Constrains a binary variable to be positive
%
% TRUE(x) returns the constraint (x>=.5).
%
% For safety, it is advised to use the TRUE operator when working with
% logic constraints, instead of relying on the automatic constraints used
% by YALMIP (expression generated using AND and OR are automatically
% assumed to be constrained to be true.
%
%   See also SDPVAR/FALSE, SDPVAR/AND, SDPVAR/OR, SDPVAR/NOT, BINVAR, BINARY

x = (x>=.5);