function x = true(x)
% TRUE Constrains a variable to be positive
%
% TRUE(x) returns the constraint set(x>=1).
%
% For safety, it is advised to use the TRUE operator when working with
% logic constraints, instead of relying on the automatic constraints used
% by YALMIP (expression generated using AND and OR are automatically
% assumed to be constrained to be true.
%
% SET(a|b) is automatically changed to SET(TRUE(a|b)), or equivalently
% SET((a|b) >= 1), while SET(a) will be interpreted as SET(a>=0). To
% constrain a to be true, the user has to explicitely use SET(TRUE(a)).
%
%   See also SDPVAR/FALSE, SDPVAR/AND, SDPVAR/OR, SDPVAR/NOT, BINVAR, BINARY

x = set(x>=1);