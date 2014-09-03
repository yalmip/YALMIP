function x = true(x)
% TRUE Constrains a variable to be positive
%
% TRUE(x) returns the constraint (x>=0.5). It is assumed that x is
% binary. The reason (x>=0.5) is used instead of x==1 is that some big-M
% modelling turns out to be less sensitive to numerical issues in this
% form. Once the model enters the solver, it will be trivially presolved
% anyway.
%
% For safety, it is advised to use the TRUE operator when working with
% logic constraints, instead of relying on the automatic constraints used
% by YALMIP (expression generated using AND and OR are automatically
% assumed to be constrained to be true.
%
% (a|b) is automatically changed to (TRUE(a|b)). To constrain a to be true,
% the user has to explicitely use (TRUE(a)). 
%
%   See also SDPVAR/FALSE, SDPVAR/AND, SDPVAR/OR, SDPVAR/NOT, BINVAR, BINARY

x = (x>=.5);