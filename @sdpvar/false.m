function x = false(x)
% FALSE Constrains a variable to be negative
%
% TRUE(x) returns the constraint (x<=0.5). It is assumed that x is
% binary. The reason (x<=0.5) is used instead of x==0 is that some big-M
% modelling turns out to be less sensitive to numerical issues in this
% form. Once the model enters the solver, it will be trivially presolved
% anyway.
%
% Example
%
% (FALSE(a))    constrains binary variable a to be 0 (i.e. false)
% (FALSE(a|b))  constrains binary variables a and b to be 0 (i.e. false)
%
%   See also SDPVAR/TRUE, SDPVAR/AND, SDPVAR/OR, SDPVAR/NOT, BINVAR, BINARY

x = (x<=.5);