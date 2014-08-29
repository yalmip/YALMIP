function x = false(x)
% FALSE Constrains a binary variable to be false (0)
%
% FALSE(x) returns the constraint x<=0.5.
%
%   See also SDPVAR/TRUE, SDPVAR/AND, SDPVAR/OR, SDPVAR/NOT, BINVAR, BINARY

x = (x<=0.5);