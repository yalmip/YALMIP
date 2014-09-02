function x = integer(x)
% integer (overloaded)

x = sdpvar(x);
x = binary(x);