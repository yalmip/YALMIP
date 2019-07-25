function A=linspace(x1,x2,N)
%LINPSPACE (overloaded)

if nargin < 3
    N = 100;
elseif isa(N,'sdpvar')
    error('SDPVAR/LINSPACE does not support N as a decision variable');
end

A = x1 + linspace(0,1,N)*(x2-x1);