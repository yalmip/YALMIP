function y = allvariables(varargin)
% ALLVARIABLES Return a vector with all variables used in expressions
%
% y = allvariables(A,B,...)
%
% Returns variables used in the objects A, B, ...
%
% The arguments can be SDPVARs or constraints
%
% Example
%
% x = sdpvar(3,1);y = sdpvar(1);
% Model = [-1 <= x <=1, y == sum(x)];
% Obj = sum(x) + y^2;
% Obj2 = Obj + norm(allvariables(Model,Obj))

y = depends(varargin{:});
y = recover(y);