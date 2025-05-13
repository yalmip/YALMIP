function y = expcone(varargin)
%EXPCONE Defines a low-level exponential cone constraint
%
% Input
%    x       : Linear 3x1 or 3xn SDPVAR object
% alternatively
%    x,y,z   : Linear scalar SDPVAR objects
% Example
%
% Standard  exponential cone constraint x(2)*exp(x(1)/x(2)) <= x(3)
%    F = expcone(x)
% Alternative  y*exp(x/y)<=z
%    F = expcone(x,y,z)
%
% To quickly define several cones, the argument can be a matrix, and the
% command is then short-hand for 
% for i = 1:size(x,2);F = [F,expcone(x(:,i))];end 
%
% See also  @SDPVAR/CONE, @SDPVAR/PCONE

if nargin > 1 && numel(varargin{1}) >  1
    error('Vectorized format not allowed with multiple arguments')
end
if nargin == 1 || nargin == 3
    x = [varargin{:}];       
else
    error('EXPCONE expects either 1 or 3 arguments');
end

[n,m] = size(x);
if min([n m])==1
    x = reshape(x,3,1);  
end
[n,m] = size(x);
if n ~=3
    error('x must be a vector or matrix of height 3')
end
y = x;
if min([n m])>1
	y.typeflag = 22;
else
	y.typeflag = 21;
end
y = lmi(y);