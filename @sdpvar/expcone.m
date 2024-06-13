function y = expcone(x)
%EXPCONE Defines a low-level exponential cone constraint x(2)*exp(x(1)/x(2)) <= x(3)
%
% Input
%    x       : Linear 3x1 SDPVAR object
%
% Example
%
% Standard  exponential cone constraint x(2)*exp(x(1)/x(2)) <= x(3)
%    F = expcone(x)
%
% To quickly define several cones, the argument can be a matrix, and the
% command is then short-hand for 
% for i = 1:size(x,2);F = [F,expcone(x(:,i))];end 
%
% See also  @SDPVAR/CONE, @SDPVAR/PCONE


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