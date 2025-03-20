function y = simplex(x)
% SIMPLEX Adds a constraint that argument is member of (standard) simplex
%
% Input
%    x       : SDPVAR object
%
% Example
%
% Add constraint that sum(x)==1 and x>=0
%    F = simplex(x)

y = reshape(x,[],1);
y.typeflag = 70;
y = lmi(y);