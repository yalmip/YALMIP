function [P,x] = polytope(v,x)
% polytope  Creates a constraint defining a polytope from vertices
%
% P = polytope(v,x)
%
% Short for P = [x == v*s, sum(s==1), s >= 0]
%
% Note that it is not projected to x but lives in (x,s)
%
% See also VERTEX, BOUNDINGBOX, CHEBYBALL, POLYTOPE

if ~isa(x,'sdpvar')
    error('First argument should be an SDPVAR');
end

if size(v,2)==length(x)
    warning('Dimension of appears flipped, so I transposed it');
    v = v';
elseif ~size(v,1)==length(x)
    error('Dimensions on x and vertices do not match');
end
try
    k = convhulln(v');
    v = v(:,unique(k(:)));
catch
end
x = reshape(x,[],1);
s = sdpvar(size(v,2),1);
P = [x == v*s, sum(s) == 1, s >= 0];