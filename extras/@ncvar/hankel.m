function H = hankel(c,r)
%HANKEL (overloaded)

% Author Johan Löfberg
% $Id: hankel.m,v 1.1 2006-08-10 18:00:20 joloef Exp $

% direct 1-to-1 copy of MATLAB double code
c = reshape(c,prod(size(c)),1);
nc = length(c);
if nargin < 2,
    r = zeros(size(c));
end
r = reshape(r,prod(size(r)),1);
nr = length(r);
x = [c;extsubsref(r,2:nr)];
cidx = (1:nc)';
ridx = 0:(nr-1);
H = cidx(:,ones(nr,1)) + ridx(ones(nc,1),:);
H = extsubsref(x,H);
if isa(H,'sdpvar')
    H.conicinfo = [0 0];
end

