function t = toeplitz(c,r)
%TOEPLITZ (overloaded)

% direct 1-to-1 copy of MATLAB double code
if nargin < 2,
    c.basis(1,:) = conj(c.basis(1,:));
    r = c;
    c.basis = conj(c.basis);
end
r = reshape(r,prod(size(r)),1);
p = length(r);
m = length(c);
x = [extsubsref(r,p:-1:2) ; reshape(c,prod(size(c)),1)];
cidx = (0:m-1)';
ridx = p:-1:1;
t = cidx(:,ones(p,1)) + ridx(ones(m,1),:);
if isrow(t)
    t = extsubsref(x,t).';
else
    t = extsubsref(x,t)
end
if isa(t,'sdpvar')
    t.conicinfo = [0 0];
end

