function t = toeplitz(c,r)
%TOEPLITZ (overloaded)

% direct 1-to-1 copy of MATLAB double code
if nargin < 2,
    c.basis(1,:) = conj(c.basis(1,:));
    r = c;
    c.basis = conj(c.basis);
    %  c(1) = conj(c(1));
    %  r = c;
    %  c = conj(c);
end
r = reshape(r,prod(size(r)),1);%r(:)
p = length(r);
m = length(c);
x = [extsubsref(r,p:-1:2) ; reshape(c,prod(size(c)),1)];
cidx = (0:m-1)';
ridx = p:-1:1;
t = cidx(:,ones(p,1)) + ridx(ones(m,1),:);
t = extsubsref(x,t);
if isa(t,'sdpvar')
    t.conicinfo = [0 0];
end

