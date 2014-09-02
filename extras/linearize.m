function plin = linearize(p)
% LINEARIZE Linearize SDPVAR object
%
% h = LINEARIZE(p)
%
% Returns linearization p(double(x)) + dp(double(x))*(x-double(x))
% where x is the SDPVAR variables defining the polynomial p(x)
%
% See also SDPVAR, JACOBIAN

if isa(p,'double')
    plin = zeros(size(p));
    return
end

if is(p,'linear') & ~is(p,'compound')
    plin = p;
    return
end

x = recover(depends(p));
x0 = double(x);
p0 = double(p);

n = size(p,1);
m = size(p,2);

if min(n,m)>1
    plin = [];
    for i = 1:m
        plin = [plin p0(:,i)+double(jacobian(p(:,i),x))*(x-x0)];
    end
else
    plin = p0+double(jacobian(p,x))*(x-x0);
end