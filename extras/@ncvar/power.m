function y = power(x,d)
%POWER (overloaded)

% Sanity check
if prod(size(x))==1 & (prod(size(d))>1)
    x = x.*ones(size(d));
end
if prod(size(d))>1
    if any(size(d)~=size(x))
        error('Matrix dimensions must agree.');
    end
else
    d = ones(x.dim(1),x.dim(2))*d;
end

if any(d ~= fix(d)) | any(d<0)
    error('Only non-negative integer powers allowed in non-commuting variables');
end

% Just loop and call mpower instead
y = [];
dimx = size(x);
x = x(:);
d = d(:);
for i = 1:length(x)
    y = [y;extsubsref(x,i,1)^d(i)];
end
y = reshape(y,dimx);
