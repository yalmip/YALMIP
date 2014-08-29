function Y=sparse(varargin)
%SPARSE (overloaded)

if nargin < 3
    error('At-least 3 arguments needed');
end

data = varargin{3};
ns = varargin{1};
ms = varargin{2};

if length(ns)~=length(ms) | length(ms)~=length(data)
    error('Length of first 3 arguments must be equal');
end

if min(size(data))>1
    error('Third argument should be a vector');
end

if nargin < 4
    n = max(ns);
else
    n = varargin{4};
end
if nargin < 5
    m = max(ms);
else
    m = varargin{5};
end

if any(ms>m)
    error('Dimension mismatch')
end

if any(ns>n)
    error('Dimension mismatch')
end

Y = data;
Y.dim(1) = n;
Y.dim(2) = m;
[i1,j1,s1] = find(data.basis);
ind = ns+(ms-1)*n;
Y.basis = sparse(ind(i1),j1,s1,n*m,1+length(Y.lmi_variables));
% Reset info about conic terms
Y.conicinfo = [0 0];