function Y=diff(varargin)
%DIFF (overloaded)

% Author Johan Löfberg
% $Id: diff.m,v 1.1 2006-08-10 18:00:19 joloef Exp $

X = varargin{1};
n = X.dim(1);
m = X.dim(2);

Y = X;
x_lmi_variables = X.lmi_variables;
lmi_variables = [];

% Fast diff of large matrices. Convert to optmized case
switch nargin
    case 1
        if (min([n m]) > 1)
            Y = diff(X',1,2)';
            return
        end

    case 2
        if varargin{2} == 1 & (min([n m]) > 1)
            Y = diff(X',1,2)';
            return
        end
    case 3
        if varargin{2}==1 & varargin{3} == 2 & (min([n m]) > 1)
            shift = [-speye(m-1) spalloc(m-1,1,0)] + [spalloc(m-1,1,0) speye(m-1)];
            %            shift = shift(1:m-1,1:m);
            shift = kron(shift,speye(n));
            Y.basis = shift*X.basis;
            Y.dim(1) = n;
            Y.dim(2) = m - 1;
            % Reset info about conic terms
            Y.conicinfo = [0 0];
            Y = clean(Y);
            return
        elseif varargin{2}==1 & varargin{3} == 1  & (min([n m]) > 1)
            Y = diff(X',1,2)';
            return
        end

    otherwise
        error('To many input arguments.')
end

% Slow but safe case. Used for higher order diff
try
    j = 1;
    diffX = diff(reshape(X.basis(:,1),n,m),varargin{2:end});
    Y.basis=diffX(:);
    for i = 1:length(x_lmi_variables)
        diffX = diff(reshape(X.basis(:,i+1),n,m),varargin{2:end});
        if (norm(diffX,inf)>0)
            Y.basis(:,j+1) = diffX(:);
            lmi_variables = [lmi_variables x_lmi_variables(i)];
            j = j+1;
        end
    end
    if j~=1
        Y.lmi_variables = lmi_variables;
    else
        Y = full(reshape(Y.basis(:,1),size(diffX,1),size(diffX,2)));
        return
    end
catch
    error(lasterr)
end
Y.dim(1) = size(diffX,1);
Y.dim(2) = size(diffX,2);
% Reset info about conic terms
Y.conicinfo = [0 0];