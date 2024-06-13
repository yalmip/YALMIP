function X=cumsum(X,I)
%CUMSUM (overloaded)

if nargin == 1 & min(X.dim)==1
    X.basis = cumsum(X.basis);     
else
    if nargin == 1
        I = min(find(X.dim>1));
        if isempty(I)
            I = 1;
        end
    end    
    if I == 2
        X = (tril(ones(size(X,2)))*X')';
    else
        X = tril(ones(size(X,1)))*X;
    end
end
X.conicinfo = [0 0];
X.extra.opname = '';
X = clean(X);