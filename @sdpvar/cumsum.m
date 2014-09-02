function X=cumsum(X,I)
%CUMSUM (overloaded)

if nargin == 1 & min(X.dim)==1
    B = cumsum(X.basis);
else
    if nargin == 1
        I = min(find(X.dim>1));
        if isempty(I)
            I = 1;
        end
    end
    
    B = [];
    for i = 1:length(X.lmi_variables)+1
        C = reshape(X.basis(:,i),X.dim);
        C = cumsum(C,I);
        B = [B C(:)];
    end
end
X.basis = B;
X.conicinfo = [0 0];
X.extra.opname = '';
X = flush(X);
X = clean(X);