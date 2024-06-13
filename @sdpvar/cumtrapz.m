function X=cumtrapz(X,I)
%CUMTRAPZ (overloaded)

if nargin == 1 && min(X.dim)==1
    B = cumtrapz(X.basis);
else
    if nargin == 1
        I = find(X.dim>1,1);
        if isempty(I)
            I = 1;
        end
    end
    
    B = [];
    for i = 1:length(X.lmi_variables)+1
        C = reshape(X.basis(:,i),X.dim);
        C = cumtrapz(C,I);
        B = [B C(:)];
    end
end
X.basis = B;
X.conicinfo = [0 0];
X.extra.opname = '';
X = clean(X);