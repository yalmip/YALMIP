function Y=trapz(X,I,dummy)
%TRAPZ (overloaded)

if length(X.dim)==2
    if nargin == 1
        I = 1;
    end
    if min(X.dim) == 1 && nargin == 1
        Y = ones(1,max(X.dim))*reshape(X,[],1);
        return
    end
    switch I
        case 1
			Y = ones(1,X.dim(1)-1) * (X(1:m-1,:) + X(2:m,:))/2;
        case 2
			Y = (X(1:m-1,:) + X(2:m,:))/2 * ones(X.dim(1)-1,1);
        otherwise
            Y = X;
    end
    return
end
try
    n = X.dim(1);
    m = X.dim(2);
    if nargin==1
        Y = X;
        if n==1 || m==1
            % Spezialized code...
            Y.basis = trapz(Y.basis,1);
            Y.dim(1) = 1;
            Y.dim(2) = 1;
            temp = 1;
        else
            % Standard case

            temp = trapz(reshape(X.basis(:,1),n,m));
            Y.basis =  kron(speye(m),ones(1,n))*X.basis;

        end
    else
        Y = X;
        temp = trapz(reshape(X.basis(:,1),n,m),I);
        Y.basis = temp(:);
        for i = 1:length(Y.lmi_variables)
            temp = trapz(reshape(X.basis(:,i+1),n,m),I);
            Y.basis(:,i+1) = temp(:);
        end
    end
catch
    error(lasterr)
end
Y.dim(1) = size(temp,1);
Y.dim(2) = size(temp,2);
% Reset info about conic terms
Y.conicinfo = [0 0];
Y.extra.opname = '';
Y = flush(Y);
Y = clean(Y);
