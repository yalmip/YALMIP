function Y=sum(X,I,dummy)
%SUM (overloaded)

if length(X.dim)==2
    if nargin == 1
        I = 1;
    end
    if min(X.dim) == 1 & nargin == 1
        Y = ones(1,max(X.dim))*reshape(X,[],1);
        return
    end
    if isequal(sort(I), [1 2]) || isequal(I, 'all')
        % special case if we sum all elements
        Y = ones(1,X.dim(1))*X*ones(X.dim(2),1);
        return
    end
    switch I
        case 1
            Y = ones(1,X.dim(1))*X;
        case 2
            Y = X*ones(X.dim(2),1);
        otherwise
            if isa(I,'double') && (I>=1) && (I == ceil(I))
                Y = X;
            else
                error('Dimension argument must be a positive integer scalar within indexing range.')
            end
    end
    return
end
try
    n = X.dim(1);
    m = X.dim(2);
    if nargin==1
        Y = X;
        if n==1 | m==1
            % Spezialized code...
            Y.basis = sum(Y.basis,1);
            Y.dim(1) = 1;
            Y.dim(2) = 1;
            temp = 1;
        else
            % Standard case

            temp = sum(reshape(X.basis(:,1),n,m));
            Y.basis =  kron(speye(m),ones(1,n))*X.basis;

        end
    else
        Y = X;
        temp = sum(reshape(X.basis(:,1),n,m),I);
        Y.basis = temp(:);
        for i = 1:length(Y.lmi_variables)
            temp = sum(reshape(X.basis(:,i+1),n,m),I);
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
