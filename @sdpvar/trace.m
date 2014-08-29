function Y=trace(X)
%TRACE (overloaded)

if isa(X,'blkvar')
    X = sdpvar(X);
end

if X.typeflag~=0
	error('Relational objects cannot be manipulated')
end

X = flush(X);
Y = X;
x_lmi_variables = X.lmi_variables;
lmi_variables = [];
n = X.dim(1);
m = X.dim(2);

if n == m % Standard case fast
    ind = 1:(n+1):n^2;
    Y.basis = sum(X.basis(ind,:),1);
    ind = find(Y.basis(2:end));
    if isempty(ind)
        Y = full(Y.basis(1));
    else
        Y.basis = Y.basis([1 1+ind]);
        Y.lmi_variables = X.lmi_variables(ind);
        Y.dim(1) = 1;
        Y.dim(2) = 1;
    end
else
    traceX = trace(reshape(X.basis(:,1),n,m));
    Y.basis = traceX(:);

    j = 1;
    for i = 1:length(x_lmi_variables)
        traceX = trace(reshape(X.basis(:,i+1),n,m))/1;
        if (norm(traceX,inf)>0)
            Y.basis(:,j+1) = traceX(:);
            lmi_variables = [lmi_variables x_lmi_variables(i)];
            j = j+1;
        end
    end
    if isempty(lmi_variables)
        Y = full(trace(reshape(X.basis(:,1),n,m)));
    else
        Y.dim(1) = size(traceX,1);
        Y.dim(2) = size(traceX,2);
        Y.lmi_variables = lmi_variables;
        % Reset info about conic terms
        Y.conicinfo = [0 0];
        Y.extra.opname='';
    end
end
