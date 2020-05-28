function X = subsasgn(X,I,Y)
% SUBSASGN (overloaded)

if isempty(X)
    y = Y;
    sizeY = size(Y);
    y_lmi_variables = y.lmi_variables;    
    X0 = subsasgn([],I,full(reshape(full(Y.basis(:,1)),sizeY)));
    dim = size(X0);
    y.basis = reshape(X0,prod(dim),1);    
    for i = 1:length(y_lmi_variables)
        X0 = subsasgn([],I,full(reshape(full(Y.basis(:,i+1)),sizeY)));
        y.basis(:,i+1) = reshape(X0,prod(dim),1);
    end
    y.dim = dim;
    % Reset info about conic terms
    y.conicinfo = [0 0];
    y.basis = sparse(y.basis);
    if length(dim)>2 && ~isa(y,'ndsdpvar')
        y = ndsdpvar(y);
    end
    y = flush(y);    
    X = y;
    return
end
base = reshape(1:size(X.basis,1),X.dim);
base = subsref(base,I);

if isa(Y,'ndsdpvar')
    Y = sdpvar(Y);
elseif isa(Y,'double')
    Y = Y(:);
end

dim = X.dim;
X = sdpvar(X);
X(base) = Y;
X = reshape(X,dim);
X = clean(X);