function y = kron(X,Y)
%KRON (overloaded)

if isequal(X,1)
    y = Y;
    return;
end

if isequal(Y,1);
    y = X;
end

if (isa(X,'sdpvar') & isa(Y,'sdpvar'))
    y = [];
    s.type = '()';
    for i = 1:size(X,1)
        this_row = [];       
        for j = 1:size(X,2);
            s.subs = {i,j};
            this_row = [this_row subsref(X,s)*Y];
        end
        y = [y;this_row];
    end    
    return
end

if isa(X,'sdpvar')
    lmi_variables = getvariables(X);
    nv = length(lmi_variables);
    y  = X;
    sparse_Y = sparse(Y);
    % This one used also for checking size
    temp = kron(getbasematrix(X,0),sparse_Y);
    if size(Y,2)==1
        % Special case
        % [kron([N1(:) N2(:)],vec(M))]=[vec(kron(N1,M)) vec(kron(N2,M))]
        y.basis = kron(X.basis,sparse_Y);
    else
        y.basis=temp(:);
        X_base = X.basis;
        for i = 1:nv
           % temp = kron(getbasematrix(X,lmi_variables(i)),sparse_Y);            
            temp = kron(reshape(X_base(:,i+1),X.dim),sparse_Y);
            y.basis(:,i+1) = temp(:);
        end;
    end
end

if isa(Y,'sdpvar')
    lmi_variables = getvariables(Y);
    nv = length(lmi_variables);
    y  = Y;
    sparse_X = sparse(X);
    % This one used also for checking size
    temp = kron(sparse_X,getbasematrix(Y,0));
    if size(X,1)==1
        % In this special case
        %[kron(vec(M),[N1(:) N2(:)])]==[vec(kron(M,N1)) vec(kron(M,N2))]        
        y.basis = kron(sparse_X(:),Y.basis);
    else
        y.basis = temp(:);
        y.basis = spalloc(length(temp(:)),nv+1,0);
        y.basis(:,1) = temp(:);
        Y_base = Y.basis;
        for i = 1:nv
            % temp = kron(sparse_X,getbasematrix(Y,lmi_variables(i)));
            temp = kron(sparse_X,reshape(Y_base(:,i+1),Y.dim));
            y.basis(:,i+1) = temp(:);
        end
    end
end
y.dim(1) = size(temp,1);
y.dim(2) = size(temp,2);
y = clean(y);
% Reset info about conic terms
if isa(y,'sdpvar')
    y.conicinfo = [0 0];
end