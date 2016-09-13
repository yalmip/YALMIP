function varargout = max_internal(varargin)

switch class(varargin{1})
    case 'double'
        varargout{1} = max(varargin{:});
    case 'char'
        extstruct.var = varargin{2};
        extstruct.arg = {varargin{3:end}};
        [F,properties,arguments]=max_model([],varargin{1},[],extstruct);
        varargout{1} = F;
        varargout{2} = properties;
        varargout{3} = arguments;
    otherwise
end

function [F,properties,arguments]=max_model(X,method,options,extstruct)
switch method
    case 'graph'
        
        t = extstruct.var;
        X = extstruct.arg{1};
        basis = getbase(X);
        inf_row = find(basis(:,1) == -inf);
        if length(inf_row)>0
            X(inf_row) = [];
        end
        F = t-X>= 0;
        arguments = X(:);
        properties = struct('convexity','convex','monotonicity','increasing','definiteness','none');
        
    case 'exact'
        arguments = [];
        F = ([]);
        t = extstruct.var;
        X = extstruct.arg{1};
        basis = getbase(X);
        inf_row = find(basis(:,1) == -inf);
        if length(inf_row)>0
            X(inf_row) = [];
        end
        X = reshape(X,length(X),1);
        if prod(size(X)) == 1
            F = F + (X == t);
        elseif (prod(size(X)) == 2) & ((nnz(basis(1,:))==0) | (nnz(basis(2,:))==0))
            % Special case to test a particular problem max(0,y), so we
            % keep it since it is optimized
            if (nnz(basis(2,:))==0)
                X = [0 1;1 0]*X;
            end
            [M,m] = derivebounds(X);
            d = binvar(1,1);         
            F = [F, 0 <= t <= M(2)*d, X(2)<=M(2)*d];
            F = [F, -(1-d)*M(2) <= t-X(2) <= (-m(2))*(1-d),X(2)>=m(2)*(1-d)];
         elseif all(ismember(getvariables(X),yalmip('binvariables'))) & (is(X,'lpcone') | is(X,'sdpcone'))
            % Special case max(x) where x is simple binary
            F = [F, X <= t, sum(X) >= t];
        else
            F = max_integer_model(X,t);
        end
        arguments = [arguments;X(:)];
        properties = struct('convexity','convex','monotonicity','increasing','definiteness','none','model','integer');
        
    otherwise
        F = [];
        return
end