function varargout = min_internal(varargin)

switch class(varargin{1})
    case 'double'
        varargout{1} = min(varargin{:});
    case 'char'
        extstruct.var = varargin{2};
        extstruct.arg = {varargin{3:end}};
        [F,properties,arguments]=min_model([],varargin{1},[],extstruct);
        varargout{1} = F;
        varargout{2} = properties;
        varargout{3} = arguments;
    otherwise
end

function [F,properties,arguments] = min_model(X,method,options,extstruct)
switch method
    case 'graph'
        arguments=[];
        F = ([]);
        basis = getbase(extstruct.arg{1});
        inf_row = find(basis(:,1) == inf);
        if length(inf_row)>0
            extstruct.arg{1}(inf_row) = [];
        end
        F = F + (extstruct.arg{1} - extstruct.var >= 0);
        arguments= [arguments;extstruct.arg{1}(:)];
        properties = struct('convexity','concave','monotonicity','increasing','definiteness','none');
        
    case 'exact'
        arguments = [];        
        t = extstruct.var;
        X = extstruct.arg{1};
        
        basis = getbase(X);
        inf_row = find(basis(:,1) == inf);
        if length(inf_row)>0
            X(inf_row) = [];
        end
        
        X = reshape(X,length(X),1);
                        
        if prod(size(X)) == 1
            F = [X == t];
        elseif all(ismember(getvariables(X),yalmip('binvariables'))) & (is(X,'lpcone') | is(X,'sdpcone'))
            % Special case min(x) where x is simple binary
            F = [X >= t, sum(X) <= t+length(X)-1];
        else
            F = max_integer_model(-X,-t);
        end
                
        arguments = [arguments;X(:)];
        properties = struct('convexity','concave','monotonicity','increasing','definiteness','none','model','integer');
        
    otherwise
        F = [];
        return
end