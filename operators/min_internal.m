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
        F = set([]);
        basis = getbase(extstruct.arg{1});
        inf_row = find(basis(:,1) == inf);
        if length(inf_row)>0
            extstruct.arg{1}(inf_row) = [];
        end
        F = F + set(extstruct.arg{1} - extstruct.var);
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
        
        
%         
%         if 1
%             F = max_integer_model(-X,-t);
%         else
%             [M,m] = derivebounds(X);
%             n = length(X);
%             d = binvar(n,1);
%             F = F + set(sum(d)==1);
%             F = F + set(-(max(M)-min(m))*(1-d) <= t-X <= (max(M)-min(m))*(1-d));
%             
%             kk = [];
%             ii = [];
%             for i = 1:n
%                 k = [1:1:i-1 i+1:1:n]';
%                 ii = [ii;repmat(i,n-1,1)];
%                 kk = [kk;k];
%                 Mm = M(k)-m(i);
%             end
%             xii = extsubsref(X,ii);
%             dii = extsubsref(d,ii);
%             xkk = extsubsref(X,kk);
%             F = F + set(xii <= xkk+(M(ii)-m(kk)).*(1-dii));
%         end
        arguments = [arguments;X(:)];
        properties = struct('convexity','concave','monotonicity','increasing','definiteness','none','model','integer');
        
    otherwise
        F = [];
        return
end