function varargout = vectorizedlogic(fun,varargin);
if nargin-1 == 1
    if length(varargin{1})==1
       % What should be returned on operator(scalar)?
       varargout{1} = varargin{1};
    else
        x = varargin{1};        
        temp = fun(x(1),x(2));
        for i = 3:length(x)
            temp = fun(temp,x(i));
        end
        varargout{1} = temp;
    end
else
    dim = size(varargin{1});
    % Non-scalar element?
    for i = 2:length(varargin)
        if prod(size(varargin{i})) > prod(size(dim)) 
            dim = size(varargin{i});
        end
    end
    
    if max(dim) > 1
        % Vectorized operator
        % First normalize to vectors
        for i = 1:nargin-1
            if numel(varargin{i}) == 1
                varargin{i} = repmat(varargin{i},prod(dim),1);
            else
                varargin{i} = reshape(varargin{i},prod(dim),1);
            end
        end
        % And now for every element, create an operator
        result = [];
        for j = 1:prod(dim)
            for i = 1:nargin-1
                temp{i} = varargin{i}(j);
            end
            result = [result; yalmip('define',func2str(fun),temp{:})];
        end
        % Devectorize to original shape
        varargout{1} = reshape(result,dim);
    else
        % Scalar operator(x,y,z,...)
        varargout{1} = yalmip('define',func2str(fun),varargin{:});
    end
end