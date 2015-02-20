function varargout = vectorizedlogic(fun,varargin);
if nargin-1 == 1
    if length(varargin{1})==1
       % What should be returned on OR(1) etc?
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
    if max(dim) > 1
        % Vectorized xor
        % First normalize to vectors
        for i = 1:nargin-1
            varargin{i} = reshape(varargin{i},prod(dim),1);
        end
        % And now for every element, create an or operator
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
        % Scalar xor(x,y,z,...)
        varargout{1} = yalmip('define',func2str(fun),varargin{:});
    end
end