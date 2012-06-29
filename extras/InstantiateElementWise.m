function y = InstantiateElementWise(funName,varargin);

if length(varargin{1}) == 1
    y = yalmip('define',funName,varargin{:});
else
    y = [];
    args = varargin{1};
    dims = size(args);
    args = args(:);
    for i = 1:length(args)
        varargin{1} = args(i);
        y = [y;yalmip('define',funName,varargin{:})];
    end
    y = reshape(y,dims);
end