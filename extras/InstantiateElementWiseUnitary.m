function y = InstantiateElementWiseUnitary(funName,varargin);

args = varargin{1};
dims = size(args);
y = yalmip('defineVectorizedUnitary',funName,varargin{:});
y = reshape(y,dims);   
