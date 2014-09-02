function Y=reshape(varargin)
%RESHAPE (overloaded)

Y = varargin{1};
% simple hack to figure out new sizes
temp = reshape(ones(size(Y.basis,1),1),varargin{2:end});
Y.conicinfo = [0 0];
Y.dim = size(temp);
Y = clean(Y);