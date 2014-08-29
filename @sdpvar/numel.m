function N=numel(varargin)
%NUMEL (overloaded)

X = varargin{1};
N = prod(size(X));
