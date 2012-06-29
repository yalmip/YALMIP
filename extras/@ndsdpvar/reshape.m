function Y=reshape(varargin)
%RESHAPE (overloaded)

% Author Johan Löfberg
% $Id: reshape.m,v 1.2 2006-07-13 19:40:59 joloef Exp $

Y = varargin{1};
% simple hack to figure out new sizes
temp = reshape(ones(size(Y.basis,1),1),varargin{2:end});
Y.conicinfo = [0 0];
Y.dim = size(temp);
Y = clean(Y);