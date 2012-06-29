function x = horzcat(varargin)
% horzcat (overloaded)

% Author Johan Löfberg
% $Id: horzcat.m,v 1.4 2006-07-28 06:27:01 joloef Exp $

x = varargin{1};

for i = 2:nargin
    y = varargin{i};

    dim1 = size(x);
    dim2 = size(y);
    m1 = size(x,1);
    m2 = size(y,1);
    if m1 == m2
        z1 = reshape(x(:),[],dim1(2));
        z2 = reshape(y(:),[],dim2(2));
        dim = dim1;
        dim(1) = dim1(2) + dim2(2);
        x =  reshape([z1;z2],dim);                
    else
        erorr('Dimension mismatch');
    end
end