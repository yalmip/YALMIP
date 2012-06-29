function x = horzcat(varargin)
% horzcat (overloaded)

% Author Johan Löfberg
% $Id: vertcat.m,v 1.3 2006-07-28 06:27:01 joloef Exp $

x = varargin{1};

for i = 2:nargin
    y = varargin{i};

    dim1 = size(x);
    dim2 = size(y);
    m1 = size(x,2);
    m2 = size(y,2);
    if m1 == m2
        z1 = reshape(x(:),dim1(1),[]);
        z2 = reshape(y(:),dim2(1),[]);
        dim = dim1;
        dim(1) = dim1(1) + dim2(1);
        x =  reshape([z1;z2],dim);                
    else
        error('Dimension mismatch');
    end
end