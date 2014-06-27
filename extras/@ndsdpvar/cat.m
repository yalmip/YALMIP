function z = cat(varargin)
% cat (overloaded)

along = varargin{1};

x = varargin{2};
if nargin > 3
    y = cat(along,varargin{3:end});
else
    y = varargin{3};
end

dimx = size(x);
dimy = size(y);

if length(dimx) < along
    error('First argument does not have sufficiently many dimensions');
end
if length(dimy) < along
    error('Second argument does not have sufficiently many dimensions');
end

aux1 = dimx;aux1(along)=[];
aux2 = dimy;aux2(along)=[];
if ~isequal(aux1,aux2)
    error('Dimension mismatch');
end

xindex = reshape(1:prod(dimx),dimx);
yindex = reshape(prod(dimx)+1:prod(dimx)+prod(dimy),dimy);
zindex = cat(along,xindex,yindex);
zdim = size(zindex);
zindex = zindex(:);
[~,locx]=ismember(zindex,xindex(:));
[~,locy]=ismember(zindex,yindex(:));

z = sdpvar(length(zindex),1);
z(find(locx))=sdpvar(x(:));
z(find(locy))=sdpvar(y(:));
z = reshape(z,zdim);
       