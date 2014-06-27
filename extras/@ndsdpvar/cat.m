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

m = max([length(dimx) length(dimy) along]);
dimx = [dimx ones(1,m-length(dimx))];
dimy = [dimy ones(1,m-length(dimy))];

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

A = sparse(find(locx),1:prod(dimx),1,prod(zdim),prod(dimx));
B = sparse(find(locy),1:prod(dimy),1,prod(zdim),prod(dimy));
z = A*sdpvar(x(:)) + B*sdpvar(y(:));
       