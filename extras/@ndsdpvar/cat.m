function z = cat(varargin)
% cat (overloaded)

along = varargin{1};
x = varargin{2};
if nargin > 3
    y = cat(along,varargin{3:end});
else
    y = varargin{3};
end

nonemptyindex = find(~cellfun('isempty',{varargin{2:end}}));
if length(nonemptyindex)==1
    z = varargin{nonemptyindex+1};
    return
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

try
    locx=ismembc2(zindex,xindex(:));
    locy=ismembc2(zindex,yindex(:));
catch
    % Octave
    [~,locx]=ismember(zindex,xindex(:));
    [~,locy]=ismember(zindex,yindex(:));
end

A = sparse(find(locx),1:prod(dimx),1,prod(zdim),prod(dimx));
B = sparse(find(locy),1:prod(dimy),1,prod(zdim),prod(dimy));

if isa(x,'sdpvar')
    x = ndsdpvar(x);
end
if isa(y,'sdpvar')
    y = ndsdpvar(y);
end
if isa(x,'ndsdpvar') && isa(y,'ndsdpvar')
    if max(x.lmi_variables) < min(y.lmi_variables)
        z = x;
        z.basis = [A*x.basis(:,1)+B*y.basis(:,1) A*x.basis(:,2:end) B*y.basis(:,2:end)];
        z.lmi_variables = [x.lmi_variables y.lmi_variables];
        z.dim = zdim;        
        return
    elseif max(y.lmi_variables) < min(x.lmi_variables)
        z = x;
        z.basis = [A*x.basis(:,1)+B*y.basis(:,1) B*y.basis(:,2:end) A*x.basis(:,2:end)];
        z.lmi_variables = [y.lmi_variables x.lmi_variables];
        z.dim = zdim;        
        return
    end
end
switch class(x)
    case 'double'
        x = A*x(:);
    case 'ndsdpvar'
        x.basis = A*x.basis;x.dim = [prod(zdim) 1];   
    otherwise
end
switch class(y)
    case 'double'
        y = B*y(:);
    case 'ndsdpvar'
        y.basis = B*y.basis;y.dim = [prod(zdim) 1];  
    otherwise
end
z = x+y;
z = reshape(z,zdim);