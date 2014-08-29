function X = ndsdpvar(varargin)
% NDSDPVAR Constructor for multi-dimensional SDPVAR object

% Sometimes it is convenient to cast a 2D SDPVAR variables as an nD
% variable

if nargin == 1 & isa(varargin{1},'sdpvar')
    X = varargin{1};
    X = struct(X);
    X.conicinfo = [0 0];
    X = class(X,'ndsdpvar');
    return
end

type  = 'symmetric';
field = 'real';

d = 0;
i = 1;
while i<=nargin
    if isa(varargin{i},'double')
        d = d + 1;
    end
    i = i + 1;
end

n = [varargin{1:d}];

if nargin > d
    type = varargin{d+1};
    if n(1)~=n(2) & ~isempty(strmatch(type,'symmetric'))
        error('non-square matrix cannot be symmetric');
    end
else
    if n(1)==n(2)
        type = 'symmetric';
    else
        type = 'full';
    end
end
if nargin > d+1
    field = varargin{d+2};
end

v0 = yalmip('nvars');
X = sdpvar(n(1),n(2),type,field);
vars = getvariables(X);
N = prod(n(3:end));
nNewVars = length(vars)*N;
usedNewVars = v0+(1:nNewVars);
appendYALMIPvariables((vars(end)+1):usedNewVars(end));

X = struct(X);
X.dim = n;
X.basis = [spalloc(size(X.basis,1)*N,1,0) kron(speye(N),X.basis(:,2:end))];
X.lmi_variables = usedNewVars;

X = class(X,'ndsdpvar');
X = clean(X);