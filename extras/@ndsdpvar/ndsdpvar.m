function X = ndsdpvar(varargin)
% NDSDPVAR Constructor for multi-dimensional SDPVAR object

% Author Johan Löfberg
% $Id: ndsdpvar.m,v 1.10 2007-01-04 08:52:49 joloef Exp $

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
NewVars = v0+(1:nNewVars);
appendYALMIPvariables(NewVars);

X = struct(X);
X.dim = n;
X.basis = [spalloc(size(X.basis,1)*N,1,0) kron(X.basis(:,2:end),eye(N))];
X.lmi_variables = NewVars;

X = class(X,'ndsdpvar');
X = clean(X);

function appendYALMIPvariables(lmi_variables);

[mt,variabletype,hashed_monoms,current_hash] = yalmip('monomtable');
% Update monomtable and pre-calculated variable type
n_mt = size(mt,1);
m_mt = size(mt,2);
newmt = [];
if min(lmi_variables)>m_mt % New variables
    if size(mt,1)~=size(mt,2)
        mt(size(mt,1),size(mt,1))=0;
    end
    % This was faster before. However in recent versions of matlab, there
    % is a compiled version of blkdiag available
    % fill=spalloc(size(mt,1),length(lmi_variables),0);   
    % mt=[mt fill;fill' speye(length(lmi_variables))]; 
    if isempty(mt)        
        mt = speye(length(lmi_variables));
        newmt = mt;
    else
        newmt = speye(length(lmi_variables));
        mt=blkdiag(mt,speye(length(lmi_variables)));
    end
else
    mt(lmi_variables,lmi_variables) = speye(length(lmi_variables));
end
variabletype(1,size(mt,1)) = 0;
if ~isempty(newmt)
    new_hash = 3*rand_hash(size(mt,2),size(newmt,2),1);
    hashed_monoms = [hashed_monoms;newmt*new_hash];
    current_hash = [current_hash;new_hash];
    yalmip('setmonomtable',mt,variabletype,hashed_monoms,current_hash);
else 
    yalmip('setmonomtable',mt,variabletype);
end