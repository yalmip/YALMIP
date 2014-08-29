function y = blkdiag(varargin)
%BLKDIAG (overloaded)

if nargin<2
    y=varargin{1};
    return
end

% Get dimensions
n = zeros(nargin,1);
m = zeros(nargin,1);
isasdpvar = zeros(nargin,1);
Symmetric = zeros(nargin,1);
for i = 1:nargin
    if isa(varargin{i},'sdpvar')
        isasdpvar(i) = 1;
        n(i)=varargin{i}.dim(1);
        m(i)=varargin{i}.dim(2);
    else
        [n(i) m(i)] = size(varargin{i});
    end
    Symmetric(i)= issymmetric(varargin{i});
end

% Find all free variables used
all_lmi_variables = [];
for i = 1:nargin
    all_lmi_variables = [all_lmi_variables getvariables(varargin{i})];
end
all_lmi_variables = unique(all_lmi_variables);

% Create an SDPVAR
y=sdpvar(1,1,'rect',all_lmi_variables,[]);

% Allocate a basis
nn = 0;
for i = 1:nargin
    if isa(varargin{i},'sdpvar')
        nn = nn+nnz(varargin{i}.basis);
    else
        nn = nn+nnz(varargin{i});
    end
end
% We work with transposed basis for speed reasons
% (still slow though...)
basis = spalloc(sum(m)*sum(n),1+length(all_lmi_variables),nn)';

% Some indexation tricks
msums = cumsum([0; m]);
nsums = cumsum([0; n]);
summ=sum(m);
sumn=sum(n);
indextable = reshape(1:sumn*summ,sumn,summ);

%ix = [];
%jx = [];
%sx = [];
for j = 1:nargin
    nnindex = indextable(1+nsums(j):nsums(j+1),1+msums(j):msums(j+1));
    if isasdpvar(j)
        this_uses = find(ismembc(all_lmi_variables,varargin{j}.lmi_variables));
        mindex = [1 this_uses+1];
        basis(mindex(:),nnindex(:)) = varargin{j}.basis';       
        %      basis(nnindex(:),mindex) = varargin{j}.basis;
        %[ax,bx,cx] = find(varargin{j}.basis);       
        %ix = [ix;ax + nsums(j)];
        %jx = [jx;mindex(bx)]
        %sx = [sx;cx];
    else           
        basis(1,nnindex(:)) = varargin{j}(:)';
        %[ax,bx,cx] = find(varargin{j});        
        %ix = [ix;ax + nsums(j)];
        %jx = [jx;mindex(bx)]
        %sx = [sx;cx];               
    end
end
y.basis = basis';
y.dim(1) = sumn;
y.dim(2) = summ;
% Reset info about conic terms
y.conicinfo = [0 0];


