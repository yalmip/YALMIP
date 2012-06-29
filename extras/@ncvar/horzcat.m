function y = horzcat(varargin)
%HORZCAT (overloaded)

% Author Johan Löfberg 
% $Id: horzcat.m,v 1.2 2006-08-11 11:48:15 joloef Exp $  

prenargin = nargin;
% Fast exit
if prenargin<2
    y=varargin{1};
    return
end

% Get dimensions
n = zeros(prenargin,1);
m = zeros(prenargin,1);
for i = 1:prenargin
    if isa(varargin{i},'blkvar')
        varargin{i} = sdpvar(varargin{i});
    end
    if isa(varargin{i},'sdpvar')
        varargin{i}= ncvar(struct(varargin{i}));
    end

    [n(i),m(i)]=size(varargin{i});
end

% Keep only non-empty
keep_these = find((n.*m)~=0);
if length(keep_these)<length(n)
    varargin = {varargin{keep_these}};
    n = n(keep_these);
    m = m(keep_these);
end;

% All heights should be equal
if any(n~=n(1))
    error('All matrices on a row in the bracketed expression must have the same number of rows.');
end

nblocks = size(varargin,2);

isasdpvar = zeros(nblocks,1);
for i = 1:nblocks
    isasdpvar(i) = isa(varargin{i},'ncvar');
    isachar(i)   = isa(varargin{i},'char');
end

% Finish if this is a symbolic expression
% including '?' operators
if any(isachar)
    y = blkvar;
    for i = 1:nargin
        if isachar(i)
            switch varargin{i}
                case {'i','I'}
                    y(1,i) = 1;
                case {'s'}
                case 'z'
                    y(1,i) = 0;
                otherwise
            end
        else
            y(1,i) = varargin{i};
        end
    end
    return
end

% Find all free variables used
all_lmi_variables = [];
for i = 1:nblocks
    if isasdpvar(i)
        all_lmi_variables = [all_lmi_variables varargin{i}.lmi_variables];
    end
end
all_lmi_variables = uniquestripped(all_lmi_variables);

% Pick one of the sdpvar objects to build on...
y = varargin{min(find(isasdpvar))};

% Some indexation tricks
n = n(1);

basis_i = [];
basis_j = [];
basis_s = [];
shft = 0;
for j = 1:nblocks
    if isasdpvar(j)
        in_this = find(ismembc(all_lmi_variables,varargin{j}.lmi_variables));
        dummy = [1 1+in_this];
        [i2,j2,s2] = find(varargin{j}.basis);
        j2 = dummy(j2);
        add_shift = size(varargin{j}.basis,1);
    else
        [i2,j2,s2] = find(varargin{j}(:));
        add_shift = size(varargin{j}(:),1);
    end
        basis_i = [basis_i;i2(:)+shft];
        basis_j = [basis_j;j2(:)];
        basis_s = [basis_s;s2(:)];
        shft = shft + add_shift;   
end
basis = sparse(basis_i,basis_j,basis_s,sum(m)*n,1+length(all_lmi_variables));

y.dim(1) = n;
y.dim(2) = sum(m);
y.basis = basis;
y.lmi_variables = all_lmi_variables;
% Reset info about conic terms
y.conicinfo = [0 0];

