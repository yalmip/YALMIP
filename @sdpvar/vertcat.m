function y = vertcat(varargin)
%VERTCAT (overloaded)

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
    [n(i),m(i)]=size(varargin{i});
end

% Keep only non-empty
keep_these = find((n.*m)~=0);
if length(keep_these)<length(n)
    varargin = {varargin{keep_these}};
    n = n(keep_these);
    m = m(keep_these);
end;

% All widths should be equal
if any(m~=m(1))
    error('All rows in the bracketed expression must have the same number of columns.');
end

nblocks = size(varargin,2);
isasdpvar = zeros(nblocks,1);
% Get dimensions
for i = 1:nblocks
    isasdpvar(i) = isa(varargin{i},'sdpvar');
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

% Some indexation tricks...
basis_i = [];
basis_j = [];
basis_s = [];
m = m(1);
sumn=sum(n);
shift=[0;cumsum(n)];
for j = 1:nblocks
    if isasdpvar(j)        
        if length(all_lmi_variables)==length(varargin{j}.lmi_variables) && all_lmi_variables(1)==varargin{j}.lmi_variables(1) &&  all_lmi_variables(end)==varargin{j}.lmi_variables(end)
            % Avoid call to ismember and find
            in_this = 1:length(all_lmi_variables);
        else
            members = ismembcYALMIP(all_lmi_variables,varargin{j}.lmi_variables);
            in_this = find(members);
        end
        dummy = [1 1+in_this];
        [i2,j2,s2] = find(varargin{j}.basis);
        j2 = dummy(j2);
    else
        [i2,j2,s2] = find(varargin{j}(:));
    end
    if ~isempty(j2)
    col = 1+floor((i2-1)/n(j));
    row = i2-(col-1)*n(j);
    basis_i = [basis_i;row(:)+shift(j)+sumn*(col(:)-1)];
    basis_j = [basis_j;j2(:)];
    basis_s = [basis_s;s2(:)];
    end
end
basis = sparse(basis_i,basis_j,basis_s,sum(n)*m,1+length(all_lmi_variables));

y.dim(1) = sum(n);
y.dim(2) = m;
y.basis = basis;
y.lmi_variables = all_lmi_variables;
% Reset info about conic terms
y.conicinfo = [0 0];
y.extra.opname='';
y.extra.createTime = definecreationtime;
y = unfactor(y);

% Update the factors
% But first, check to see that factors exist in all terms, if not simply
% exit
for i = 1:length(varargin)
    if isa(varargin{i},'sdpvar')
        if length(varargin{i}.leftfactors)==0
            y = flush(y);
            return
        end
    end
end
doublehere = [];
for i = 1:length(varargin)
    if isa(varargin{i},'sdpvar')      
        for j = 1:length(varargin{i}.leftfactors)
            w = size(varargin{i}.leftfactors{j},2);
            y.leftfactors{end+1} = [spalloc(sum(n(1:1:i-1)),w,0); varargin{i}.leftfactors{j}; spalloc(sum(n(i+1:1:end)),w,0)];
            y.midfactors{end+1} = varargin{i}.midfactors{j};
            y.rightfactors{end+1} = varargin{i}.rightfactors{j};
        end
    elseif isnumeric(varargin{i})
        here = length(y.midfactors)+1;
        doublehere = [doublehere here];      
        y.leftfactors{here} = [spalloc(sum(n(1:1:i-1)),size(varargin{i},1),0); speye(size(varargin{i},1)); spalloc(sum(n(i+1:1:end)),size(varargin{i},1),0)];
        y.midfactors{here}  = varargin{i};
        y.rightfactors{here}  = speye(size(varargin{i},2));
    end
end
y = cleandoublefactors(y);