function y = horzcat(varargin)
%HORZCAT (overloaded)

prenargin = nargin;
% Fast exit
if prenargin<2
    y=varargin{1};
    return
end

if nargin>20
    y = horzcat(horzcat(varargin{1:fix(nargin/2)}),horzcat(varargin{fix((nargin/2))+1:end}));
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

% All heights should be equal
if any(n~=n(1))
    error('All matrices on a row in the bracketed expression must have the same number of rows.');
end

nblocks = size(varargin,2);

isasdpvar = zeros(nblocks,1);
isachar = zeros(nblocks,1);
for i = 1:nblocks
    if isa(varargin{i},'sdpvar')
        isasdpvar(i) = 1;
    elseif isa(varargin{i},'char');
        isachar(i) = 1;
    end
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
            h = size(varargin{i}.rightfactors{j},1);
            y.rightfactors{end+1} = [spalloc(h,sum(m(1:1:i-1)),0) varargin{i}.rightfactors{j} spalloc(h,sum(m(i+1:1:end)),0)];
            y.midfactors{end+1} = varargin{i}.midfactors{j};
            y.leftfactors{end+1} = varargin{i}.leftfactors{j};
        end
    elseif isnumeric(varargin{i})
        if ~all(varargin{i}==0)
            %  if ~doublehere
            here = length(y.midfactors)+1;
            % doublehere = [doublehere here];
            %  end
            y.rightfactors{here} = [spalloc(m(i),sum(m(1:1:i-1)),0) speye(m(i)) spalloc(m(i),sum(m(i+1:1:end)),0)];
            y.midfactors{here}  = varargin{i};
            y.leftfactors{here}  = speye(size(varargin{i},1));
        end
    end
end
y = cleandoublefactors(y);
y = flushmidfactors(y);