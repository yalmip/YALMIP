function y = stackcell(dummy,blocks)
%STACKCELL Internal function for rapid concatenation
%
% [x1 x2 ... xn] is written as stackcell(sdpvar(1,1),{x1,x2,x3,...,xn})
% Why the first argument? A hack to make it a method for sdpvar...


% Author Johan Löfberg
% $Id: stackcell.m,v 1.5 2006-07-26 20:17:58 joloef Exp $

nblocks = size(blocks,2);
% Get dimensions
n = zeros(nblocks,1);
m = zeros(nblocks,1);
isasdpvar = zeros(nblocks,1);
% Get dimensions
for i = 1:nblocks
    [n(i),m(i)]=size(blocks{i});
    isasdpvar(i) = isa(blocks{i},'sdpvar');
end

if ~any(isasdpvar)
    y = blocks{1};
    for i = 2:length(blocks)
        y = [y blocks{i}];
    end
    return;
end

% Keep only non-empty
keep_these = find((n.*m)~=0);
if length(keep_these)<length(n)
    blocks = {blocks{keep_these}};
    n = n(keep_these);
    m = m(keep_these);
    isasdpvar=isasdpvar(keep_these);
    nblocks = length(n);
end;

% Find all free variables used
all_lmi_variables = [];
for i = 1:nblocks
    if isasdpvar(i)
        all_lmi_variables = [all_lmi_variables blocks{i}.lmi_variables];
    end
end
all_lmi_variables = uniquestripped(all_lmi_variables);

% Pick one of the sdpvar objects to build on...
y = blocks{min(find(isasdpvar))};

% Some indexation tricks
n = n(1);
allindextopos = (sparse(1,[1;1+all_lmi_variables(:)],1:1+length(all_lmi_variables),1,1+all_lmi_variables(end)));

basis = [];
basis_i = [];
basis_j = [];
basis_s = [];
shft = 0;
for j = 1:nblocks
    if isasdpvar(j)
        in_this = find(ismembc(all_lmi_variables,blocks{j}.lmi_variables));
        dummy = [1 1+in_this];
        [i2,j2,s2] = find(blocks{j}.basis);
        j2 = dummy(j2);
        add_shift = size(blocks{j}.basis,1);
    else
        [i2,j2,s2] = find(blocks{j}(:));
        add_shift = size(blocks{j}(:),1);
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

