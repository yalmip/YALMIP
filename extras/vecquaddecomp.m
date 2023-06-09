function [Q,c,f,x,info] = vecquaddecomp(p,z)
%QUADDECOMP Internal function to decompose quadratic expression

[n,m]=size(p);
info = 0;

Q = cell(n*m,1);
c = cell(n*m,1);
f = cell(n*m,1);

% Involved in polynomial expression
[mt,variabletype] = yalmip('monomtable');
x_lin = getvariables(p);
x_var = find(any(mt(x_lin,:),1));
if nargin==2
    x_var = union(x_var,depends(z));
end
x = recover(x_var);

basis = getbase(p);
vars = getvariables(p);
vt = variabletype;

mt = mt(x_lin,x_var);
variabletype = variabletype(x_lin);

% If all monomials are different. we can do things faster
hash = mt*randn(size(mt,2),1);
simple =  length(unique(hash))==length(hash);

linear    = find(variabletype==0);
quadratic = find(variabletype==2);
bilinear  = find(variabletype==1);
[il,jl,kl] = find(mt(linear,:));
[jb,ib,kb] = find(mt(bilinear,:)');
[iq,jq,kq] = find(mt(quadratic,:));
bilinindex = sub2ind([length(x) length(x)],jb(1:2:end),jb(2:2:end));
quadindex = sub2ind([length(x) length(x)],jq,jq);

for index = 1:n*m
    base = basis(index,:);
    if all(vt(vars(find(base(2:end)))) == 0)
        n = length(x);
        Q{index} = spalloc(n,n,0);
        fc = basis(index,:);
        f{index} = fc(1);
        c{index} = zeros(length(x),1);
        
        [aux,loc] = ismember(vars,x_var);
        i = 1+(1:length(vars));
        nz = find(loc);
        loc = loc(nz);
        i = i(nz);
        c{index}(loc) = fc(i);
    else
        if all(vt(vars(find(base(2:end)))) <= 2)
            Qtemp = spalloc(length(x),length(x),2*nnz(base));
            ctemp = spalloc(length(x),1,0);            
            f{index} = base(1);
            base = base(2:end);                        
            % Old version of code
            [jj,ii,kk] = find(mt');ii = [ii(:) ;0];
            bilinear_places = find(diff(ii)==0);
            quadratic_places   = find(diff(ii)~=0 & kk==2);
            Qtemp(sub2ind(size(Qtemp),jj(bilinear_places),jj(bilinear_places+1)))= base(ii(bilinear_places))/2;
            Qtemp(sub2ind(size(Qtemp),jj(quadratic_places),jj(quadratic_places)))= base(ii(quadratic_places))/2;
            linear_places = find(variabletype==0);
            [iii,jjj] = find(mt(linear_places,:));
            ctemp(jjj) = base(linear_places);
            
            Qtemp =  Qtemp+Qtemp';
            
            Q{index} = Qtemp;
            c{index} = ctemp;
        else
            info = 1;
            x = [];
        end
    end
end

