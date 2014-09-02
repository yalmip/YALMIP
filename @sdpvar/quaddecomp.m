function [Q,c,f,x,info] = quaddecomp(p,z)
%QUADDECOMP Internal function to decompose quadratic expression

[n,m]=size(p);
info = 0;

% Is it a scalar
if (n*m==1)
    % Involved in polynomial expression
    [mt,variabletype] = yalmip('monomtable');
    x_lin = getvariables(p);
    x_var = find(any(mt(x_lin,:),1));
    if nargin==2
        x_var = union(x_var,depends(z));
    end
    x = recover(x_var);
    if all(variabletype(x_lin) ==0)% is(p,'linear')
        n = length(x);
        Q = spalloc(n,n,0);
        fc = getbase(p);
        f = fc(1);
        if nargin==2
            vars = getvariables(p);
            c = zeros(length(x),1);
            for i = 1:length(vars)
                c(find(vars(i)==x_var)) = fc(1+i);
            end
        else
            c = fc(2:end);c=c(:);
        end
        return
    end
    variabletype = variabletype(x_lin);
    if all(variabletype<=2)

        base = getbase(p);
        if nnz(base(1))==0
            f = 0;
            base = base(2:end);
        else
            f = base(1);
            base = base(2:end);
        end
        mt = mt(x_lin,x_var);     
        quads   = find (variabletype == 2);
        bilins  = find (variabletype == 1);
        linears  = find (variabletype == 0);
        [varsC,aux1,aux2] = find(mt(linears,:)');
        [varsQ,aux1,aux2] = find(mt(quads,:)');
        [varsB,aux1,aux2] = find(mt(bilins,:)');
        if isempty(varsQ)
            varsQ = [];
        end
        if isempty(varsB)
            varsB = [];
        end
        if isempty(varsC)
            varsC = [];
        end
        c = sparse(varsC,1,base(linears),length(x),1);
        ii = [varsQ ; varsB(1:2:end) ; varsB(2:2:end)];
        jj = [varsQ ; varsB(2:2:end) ; varsB(1:2:end)];
        kk = [base(quads)  base(bilins)/2  base(bilins)/2];
        Q = sparse(ii,jj,kk,length(x),length(x));
    else
        if nargout==5
            info = 1;
            Q = [];
            c = [];
            f = [];
            x = [];
        else
            error('Function is not quadratic');
        end
    end

else
    if nargout==5
        info = 1;
        Q = [];
        c = [];
        f = [];
        x = [];
    else
        error('Function is not scalar');
    end
end
