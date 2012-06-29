function [Q,c,f,x,info] = quaddecomp(p,z)
%QUADDECOMP Internal function to decompose quadratic expression

% Author Johan Löfberg 
% $Id: quaddecomp.m,v 1.1 2006-08-10 18:00:22 joloef Exp $   

[n,m]=size(p);
info = 0;

% Is it a scalar
if (n*m==1)
    % Involved in polynomial expression
    mt = yalmip('monomtable');
    x_lin = getvariables(p);
    x_var = find(any(mt(x_lin,:),1));    
    if nargin==2
        x_var = union(x_var,depends(z));
    end
    x = recover(x_var);  
    if is(p,'linear')
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
    degrees = sum(mt(x_lin,:),2);
    if all(degrees<=2)
        Q = spalloc(length(x),length(x),2*nnz(getbase(p)));
        c = zeros(length(x),1);
        base = getbase(p);
        if nnz(base(1))==0
            f = 0;
            base = base(2:end);
        else
            f = base(1);
            base = base(2:end);
        end
        
        mt = mt(x_lin,x_var);
        
        if 1
            [jj,ii,kk] = find(mt');ii = [ii(:) ;0];
            top = 1;
            for i = 1:length(x_lin)
                if ii(top) == ii(top+1)
                    j = [jj(top) jj(top+1)];
                    top = top + 2;
                else
                    j = [jj(top)];
                    top = top + 1;
                end
                if length(j)==1    % One variable
                    if kk(top-1)==1
                        c(j)=base(i);
                    else
                        Q(j,j)=Q(j,j) + base(i)/2;
                    end
                else
                    Q(j(1),j(2))=Q(j(1),j(2)) + base(i)/2;
                end
            end

        else
            for i = 1:length(x_lin)
                j = find(mt(i,:)); % What variables are used
                if length(j)==1    % One variable
                    if mt(i,j)==1
                        c(j)=base(i);
                    else
                        Q(j,j)=Q(j,j) + base(i)/2;
                    end
                else
                    Q(j(1),j(2))=Q(j(1),j(2)) + base(i)/2;
                end
            end
        end
        Q = Q+Q';
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
