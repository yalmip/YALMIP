function p = add_convex_upper_bound_cut(p,upper,x_min)

if any(p.variabletype == 2) && ~any(p.variabletype > 2) && ~isinf(upper) && ~any(p.c(p.evalVariables))
    x = x_min(p.linears);
    n = length(x);
    if ~isempty(p.shiftedQP) && ~any(p.shiftedQP.c(find(p.variabletype)))
        Q = p.shiftedQP.Q(p.linears,p.linears);
        c = p.shiftedQP.c(p.linears);
        f = p.shiftedQP.f;
        fail = 0;
    else
        % Why would there be no shifted QP then?
        % Original polynomial lifted to QCQP
        [Q, c] =  compileQuadratic(p.c,p,0);
        Q = Q(p.linears,p.linears);
        c = c(p.linears);
        f = p.f;
        if nnz(Q)==0 || min(eig(full(Q)))<-1e-15
            c_temp = p.originalModel.c;
            f_temp = p.originalModel.f;
            p = add_convex_upper_bound_original(p,upper,x_min);
            p.originalModel.c = c_temp;
            p.originalModel.f = f_temp;
            fail = 1;
        else
            fail = 0;
        end
    end
    if ~fail 
        f0 = x'*Q*x + c'*x + p.f;
        df = 2*Q*x + c;
        % f = x'*Q*x + c'*x + f
        % f0 + df(x-x0) <= upper
        N = emptyNumericalModel;
        N.F_struc = zeros(1,1+length(p.c));
        N.K.l = 1;
        N.F_struc(1) = upper-f0+df'*x;
        N.F_struc(1+p.linears) = -df';       
        p = mergeNumericalModels(p,N);
    end
end

function p = add_convex_upper_bound_original(p,upper,x_min)

if any(p.originalModel.variabletype == 3) && ~any(p.originalModel.variabletype > 3 ) && ~any(p.originalModel.c(p.originalModel.evalVariables))
   % We might have a polynomial which has been rewritten
   % for the relaxations as a quadratic model
   % Maybe it is convex in original form
   
   higher = find(p.originalModel.variabletype>2);
   if all(ismember(higher, p.originalModel.Univariates))
       % All non-quadratic polynomials are univariate
       % Easy to check convexity
       % Which variable and what degree
       z = p.originalModel.UnivariatesList(higher,:);       
       for i = 1:length(higher)
           if even(z(i,2))
               % Even monomial must enter with postive coeff
               if p.originalModel.c(higher(i),1)<0
                   %Does not, so lower bound with worst-case
                   xmax = p.ub(higher(i));
                   p.originalModel.f=p.originalModel.f+xmax* p.originalModel.c(higher(i),1);
                   p.originalModel.c(higher(i),1) = 0;                  
               end
           else
               % Positive coeff and lb>=0 ok
               % Negative coeff and ub<=0 ok
               if (p.lb(z(i,1))>=0) && (p.originalModel.c(higher(i))>=0)
                   % ok
               elseif (p.ub(z(i,1))<=0) && (p.originalModel.c(higher(i))<=0)
                   % ok
               else
                   return
               end
           end       
       end
       % Ok, either eve, or odd with non-negative domain!
       % Stll, we have to check the quadratic part
        c = p.originalModel.c;
        d = c;d(p.originalModel.variabletype<=2)=0;
        c(p.originalModel.variabletype > 2)=0;
        [Q, c] =  compileQuadratic(c,p.originalModel,0);
        Q = Q(p.originalModel.linears,p.originalModel.linears);
        c = c(p.originalModel.linears);
        f = p.originalModel.f;
        if min(eig(full(Q)))<-1e-15
            return        
        end
        x = x_min(p.originalModel.linears);
        f0 = x'*Q*x + c'*x + f;
        f0 = f0 + sum(d(higher).*x(z(:,1)).^z(:,2));
        df = 2*Q*x + c;
        df(z(:,1)) = df(z(:,1)) + d(higher).*z(:,2).*x(z(:,1)).^(z(:,2)-1);
        % f = x'*Q*x + c'*x + f
        % f0 + df(x-x0) <= upper
        N = emptyNumericalModel;
        N.F_struc = zeros(1,1+length(p.c));
        N.K.l = 1;
        N.F_struc(1) = upper-f0+df'*x;
        N.F_struc(1+p.originalModel.linears) = -df';       
        p = mergeNumericalModels(p,N);
        
   end
end