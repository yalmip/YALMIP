function  output = bnb_solvelower(lowersolver,relaxed_p,upper,lower,x_min,aggresiveprune,allSolutions)

if all(relaxed_p.lb==relaxed_p.ub)
    x = relaxed_p.lb;
    if checkfeasiblefast(relaxed_p,relaxed_p.lb,relaxed_p.options.bnb.feastol)
        output.problem = 0;
    else
        output.problem = 1;
    end
    output.Primal = x;
    return
end

if ~(relaxed_p.all_integers && all(relaxed_p.c == fix(relaxed_p.c)) && nnz(relaxed_p.Q)==0)
     % Objective contains floating-point numbers, so add some margin
     % if we add an  upper bound cut
     upper = upper + 1e-4;
end

p = relaxed_p;
p.solver.tag = p.solver.lower.tag;

if ~isinf(upper) & nnz(p.Q)==0 & isequal(p.K.m,0) && ~any(p.variabletype)
    if p.all_integers && all(p.c == fix(p.c))
        % All integer objective coefficients and all integer
        % variables, we must find a solution which is at least
        % 1 better than current upper bound
        p.F_struc = [p.F_struc(1:p.K.f,:);upper-1-p.f -p.c';p.F_struc(1+p.K.f:end,:)];
        p.K.l=p.K.l+1;       
    else
       p.F_struc = [p.F_struc(1:p.K.f,:);upper-p.f -p.c';p.F_struc(1+p.K.f:end,:)];      
       p.K.l=p.K.l+1;
    end
end

% Exclusion cuts for negated binaries based on some optimal solutions
if ~isinf(upper) && p.all_integers && all(p.ub <= 0) && all(p.lb >= -1)   
    for i = 1:min(size(allSolutions,2),10);
        [b,a] = exclusionCut(allSolutions(:,end-i+1),-1);
        p.F_struc = [p.F_struc(1:p.K.f,:);b a;p.F_struc(1+p.K.f:end,:)];
        p.K.l=p.K.l+1;
    end    
end

removethese = p.lb==p.ub;
if nnz(removethese)>0 & all(p.variabletype == 0) & isempty(p.evalMap)% ~isequal(lowersolver,'callfmincongp') & ~isequal(lowersolver,'callgpposy')
 
    if ~isempty(p.F_struc)
        if ~isequal(p.K.l,0) & p.options.bnb.ineq2eq
            affected = find(any(p.F_struc(:,1+find(removethese)),2));
        end
        p.F_struc(:,1)=p.F_struc(:,1)+p.F_struc(:,1+find(removethese))*p.lb(removethese);
        p.F_struc(:,1+find(removethese))=[];       
    end
        
    idx = find(removethese);
    p.f = p.f + p.c(idx)'*p.lb(idx);
    p.c(idx)=[];
    if nnz(p.Q)>0     
        p.c = p.c + 2*p.Q(find(~removethese),idx)*p.lb(idx);
        p.f = p.f + p.lb(idx)'*p.Q(idx,idx)*p.lb(idx);
        p.Q(:,find(removethese))=[];
        p.Q(find(removethese),:)=[];
    else
        p.Q = spalloc(length(p.c),length(p.c),0);
    end
    p.lb(removethese)=[];
    p.ub(removethese)=[];
    p.x0(removethese)=[];
    p.monomtable(:,find(removethese))=[];
    p.monomtable(find(removethese),:)=[];
    
    % This is not necessarily correct!! x*y^2, fix y and we have a linear!
    p.variabletype(removethese) = [];
    % p.variabletype = []; % Reset, to messy to recompute
    if ~isequal(p.K.l,0) & p.options.bnb.ineq2eq
        Beq = p.F_struc(1:p.K.f,1);
        Aeq = -p.F_struc(1:p.K.f,2:end);
        B   = p.F_struc(1+p.K.f:p.K.l+p.K.f,1);
        A   = -p.F_struc(1+p.K.f:p.K.l+p.K.f,2:end);
        affected = affected(affected <= p.K.f + p.K.l);
        affected = affected(affected > p.K.f) - p.K.f;
        aaa = zeros(p.K.l,1);aaa(affected) = 1;
        A1 = A(find(~aaa),:);
        B1 = B(find(~aaa),:);
        [A,B,Aeq2,Beq2,index] = ineq2eq(A(affected,:),B(affected));
        if ~isempty(index)
            actuallyused = find(any([Aeq2,Beq2],2));
            Beq2 = Beq2(actuallyused);if size(Beq2,1)==0;Beq2 = [];end
            Aeq2 = Aeq2(actuallyused,:);if size(Aeq2,1)==0;Aeq2 = [];end
            p.F_struc = [Beq -Aeq;Beq2 -Aeq2;B1 -A1;B -A;p.F_struc(1+p.K.f + p.K.l:end,:)];
            p.K.f = length(Beq) + length(Beq2);
            p.K.l = length(B) + length(B1);
        end
    end
    
    % Find completely empty rows
    zero_row = find(~any(p.F_struc,2));    
    zero_row = zero_row(zero_row <= p.K.f + p.K.l);
    if ~isempty(zero_row)
        p.F_struc(zero_row,:) = [];
        p.K.l =  p.K.l - nnz(zero_row > p.K.f);
        p.K.f =  p.K.f - nnz(zero_row <= p.K.f);
    end
    if p.K.l > 0
         zero_row = find(~any(p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end),2));
         if ~isempty(zero_row)
             lhs = p.F_struc(p.K.f + zero_row,1);
             zero_row_pos = find(lhs >= 0);
             remove_these = zero_row(zero_row_pos);
             p.F_struc(p.K.f + remove_these,:) = [];
             p.K.l = p.K.l - length(remove_these);
             zero_row_neg = find(lhs < -p.options.bnb.feastol);
             if ~isempty(zero_row_neg)
                output.problem = 1;    
                output.Primal = p.lb; 
             end
         end
    end
    
%     if relaxed_p.K.q(1) > 0
%         top = 1 + p.K.f + p.K.l;
%         dels = zeros(sum(p.K.q),1);
%         deli = zeros(length(p.K.q),1);
%         for j = 1:length(p.K.q)
%             m = p.K.q(j);
%             M = p.F_struc(top:top+1,:);
%             if nnz(M)==4
%                 if all(M(:,1)==0)
%                     [ii,jj,s] = find(M(1,:));
%                     if all(s==1)
%                         [ii2,jj2,s2] = find(M(2,:));
%                         if isequal(s2,[1 -1]) && isequal(ii,ii2) && isequal(jj,jj2)
%                            if any(isinf(p.ub(jj-1)))
%                                deli(j)=1;
%                                dels(top:top+p.K.q(j)-1)=1;
%                            end
%                         end
%                     end
%                 end
%             end
%             top = top + p.K.q(j);
%         end
%         if any(deli)
%              p.F_struc(find(dels),:)=[];
%              p.K.q(find(deli))=[];
%              if length(p.K.q) == 0
%                  p.K.q = 0;
%              end
%         end
%     end
    
    if p.K.s(1) > 0
        top = 1 + p.K.f + p.K.l + sum(p.K.q);
        for j = 1:length(p.K.s)
            X = p.F_struc(top:top + p.K.s(j)^2-1,:); 
            X = reshape(any(X,2),p.K.s(j),p.K.s(j));
            e = find(~any(X,2));
            if any(e)
                Z = spalloc(p.K.s(j),p.K.s(j),length(e)*2*p.K.s(j));
                for k = 1:length(e);
                    Z(:,e(k))=1;
                    Z(e(k),:)=1;
                end            
                m = find(Z(:));
                p.F_struc(top + m - 1,:)=[];
                p.K.s(j) = p.K.s(j) - length(e);
            end
            top = top + p.K.s(j)^2;
        end
    end
                          
    % Derive bounds from this model, and if we fix more variables, apply
    % recursively  
    if isempty(p.F_struc)
        lb = p.lb;
        ub = p.ub;
    else
        [lb,ub] = findulb(p.F_struc,p.K);
    end
    newub = min(ub,p.ub);
    newlb = max(lb,p.lb);
    if any(newub == newlb)
        dummy = p;
        dummy.lb = newlb;
        dummy.ub = newub;
        output = bnb_solvelower(lowersolver,dummy,inf,lower,x_min,aggresiveprune,[]);        
    else
        if any(p.lb>p.ub+0.1)
            output.problem = 1;
            output.Primal = zeros(length(p.lb),1);
        else
            p.solver.version = p.solver.lower.version;
            p.solver.subversion = p.solver.lower.subversion;                                                                          
            output = feval(lowersolver,p);                       
        end
    end
    x=relaxed_p.c*0;
    x(removethese)=relaxed_p.lb(removethese);
    x(~removethese)=output.Primal;
    output.Primal=x;
else
    p.solver = p.solver.lower;
    output = feval(lowersolver,p);    
end



function [A, B, Aeq, Beq, ind_eq] = ineq2eq(A, B)
% Copyright is with the following author(s):
%
% (C) 2006 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%          loefberg@control.ee.ethz.ch
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

[ne, nx] = size(A);
Aeq = [];
Beq = [];
ind_eq = [];
if isempty(A)
    return
end
sumM = sum(A, 2) + B;
for ii = 1:ne-1,   
    s = sumM(1);
    
    % get matrix which contains all rows starting from ii+1 
    sumM = sumM(2:end,:);
    
    % possible candidates are those rows whose sum is equal to the sum of the
    % original row
    possible_eq = find(abs(sumM + s) < 1e-12);
    if isempty(possible_eq),
        continue
    end
    possible_eq = possible_eq + ii;
    b1 = B(ii);    
    a1 = A(ii,  :) ;
    
    % now compare if the two inequalities (the second one with opposite
    % sign) are really equal (hence they form an equality constraint)
    for jj = possible_eq',
        % first compare the B part, this is very cheap
        if abs(b1 + B(jj)) < 1e-12,
            % now compare the A parts as well
            if norm(a1 + A(jj,  :) , Inf) < 1e-12,
                % jj-th inequality together with ii-th inequality forms an equality
                % constraint
                ind_eq = [ind_eq; ii jj];
                break
            end
        end
    end
end

if isempty(ind_eq),
    % no equality constraints
    return
else
    % indices of remaining constraints which are inequalities
    ind_ineq = setdiff(1:ne, ind_eq(:));
    Aeq = A(ind_eq(:,1),  :) ;
    Beq = B(ind_eq(:,1));
    A = A(ind_ineq,  :) ;
    B = B(ind_ineq);
end