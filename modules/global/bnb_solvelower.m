function  output = bnb_solvelower(lowersolver,relaxed_p,upper,lower,x_min,allSolutions)

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
     % if we add an upper bound cut
     upper = upper + 1e-4;
end

p = relaxed_p;
p.solver.tag = p.solver.lower.tag;

if ~isinf(upper) && nnz(p.Q)==0 && isequal(p.K.m,0) && ~any(p.variabletype)
    if p.all_integers && all(p.c == fix(p.c))
        % All integer objective coefficients and all integer
        % variables, we must find a solution which is at least
        % 1 better than current upper bound
        p.F_struc = [p.F_struc(1:p.K.f,:);upper-1-p.f -p.c';p.F_struc(1+p.K.f:end,:)];
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
if nnz(removethese)>0 && all(p.variabletype == 0) && isempty(p.evalMap)
 
    if ~isempty(p.F_struc)    
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
    if ~isempty(p.x0)
        p.x0(removethese)=[];
    end
    p.monomtable(:,find(removethese))=[];
    p.monomtable(find(removethese),:)=[];
    
    % This is not necessarily correct!! x*y^2, fix y and we have a linear!
    p.variabletype(removethese) = [];    
    % Find completely empty rows
    zero_row = find(~any(p.F_struc,2));    
    zero_row = zero_row(zero_row <= p.K.f + p.K.l);
    if ~isempty(zero_row)
        p.F_struc(zero_row,:) = [];
        p.K.l =  p.K.l - nnz(zero_row > p.K.f);
        p.K.f =  p.K.f - nnz(zero_row <= p.K.f);
    end
    if any(p.K.l)
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
    
    % Remove zero rows in SOCP cone
	if any(p.K.q)
         top = startofSOCPCone(p.K);
         for j = 1:length(p.K.q)
             m = p.K.q(j);
             M = p.F_struc(top:top+m-1,:);
             remove = any(M,2)==0;
             if any(remove)
                 remove = find(remove);
                 p.F_struc(top + remove-1,:)=[];
                 p.K.q(j) = p.K.q(j)-length(remove);         
             end
              top = top + p.K.q(j);
         end
    end
                 
    if any(p.K.s)
        top = startofSDPCone(p.K);
        newEqualities = [];
        for j = 1:length(p.K.s)
            X = p.F_struc(top:top + p.K.s(j)^2-1,:); 
            X = reshape(any(X,2),p.K.s(j),p.K.s(j));
            e = find(~any(X,2));
            if any(e)
                % Not a single element is used, so simply and reduce SDP
                Z = spalloc(p.K.s(j),p.K.s(j),length(e)*2*p.K.s(j));
                for k = 1:length(e);
                    Z(:,e(k))=1;
                    Z(e(k),:)=1;
                end            
                m = find(Z(:));
                p.F_struc(top + m - 1,:)=[];
                p.K.s(j) = p.K.s(j) - length(e);
            else
                % Look for zero diagonal. This means we can move all
                % nonzero elements to a zero equality                
                e = find(diag(X)==0);
                if length(e)>0
                    Z = spalloc(p.K.s(j),p.K.s(j),length(e)*2*p.K.s(j));
                    for k = 1:length(e)
                        Z(:,e(k))=1;
                        Z(e(k),:)=1;
                    end            
                    m1 = find(Z(:)); % To be removed
                    m2 = find(triu(Z,1)); % To be moved
                    equalityRows = p.F_struc(top + m2 - 1,:);
                    p.F_struc(top + m1 - 1,:) = [];
                    p.K.s(j) = p.K.s(j) - length(e);   
                    equalityRows = equalityRows(find(any(equalityRows,2)),:);
                    newEqualities = [newEqualities;equalityRows];
                end
            end
            top = top + p.K.s(j)^2;
        end
        if ~isempty(newEqualities)
            p.F_struc = [newEqualities;p.F_struc];
            p.K.f = p.K.f + size(newEqualities,1);
        end
    end
                          
    % Derive bounds from this model, and if we fix more variables, apply
    % recursively  
    if isempty(p.F_struc)
        lb = p.lb;
        ub = p.ub;
    else
        [lb,ub] = find_lp_bounds(p.F_struc,p.K);
    end
    newub = min(ub,p.ub);
    newlb = max(lb,p.lb);
    if any(newub == newlb)
        dummy = p;
        dummy.lb = newlb;
        dummy.ub = newub;
        output = bnb_solvelower(lowersolver,dummy,inf,lower,x_min,[]);        
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
