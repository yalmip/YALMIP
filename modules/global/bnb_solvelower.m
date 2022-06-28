function output = bnb_solvelower(lowersolver,relaxed_p,upper,lower,x_min,allSolutions)

if all(relaxed_p.lb==relaxed_p.ub)
    x = relaxed_p.lb;
    if checkfeasiblefast(relaxed_p,relaxed_p.lb,relaxed_p.options.bnb.feastol)
        output.problem = 0;
    else
        output.problem = 24;
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
        if upper == lower + 1
            p = addEquality(p,[upper-1-p.f -p.c']);
        else
            p = addInequality(p,[upper-1-p.f -p.c']);        
        end        
    end
end

if ~isinf(upper) && p.all_integers && all(p.ub <= 0) && all(p.lb >= -1)   
    % Exclusion cuts for negated binaries based on some optimal solutions
    % kept for historical reason on malformed model
    for i = 1:min(size(allSolutions,2),10)
        [b,a] = exclusionCut(allSolutions(:,end-i+1),-1);
        p = addInequality(p,[b a]);        
    end    
elseif ~isinf(upper) && p.all_integers && all(p.ub <= 1) && all(p.lb >= 0)   
    % Normal exclusion on binary
     for i = 1:min(size(allSolutions,2),10)
         [b,a] = exclusionCut(allSolutions(:,end-i+1),1);
         p = addInequality(p,[b a]);        
     end    
end

% for i = 1:length(p.cardinalityvariables)
%     used = p.cardinalityvariables{i};
%     L = p.lb(used);
%     if sum(L) == p.cardinalitysize{i}
%         p.ub(used(p.lb(used) < p.ub(used) )) = 0;  
%     end
% end

removethese = p.lb==p.ub;
[~,map] = ismember(find(~p.lb==p.ub),1:length(p.c));
if nnz(removethese)>0 && all(p.variabletype == 0) && isempty(p.evalMap) && p.options.allowsmashing
 
    % Fixed variables, so let us try to presolve as muh as possible
    % (needed often in SDPs etc where solvers are less robust to weird
    % models having no interior etc)    
    p = smashFixed(p,'delete');
    p = smashQPOjective(p,removethese);    
    idx = find(removethese);    
    
    % FIX: should be updated
    % safe fix now where we just remove
    p.cardinalityvariables = [];
    p.cardinalitygroups = [];
    %p.atmost.groups = [];
    %p = smashAtmost(p,idx);
    
    p.lb(idx)=[];
    p.ub(idx)=[];
    if ~isempty(p.x0)
        p.x0(idx)=[];
    end
    p.monomtable(:,idx)=[];
    p.monomtable(idx,:)=[];
    p.variabletype(idx) = [];  
       
    [p,infeasible] = detectRedundantInfeasibleSOCPRows(p);
    if infeasible
        output = createOutputStructure(24);
        return
    end
    
    [p,infeasible] = detectRedundantInfeasibleSDPRows(p);
    if infeasible
        output = createOutputStructure(24);
        return
    end
    
    [p,infeasible] = detectRedundantInfeasibleEXPRows(p);
    if infeasible
        output = createOutputStructure(24);
        return
    end
    
    [p,infeasible] = detectRedundantInfeasiblePOWRows(p);
    if infeasible
        output = createOutputStructure(24);
        return
    end
        
    % We do this last, as the SOCP/EXP/SDP presolve might add trivial
    % equalities from presolving 0 >= norm(z) etc
    p = removeEmptyLPRows(p);
    [p,infeasible] = detectRedundantInfeasibleLPRows(p);
    if infeasible
        output = createOutputStructure(24);
        return
    end
                     
    % Derive bounds from this presolved model, and if we detect new fixed
    % variables, apply recursively   
    if ~isempty(p.F_struc)
        [lb,ub] = find_lp_bounds(p.F_struc,p.K,p.lb,p.ub,1);    
        p.ub = min(ub,p.ub);
        p.lb = max(lb,p.lb);
    end
    
    if any(p.lb > p.ub)
        output = createOutputStructure(24);
        return
    end
    
    if any(p.lb == p.ub) 
        % Recurse!
        output = bnb_solvelower(lowersolver,p,inf,lower,x_min,[]);
    elseif any(p.lb > p.ub-p.options.bnb.feastol)
        % Infeasible
        output = createOutputStructure(24);
        return
    else
        % Solve relaxation
        p.solver.version = p.solver.lower.version;
        p.solver.subversion = p.solver.lower.subversion;        
        output = feval(lowersolver,p);        
    end
    if output.problem == 1 || output.problem == 24
        output.Primal = [];
        return
    end
    % Recover
    x=relaxed_p.c*0;
    x(removethese)=relaxed_p.lb(removethese);
    x(~removethese)=output.Primal;
    output.Primal=x;
else
    p.solver = p.solver.lower;
    output = feval(lowersolver,p);        
end