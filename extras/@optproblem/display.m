function display(X)
%display           Overloaded

%disp(['Optimization problem']);
X.Options.verbose=max(0,X.Options.verbose - 1);
if ~isempty(X.Constraints)
    if any(is(X.Constraints,'uncertain'))
        X = robustify(X);
    end
end
[p,recoverdata,solver,diagnostic,F,Fremoved] = compileinterfacedata(X.Constraints,[],[],X.Objective,X.Options,1,0);

disp(' ');
type = 'Optimization';
if isempty(X.Objective)
    type = 'Feasibility';
end
if length(p.aux_variables > 0)
    disp([type ' problem with ' num2str(nnz(p.variabletype==0)) ' variables (' num2str( length(p.aux_variables)) ' introduced by YALMIP)']);
else
     disp([type ' problem with ' num2str( nnz(p.variabletype==0)) ' variables ']);
end
    
disp(['                          ' num2str((p.K.f)) ' equality constraints ']);
disp(['                          ' num2str((p.K.l)) ' scalar inequalities ']);
if any(p.K.s)
    disp(['                          ' num2str(length(p.K.s)) ' semidefinite constraints ']);
else
     disp(['                          ' num2str(p.K.s) ' semidefinite constraints ']);
end

if p.ProblemClass.objective.linear 
    if any( p.c(p.evalVariables))
        disp(['Objective is general nonlinear']);        
    else
        disp(['Objective is linear']);
    end
elseif p.ProblemClass.objective.quadratic.convex 
    if any( p.Q(p.evalVariables),2)
    else
    disp(['Objective is convex quadratic']);
    end
elseif p.ProblemClass.objective.quadratic.nonconvex 
    disp(['Objective is nonconvex quadratic']);
elseif p.ProblemClass.objective.polynomial
    disp(['Objective is polynomial']);
elseif p.ProblemClass.objective.sigmonial
    disp(['Objective is sigmonial']);
end

if length(p.binary_variables)> 0 |  length(p.integer_variables)> 0
    disp([num2str( length(p.binary_variables)) ' binary variables and ' num2str( length(p.integer_variables)) ' integer variables'])
end





