function [solver,problem,forced_choice,failureMode] = selectsolver(options,ProblemClass,solvers,socp_are_really_qc,allsolvers);
%SELECTSOLVER Internal function to select solver based on problem category

problem = 0;
failureMode = '';

% UNDOCUMENTED
force_solver = yalmip('solver');
if length(force_solver)>0
    options.solver = force_solver;
end

% YALMIP has discovered in an previous call that the model isn't a GP, and
% now searches for a non-GP solver
if options.thisisnotagp
    ProblemClass.gppossible = 0;
end

% ***************************************************
% Maybe the user is stubborn and wants to pick solver
% ***************************************************
forced_choice = 0;
if length(options.solver)>0 && isempty(strfind(options.solver,'*'))
    
    if isequal(strfind(options.solver,'-'),1)
        % No, but avoid!
        options.solver = '';
        keep = find(~strcmpi({solvers.tag},'fmincon'));
        solvers = solvers(keep);
    else
        if strfind(options.solver,'+')
            forced_choice = 1;
            options.solver = strrep(options.solver,'+','');
        end
        % Create tags with version also
        temp = expandSolverName(solvers);
        
        opsolver = lower(options.solver);
        splits = strfind(opsolver,',');
        if isempty(splits)
            names{1} = opsolver;
        else
            start = 1;
            for i  = 1:length(splits)
                names{i} = opsolver(start:splits(i)-1);
                start = splits(i)+1;
            end
            names{end+1} = opsolver(start:end);
        end
        
        index1 = [];
        index2 = [];
        for i = 1:length(names)
            index1 = [index1 find(strcmpi({solvers.tag},names{i}))];
            index2 = [index1 find(strcmpi({temp.tag},names{i}))];
        end
        if isempty(index1) & isempty(index2)
            % Specified solver not found among available solvers
            % Is it even a supported solver
            temp = expandSolverName(allsolvers);
            for i = 1:length(names)
                index1 = [index1 find(strcmp(lower({allsolvers.tag}),names{i}))];
                index2 = [index1 find(strcmp(lower({temp.tag}),names{i}))];
            end
            if isempty(index1) & isempty(index2)
                problem = -9;
            else
                problem = -3;
            end
            solver = [];
            return;
        else
            solvers = solvers(union(index1,index2));
        end
    end
end

% ************************************************
% Prune based on forced global solver
% ************************************************
if  ~forced_choice && options.forceglobal
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).global;                         
    end    
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'global optimization');
end
% Prune based on avoiding global solver (in BMIBNB we don't want that as
% upper solver)
if isfield(options,'avoidGlobalSolver')
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = ~solvers(i).global;                         
    end    
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'global optimization');
end    

% ************************************************
% Prune based on objective
% ************************************************
if ProblemClass.objective.sigmonial & ~forced_choice & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.sigmonial;                         
    end    
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'signomial terms in objective');
end   

% ************************************************
% Prune based on objective
% ************************************************
if ProblemClass.objective.sigmonial & ~forced_choice & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.sigmonial;                         
    end    
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'signomial terms in objective');
end    
if ProblemClass.objective.polynomial & ~forced_choice & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.polynomial | solvers(i).objective.sigmonial;            
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'polynomial terms in objective');
end  
if ProblemClass.objective.quadratic.nonconvex & ~forced_choice & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.polynomial | solvers(i).objective.sigmonial | solvers(i).objective.quadratic.nonconvex;        
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'nonconvex quadratic terms in objective');
end  
if ProblemClass.objective.quadratic.convex & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        direct = solvers(i).objective.polynomial | solvers(i).objective.sigmonial | solvers(i).objective.quadratic.nonconvex | solvers(i).objective.quadratic.convex;
        indirect = solvers(i).constraint.inequalities.semidefinite.linear | solvers(i).constraint.inequalities.secondordercone.linear;
        if direct | indirect
            keep(i)=1;
        else
            keep(i)=0;
        end
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'quadratic terms in objective');
end  
if ProblemClass.objective.linear & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.polynomial | solvers(i).objective.sigmonial | solvers(i).objective.quadratic.nonconvex | solvers(i).objective.quadratic.convex | solvers(i).objective.linear;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'linear terms in objective');
end  

if ProblemClass.objective.maxdet.convex & ~ProblemClass.objective.linear & ~ProblemClass.objective.quadratic.convex & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.maxdet.convex | solvers(i).constraint.inequalities.semidefinite.linear;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'maxdet terms in objective');
end  

if ProblemClass.objective.maxdet.convex & (ProblemClass.objective.linear | ProblemClass.objective.quadratic.convex) & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.maxdet.convex;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'maxdet terms in objective');
end

if ProblemClass.objective.maxdet.nonconvex & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.maxdet.nonconvex;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'nonconvex maxdet terms in objective');
end


% ******************************************************
% Prune based on rank constraints
% ******************************************************
if ProblemClass.constraint.inequalities.rank & ~forced_choice & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.rank;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'rank constraints');
end  

% ******************************************************
% Prune based on semidefinite constraints
% ******************************************************
if ProblemClass.constraint.inequalities.semidefinite.sigmonial & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'signomial semidefinite constraints');
end  
if ProblemClass.constraint.inequalities.semidefinite.polynomial & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial |  solvers(i).constraint.inequalities.semidefinite.polynomial;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'polynomial semidefinite constraints');
end  
if ProblemClass.constraint.inequalities.semidefinite.quadratic & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial |  solvers(i).constraint.inequalities.semidefinite.polynomial | solvers(i).constraint.inequalities.semidefinite.quadratic;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'quadratic semidefinite constraints');
end  
if ProblemClass.constraint.inequalities.semidefinite.linear & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial |  solvers(i).constraint.inequalities.semidefinite.polynomial | solvers(i).constraint.inequalities.semidefinite.quadratic | solvers(i).constraint.inequalities.semidefinite.linear;
    end        
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'semidefinite constraints');
end

% If user has specified a, e.g., LP solver for an SDP when using OPTIMIZER,
% we must bail out, as there is no chance this model instantiates as an LP.
if forced_choice &  (ProblemClass.constraint.inequalities.semidefinite.linear | ProblemClass.constraint.inequalities.semidefinite.quadratic | ProblemClass.constraint.inequalities.semidefinite.polynomial | ProblemClass.constraint.inequalities.semidefinite.sigmonial)  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial |  solvers(i).constraint.inequalities.semidefinite.polynomial | solvers(i).constraint.inequalities.semidefinite.quadratic | solvers(i).constraint.inequalities.semidefinite.linear;
    end        
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'semidefinite constraints');
end
% Similarily, we have a SOCP by definition. We must support that
if forced_choice & ~socp_are_really_qc & (ProblemClass.constraint.inequalities.secondordercone.linear | ProblemClass.constraint.inequalities.secondordercone.nonlinear)  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.secondordercone.linear;
    end        
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'second-order cone constraints');
end

% ******************************************************
% Prune based on cone constraints
% ******************************************************
if ProblemClass.constraint.inequalities.secondordercone.linear & ~socp_are_really_qc & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.inequalities.secondordercone.linear | solvers(i).constraint.inequalities.semidefinite.linear | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'second-order cone constraints');
end  
if ProblemClass.constraint.inequalities.secondordercone.nonlinear & ~socp_are_really_qc & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.inequalities.secondordercone.nonlinear | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'nonlinear second-order cone constraints');
end  
if ProblemClass.constraint.inequalities.rotatedsecondordercone.linear & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.rotatedsecondordercone.linear | solvers(i).constraint.inequalities.secondordercone.linear | solvers(i).constraint.inequalities.semidefinite.linear;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'rotated second-order cone constraints');
end
if ProblemClass.constraint.inequalities.rotatedsecondordercone.nonlinear & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.rotatedsecondordercone.nonlinear | solvers(i).constraint.inequalities.secondordercone.linear | solvers(i).constraint.inequalities.semidefinite.linear;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'nonlinear rotated second-order cone constraints');
end
% if ProblemClass.constraint.inequalities.powercone && ~forced_choice  && length(solvers)>0
%     keep = ones(length(solvers),1);
%     for i = 1:length(solvers)                      
%         keep(i) = solvers(i).powercone;
%     end
%     [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'power cone constraints');
% end  

% ******************************************************
% Prune based on element-wise inequality constraints
% ******************************************************
if ProblemClass.constraint.inequalities.elementwise.sigmonial & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial;            
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'signomial constraints');
end  
if ProblemClass.constraint.inequalities.elementwise.polynomial & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'polynomial constraints');
end  
if ProblemClass.constraint.inequalities.elementwise.quadratic.nonconvex & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'nonconvex quadratic constraints');
end 
if ProblemClass.constraint.inequalities.elementwise.quadratic.convex | (ProblemClass.constraint.inequalities.secondordercone.linear & socp_are_really_qc) & ~forced_choice & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex | solvers(i).constraint.inequalities.elementwise.quadratic.convex | solvers(i).constraint.inequalities.secondordercone.linear | solvers(i).constraint.inequalities.semidefinite.linear;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'quadratic constraints');
end 
if ProblemClass.constraint.inequalities.elementwise.linear & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial | solvers(i).constraint.inequalities.semidefinite.quadratic | solvers(i).constraint.inequalities.elementwise.linear;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'constraints');
end  

% ******************************************************
% Prune based on element-wise constraints
% ******************************************************
if ProblemClass.constraint.equalities.sigmonial & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial | solvers(i).constraint.equalities.sigmonial;            
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'signomial equality constraints');
end  
if ProblemClass.constraint.equalities.polynomial & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)    
        indirect = solvers(i).constraint.inequalities.elementwise.sigmonial | solvers(i).constraint.inequalities.elementwise.polynomial;
        indirect = indirect | solvers(i).constraint.inequalities.elementwise.sigmonial | solvers(i).constraint.inequalities.elementwise.polynomial;
        direct = solvers(i).constraint.equalities.sigmonial | solvers(i).constraint.equalities.polynomial;
        keep(i) = direct | indirect;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'polynomial equality constraints');
end  
if ProblemClass.constraint.equalities.quadratic & ~forced_choice & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        indirect = solvers(i).constraint.inequalities.elementwise.sigmonial | solvers(i).constraint.inequalities.elementwise.polynomial | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex;
        direct = solvers(i).constraint.equalities.sigmonial |  solvers(i).constraint.equalities.polynomial |  solvers(i).constraint.equalities.quadratic;
        keep(i) = direct | indirect;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'quadratic equality constraints');
end 
if ProblemClass.constraint.equalities.linear & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers) 
        indirect = solvers(i).constraint.inequalities.elementwise.linear | solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial;
        direct = solvers(i).constraint.equalities.linear | solvers(i).constraint.equalities.sigmonial |  solvers(i).constraint.equalities.polynomial;
        keep(i) = direct | indirect;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'equality constraints');
end  

% ******************************************************
% Discrete data
% ******************************************************
if ProblemClass.constraint.integer & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.integer;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'integer variables');
end  
if ProblemClass.constraint.binary & ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.integer | solvers(i).constraint.binary;            
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'binary variables');
end  
if ProblemClass.constraint.sos1 & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         %keep(i) = solvers(i).constraint.integer | solvers(i).constraint.binary | solvers(i).constraint.sos2;            
         keep(i) = solvers(i).constraint.sos2;            
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'sos1 constraints');
end 
if ProblemClass.constraint.sos2 & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.integer | solvers(i).constraint.binary | solvers(i).constraint.sos2;                       
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'sos2 constraints');
end  
if ProblemClass.constraint.semicont  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         %keep(i) = solvers(i).constraint.integer | solvers(i).constraint.binary | solvers(i).constraint.sos2;            
         keep(i) = solvers(i).constraint.semivar;            
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'semicontinuous variables');
end  

% ******************************************************
% Equalities with multiple monomoials (rule out GP)
% ******************************************************
if ProblemClass.constraint.equalities.multiterm  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.equalities.multiterm;
    end    
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'multi-term monomial equalities');
end  
% FIXME
% No support for multiterm is YALMIPs current way of saying "GP solver". We
% use this flag to prune GPs based on objective too
if ~ProblemClass.gppossible  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.equalities.multiterm;
    end   
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'multi-term monomial equalities');
end  

% ******************************************************
% Complementarity constraints
% ******************************************************
if ProblemClass.constraint.complementarity.linear & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.complementarity.linear | solvers(i).constraint.integer | solvers(i).constraint.binary | solvers(i).constraint.equalities.polynomial | solvers(i).constraint.equalities.quadratic;            
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'complementarity constraints');
end  

% ******************************************************
% Interval data
% ******************************************************
if ProblemClass.interval  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).interval;
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'interval data');
end  

% ******************************************************
% Parametric problem
% ******************************************************
if ~forced_choice  & length(solvers)>0
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = (ProblemClass.parametric == solvers(i).parametric);
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'multiparametric models');
end

% ******************************************************
% Exponential cone representable (exp, log,...)
% ******************************************************
keep = ones(length(solvers),1);
if ~forced_choice  & length(solvers)>0
    for i = 1:length(solvers)
        keep(i) = (ProblemClass.exponentialcone <= solvers(i).exponentialcone) || (ProblemClass.exponentialcone <= solvers(i).evaluation);
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'exponential cones');
end
% ******************************************************
% Power cone representable 
% ******************************************************
keep = ones(length(solvers),1);
if ~forced_choice  & length(solvers)>0
    for i = 1:length(solvers)
        keep(i) = (ProblemClass.constraint.inequalities.powercone <= solvers(i).powercone) || (ProblemClass.constraint.inequalities.powercone <= solvers(i).evaluation);
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'exponential cones');
end

% ******************************************************
% General functions (sin, cos,...)
% ******************************************************
keep = ones(length(solvers),1);
if ~forced_choice  & length(solvers)>0
    for i = 1:length(solvers)
        keep(i) = (ProblemClass.evaluation <= solvers(i).evaluation) || (ProblemClass.exponentialcone && solvers(i).exponentialcone);
    end
    [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,'general nonlinearities');
end

% FIX : UUUUUUGLY
if isempty(solvers)
    solver = [];
else
    if length(options.solver)>0
        solver = [];

        % FIX : Re-use from above
        opsolver = lower(options.solver);
        splits = strfind(opsolver,',');
        if isempty(splits)
            names{1} = opsolver;
        else
            names = {};
            start = 1;
            for i  = 1:length(splits)
                names{i} = opsolver(start:splits(i)-1);
                start = splits(i)+1;
            end
            names{end+1} = opsolver(start:end);
        end

        temp = solvers;
        for i = 1:length(temp)
        if length(temp(i).version)>0
            temp(i).tag = lower([temp(i).tag '-' temp(i).version]);
        end
        end
    
        for i = 1:length(names)
            if isequal(names{i},'*')
                solver = solvers(1);
                break
            else
                j = find(strcmpi(lower({solvers.tag}),names{i}));
                if ~isempty(j)
                    solver = solvers(j(1));
                    break
                end
                j = find(strcmpi(lower({temp.tag}),names{i}));
                if ~isempty(j)
                    solver = solvers(j(1));
                    break
                end                
            end
        end
    else
        solver = solvers(1);
    end
end

if isempty(solver)
    if length(options.solver)>0 % User selected available solver, but it is not applicable
        problem = -4;
    else
        problem = -2;
    end
end

% FIX : Hack when chosing the wrong fmincon thingy
if ~isempty(solver)
    c1 = (length(options.solver)==0 | isequal(lower(options.solver),'fmincon')) & isequal(lower(solver.tag),'fmincon') & isequal(solver.version,'geometric');
    c2 = (length(options.solver)==0 | isequal(lower(options.solver),'snopt')) & isequal(lower(solver.tag),'snopt') & isequal(solver.version,'geometric');
    if c1 | c2
        if ~(ProblemClass.objective.sigmonial | ProblemClass.constraint.inequalities.elementwise.sigmonial)
            solver.version = 'standard';
            solver.call    = strrep(solver.call,'gp','');
            solver.objective.linear = 1;
            solver.objective.quadratic.convex = 1;
            solver.objective.quadratic.nonconvex = 1;
            solver.objective.polynomial = 1;
            solver.objective.sigmonial = 1;
            solver.constraint.equalities.elementwise.linear = 1;
            solver.constraint.equalities.elementwise.quadratic.convex = 1;
            solver.constraint.equalities.elementwise.quadratic.nonconvex = 1;
            solver.constraint.equalities.elementwise.polynomial = 1;
            solver.constraint.equalities.elementwise.sigmonial = 1;
            solver.constraint.inequalities.elementwise.linear = 1;
            solver.constraint.inequalities.elementwise.quadratic.convex = 1;
            solver.constraint.inequalities.elementwise.quadratic.nonconvex = 1;
            solver.constraint.inequalities.elementwise.polynomial = 1;
            solver.constraint.inequalities.elementwise.sigmonial = 1;
            solver.constraint.inequalities.secondordercone.linear = 1;
            solver.constraint.inequalities.semidefinite.linear = 1;
            solver.constraint.inequalities.semidefinite.quadratic = 1;
            solver.constraint.inequalities.semidefinite.polynomial = 1;
            solver.constraint.inequalities.semidefinite.nonlinear = 1;
            solver.exponentialcone = 1;
            solver.powercone = 1;
            solver.dual = 1;
            solver.evaluation = 1;
        end
    end
end

function temp = expandSolverName(temp)
for i = 1:length(temp)
    if length(temp(i).version)>0
        temp(i).tag = lower([temp(i).tag '-' temp(i).version]);
    end
end

function [solvers,failureMode] = pruneBasedOn(solvers,keep,failureMode,reason)

solvers = solvers(find(keep));
if ~any(keep)
    failureMode = reason;
end