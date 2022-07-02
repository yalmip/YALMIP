function [interfacedata,recoverdata,solver,diagnostic,F,Fremoved,ForiginalQuadratics] = compileinterfacedata(F,aux_obsolete,logdetStruct,h,options,findallsolvers,parametric)

persistent CACHED_SOLVERS
persistent allsolvers
persistent EXISTTIME
persistent NCHECKS

%% Initilize default empty outputs
diagnostic = [];
interfacedata = [];
recoverdata = [];
solver = [];
Fremoved  = [];
ForiginalQuadratics = [];

%% Did we make the call from SOLVEMP
if nargin<7
    parametric = 0;
end

%% Clean objective to default empty
if isa(h,'double')
    h = [];
end

% *************************************************************************
%% Exit if LOGDET objective is nonlinear
% *************************************************************************
if ~isempty(logdetStruct)
    for i = 1:length(logdetStruct.P)
        if ~is(logdetStruct.P{i},'linear')
            diagnostic.solvertime = 0;
            diagnostic.problem = -2;
            diagnostic.info = yalmiperror(diagnostic.problem,'');
            return
        end
    end
end

% *************************************************************************
%% EXTRACT LOW-RANK DESCRIPTION
% *************************************************************************
lowrankdetails = getlrdata(F);
if ~isempty(lowrankdetails)
    F = F(~is(F,'lowrank'));
end

% *************************************************************************
%% PERTURB STRICT INEQULAITIES
% *************************************************************************
if isa(options.shift,'sdpvar') | (options.shift~=0)
    F = shift(F,options.shift);
end

% *************************************************************************
%% ADD RADIUS CONSTRAINT
% *************************************************************************
if isa(options.radius,'sdpvar') | ~isinf(options.radius)
    x = recover(unique(union(depends(h),depends(F))));
    if length(x)>1
        F = F + (cone(x,options.radius));
    else
        F = F + (-options.radius <= x <= options.radius);
    end
    F = flatten(F);
end

% *************************************************************************
%% CONVERT LOGIC CONSTRAINTS
% *************************************************************************
[F,changed] = convertlogics(F);
if changed
    F = flatten(F);
    options.saveduals = 0; % Don't calculate duals since we changed the problem
end

% *************************************************************************
%% Take care of the nonlinear operators by converting expressions such as
% t = max(x,y) to standard conic models and mixed integer models
% This part also adds parts from logical expressions and mpower terms
% *************************************************************************
if options.expand
    % Experimental hack due to support for the PWQ function used for
    % quadratic dynamic programming with MPT.
    % FIX: Clean up and generalize
    try
        h1v = depends(h);
        h2v = getvariables(h);
        if ~isequal(h1v,h2v)            
            variables = uniquestripped([h1v h2v]);
        else
            variables = h1v;
        end
        extendedvariables = yalmip('extvariables');
        index_in_extended = find(ismembcYALMIP(variables,extendedvariables));
        if ~isempty(index_in_extended)
            extstruct = yalmip('extstruct',variables(index_in_extended));
            if ~isa(extstruct,'cell')
                extstruct = {extstruct};
            end
            for i = 1:length(extstruct)
                if isequal(extstruct{i}.fcn ,'pwq_yalmip')
                    [properties,Fz,arguments]=model(extstruct{i}.var,'integer',options,extstruct{i});
                    if iscell(properties)
                        properties = properties{1};
                    end
                    gain = getbasematrix(h,getvariables(extstruct{i}.var));
                    h = replace(h,extstruct{i}.var,0);
                    h = h + gain*properties.replacer;
                    F = F + Fz;
                end
            end
        end
    catch
    end
    [F,failure,cause,operators] = expandmodel(F,h,options);
    % Models temporarily added for analysis phase. Can be deleted now
    try
        F('placeholder') = [];
    catch
    end    
    % Some operators are temporarily modelled using linear hypograph
    % models, but might benefit from some other more direct form, and this
    % has been estabslished by the operator
    for i = 1:length(operators)
        if ~isempty(operators{i}.properties.replace)            
            F = replace(F, operators{i}.properties.replace.this,operators{i}.properties.replace.with);
            h = replace(h, operators{i}.properties.replace.this,operators{i}.properties.replace.with);
        end
    end
    F = flatten(F);
    if failure % Convexity propgation failed
        interfacedata = [];
        recoverdata = [];
        solver = '';
        diagnostic.solvertime = 0;
        diagnostic.problem = 14;
        diagnostic.info = yalmiperror(14,cause);
        return
    end
    evalVariables = unique(determineEvaluationBased(operators));
    if isempty(evalVariables)
        evaluation_based = 0;
        exponential_cone = 0;
    else
        used = [depends(h) depends(F)];
        usedevalvariables = intersect(used,evalVariables);
        evaluation_based = ~isempty(usedevalvariables);
        exponential_cone = isempty(setdiff(usedevalvariables,yalmip('expvariables')));
    end
else
    evalVariables = [];
    evaluation_based = 0;    
    exponential_cone = 0;
end

% *************************************************************************
%% LOOK FOR AVAILABLE SOLVERS
% Finding solvers can be very slow on some systems. To alleviate this
% problem, YALMIP can cache the list of available solvers.
% *************************************************************************
if (options.cachesolvers==0) | isempty(CACHED_SOLVERS)
    getsolvertime = clock;
    [solvers,kept,allsolvers] = getavailablesolvers(findallsolvers,options);
    getsolvertime = etime(clock,getsolvertime);
    % CODE TO INFORM USERS ABOUT SLOW NETWORKS!
    if isempty(EXISTTIME)
        EXISTTIME = getsolvertime;
        NCHECKS = 1;
    else
        EXISTTIME = [EXISTTIME getsolvertime];
        NCHECKS = NCHECKS + 1;
    end
    if (options.cachesolvers==0)
        if ((NCHECKS >= 3 & (sum(EXISTTIME)/NCHECKS > 1)) | EXISTTIME(end)>2)
            if warningon
                info = 'Warning: YALMIP has detected that your drive or network is unusually slow.\nThis causes a severe delay in OPTIMIZE when I try to find available solvers.\nTo avoid this, use the options CACHESOLVERS in SDPSETTINGS.\nSee the FAQ for more information.\n';
                fprintf(info);
            end
        end
    end
    if length(EXISTTIME) > 5
        EXISTTIME = EXISTTIME(end-4:end);
        NCHECKS = 5;
    end
    CACHED_SOLVERS = solvers;
else
    solvers = CACHED_SOLVERS;
end

% *************************************************************************
%% NO SOLVER AVAILABLE
% *************************************************************************
if isempty(solvers)
    diagnostic.solvertime = 0;
    if isempty(options.solver)
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
    else
        diagnostic.info = yalmiperror(-3,'YALMIP');
        diagnostic.problem = -3;
    end
    if warningon & options.warning & isempty(strfind(diagnostic.info,'No problems detected'))
        disp(['Warning: ' diagnostic.info]);
    end
    return
end

% *************************************************************************
%% CONVERT CONVEX QUADRATIC CONSTRAINTS
% We do not convert quadratic constraints to SOCPs if we have have
% sigmonial terms (thus indicating a GP problem), if we have relaxed
% nonlinear expressions, or if we have specified a nonlinear solver.
% Why do we convert them already here? Don't remember, should be cleaned up
% *************************************************************************
[monomtable,variabletype] = yalmip('monomtable');
F_vars = getvariables(F);
do_not_convert = any(variabletype(F_vars)==4);
% Handle case where optimizer adds a '+'
selected_solver = options.solver;
if ~isempty(selected_solver) && selected_solver(1) == '+'
    selected_solver(1)=[];
end
%do_not_convert = do_not_convert | ~solverCapable(solvers,options.solver,'constraint.inequalities.secondordercone');
do_not_convert = do_not_convert | strcmpi(selected_solver,'bmibnb');
do_not_convert = do_not_convert | strcmpi(selected_solver,'scip');
do_not_convert = do_not_convert | strcmpi(selected_solver,'snopt');
do_not_convert = do_not_convert | strcmpi(selected_solver,'knitro');
do_not_convert = do_not_convert | strcmpi(selected_solver,'knitro-standard');
do_not_convert = do_not_convert | strcmpi(selected_solver,'knitro-geometric');
do_not_convert = do_not_convert | strcmpi(selected_solver,'snopt-geometric'); 
do_not_convert = do_not_convert | strcmpi(selected_solver,'snopt-standard');
do_not_convert = do_not_convert | strcmpi(selected_solver,'bonmin');
do_not_convert = do_not_convert | strcmpi(selected_solver,'nomad');
do_not_convert = do_not_convert | strcmpi(selected_solver,'ipopt');
do_not_convert = do_not_convert | strcmpi(selected_solver,'ipopt-standard');
do_not_convert = do_not_convert | strcmpi(selected_solver,'ipopt-geometric');
do_not_convert = do_not_convert | strcmpi(selected_solver,'filtersd');
do_not_convert = do_not_convert | strcmpi(selected_solver,'filtersd-dense');
do_not_convert = do_not_convert | strcmpi(selected_solver,'filtersd-sparse');
do_not_convert = do_not_convert | strcmpi(selected_solver,'pennon');
do_not_convert = do_not_convert | strcmpi(selected_solver,'pennon-geometric');
do_not_convert = do_not_convert | strcmpi(selected_solver,'pennon-standard');
do_not_convert = do_not_convert | strcmpi(selected_solver,'pennlp');
do_not_convert = do_not_convert | strcmpi(selected_solver,'penbmi');
do_not_convert = do_not_convert | strcmpi(selected_solver,'fmincon');
do_not_convert = do_not_convert | strcmpi(selected_solver,'fmincon-standard');
do_not_convert = do_not_convert | strcmpi(selected_solver,'fmincon-geometric');
do_not_convert = do_not_convert | strcmpi(selected_solver,'sqplab');
do_not_convert = do_not_convert | strcmpi(selected_solver,'bmibnb');
do_not_convert = do_not_convert | strcmpi(selected_solver,'moment');
do_not_convert = do_not_convert | strcmpi(selected_solver,'sparsepop');
do_not_convert = do_not_convert | strcmpi(selected_solver,'baron');
do_not_convert = do_not_convert | strcmpi(selected_solver,'penlab');
do_not_convert = do_not_convert | strcmpi(selected_solver,'scip-nl');
do_not_convert = do_not_convert | (options.convertconvexquad == 0);
do_not_convert = do_not_convert | (options.relax == 1);
if ~do_not_convert & any(variabletype(F_vars))
    [F,socp_changed,infeasible,ForiginalQuadratics] = convertquadratics(F);
    if infeasible
        diagnostic.solvertime = 0;
        diagnostic.problem = 1;
        diagnostic.info = yalmiperror(diagnostic.problem,'YALMIP');
        return        
    end
    if socp_changed % changed holds the number of QC -> SOCC conversions
        options.saveduals = 0; % We cannot calculate duals since we changed the problem
        F_vars = []; % We have changed model so we cannot use this in categorizemodel
    end
else
    socp_changed = 0;
end

% CHEAT FOR QC
if socp_changed>0 & length(find(is(F,'socc')))==socp_changed
    socp_are_really_qc = 1;
else
    socp_are_really_qc = 0;
end

% *************************************************************************
%% WHAT KIND OF PROBLEM DO WE HAVE NOW?
% *************************************************************************
[ProblemClass,integer_variables,binary_variables,parametric_variables,uncertain_variables,semicont_variables,quad_info] = categorizeproblem(F,logdetStruct,h,options.relax,parametric,evaluation_based,F_vars,exponential_cone);

% Ugly fix to short-cut any decision on GP. min -x-y cannot be cast as GP,
% while min -x can, as we can invert the objective
ProblemClass.gppossible = 1;
if ~isempty(h)
    c = getbase(h);c = c(2:end);
    if nnz(c)>1
     if any(c<0)
         ProblemClass.gppossible = 0;
     end
    end
end
   
% *************************************************************************
%% SELECT SUITABLE SOLVER
% *************************************************************************
[solver,problem,~,failureMode] = selectsolver(options,ProblemClass,solvers,socp_are_really_qc,allsolvers);
if isempty(solver)
    diagnostic.solvertime = 0;
    if problem == -4 
        s = [options.solver ' does not support ' failureMode];
        diagnostic.info = yalmiperror(problem,s);
    elseif problem == -3 || problem == -9 
        s = options.solver;
        diagnostic.info = yalmiperror(problem,s);
    else
        diagnostic.info = yalmiperror(problem,'YALMIP');
    end
    diagnostic.problem = problem;

    if warningon & options.warning
        disp(['Warning: ' diagnostic.info]);
    end
    return
end
if length(solver.version)>0
    solver.tag = [solver.tag '-' solver.version];
end

if ProblemClass.constraint.complementarity.variable | ProblemClass.constraint.complementarity.linear | ProblemClass.constraint.complementarity.nonlinear
    if ~(solver.constraint.complementarity.variable | solver.constraint.complementarity.linear | solver.constraint.complementarity.nonlinear)               
        % Extract the terms in the complementarity constraints x^Ty==0,
        % x>=0, y>=0, since these involves bounds that should be appended
        % to the list of constraints from which we do bound propagation
        Fc = F(find(is(F,'complementarity')));      
        Ftemp = F;
        for i = 1:length(Fc)
            [Cx,Cy] = getComplementarityTerms(Fc(i));
            Ftemp = [Ftemp, Cx>=0, Cy >=0];
        end
        % FIXME: SYNC with expandmodel       
        setupBounds(Ftemp,options,extendedvariables);
                        
        [F] = modelComplementarityConstraints(F,solver,ProblemClass);  
        % FIXME Reclassify should be possible to do manually!
        oldProblemClass = ProblemClass;
        [ProblemClass,integer_variables,binary_variables,parametric_variables,uncertain_variables,semicont_variables,quad_info] = categorizeproblem(F,logdetStruct,h,options.relax,parametric,evaluation_based,F_vars,exponential_cone);
        ProblemClass.gppossible = oldProblemClass.gppossible;
    elseif solver.constraint.complementarity.variable
        % Solver supports x(i)*x(j)==0
        Fok = [];
        Freform = [];
        Fc = F(find(is(F,'complementarity')));              
        for i = 1:length(Fc)
            [Cx,Cy] = getComplementarityTerms(Fc(i));
            if (islinear(Cx) & islinear(Cy) & is(Cx,'lpcone') &  is(Cy,'lpcone'))
             Fok = [Fok, Fc(i)];
            else
                s1 = sdpvar(length(Cx),1);
                s2 = sdpvar(length(Cy),1);                
                Freform = [Freform,complements(s1>=0,s2>=0),s1 == Cx, s2 == Cy];
            end            
        end   
        F = F-Fc;
        F = F + Fok + Freform;
    end
end

% *************************************************************************
%% DID WE SELECT THE INTERNAL BNB SOLVER
% IN THAT CASE, SELECT LOCAL SOLVER
% (UNLESS ALREADY SPECIFIED IN OPTIONS.BNB)
% *************************************************************************
localsolver.qc = 0;
localsolver = solver;
if strcmpi(solver.tag,'bnb')
    [solver,diagnostic] = setupBNB(solver,ProblemClass,options,solvers,socp_are_really_qc,F,h,logdetStruct,parametric,evaluation_based,F_vars,exponential_cone,allsolvers);
    if ~isempty(diagnostic)
        return
    end 
end

if strfind(lower(solver.tag),'sparsecolo')
    temp_options = options;
    temp_options.options.forceglobal = 0;
    temp_options.solver = options.sparsecolo.SDPsolver;
    tempProblemClass = ProblemClass;   
    localsolver = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
    if isempty(localsolver) | strcmpi(localsolver.tag,'sparsecolo')
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
        return
    end
    solver.sdpsolver = localsolver;
end

if strfind(lower(solver.tag),'frlib')
    temp_options = options;
    temp_options.forceglobal = 0;
    temp_options.solver = options.frlib.solver;
    tempProblemClass = ProblemClass;   
    localsolver = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
    if isempty(localsolver) | strcmpi(localsolver.tag,'frlib')
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
        return
    end
    solver.solver = localsolver;
end

% *************************************************************************
%% DID WE SELECT THE MPCVX SOLVER
% IN THAT CASE, SELECT SOLVER TO SOLVE BOUND COMPUTATIONS
% *************************************************************************
localsolver.qc = 0;
localsolver = solver;
if strcmpi(solver.tag,'mpcvx')
    temp_options = options;
    temp_options.solver = options.mpcvx.solver;
    tempProblemClass = ProblemClass;    
    tempProblemClass.objective.quadratic.convex = tempProblemClass.objective.quadratic.convex | tempProblemClass.objective.quadratic.nonconvex;
    tempProblemClass.objective.quadratic.nonconvex = 0;
    tempProblemClass.parametric = 0;
    localsolver = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
    if isempty(localsolver) | strcmpi(localsolver.tag,'bnb') | strcmpi(localsolver.tag,'kktqp')
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
        return
    end
    solver.lower = localsolver;
end

% *************************************************************************
%% DID WE SELECT THE INTERNAL EXPERIMENTAL KKT SOLVER
% IN THAT CASE, SELECT SOLVER TO SOLVE THE MILP PROBLEM
% *************************************************************************
localsolver.qc = 0;
localsolver = solver;
if strcmpi(solver.tag,'kktqp')
    temp_options = options;
    temp_options.forceglobal = 0;
    temp_options.solver = '';
    tempProblemClass = ProblemClass;
    tempProblemClass.constraint.binary = 1;
    tempProblemClass.objective.quadratic.convex = 0;
    tempProblemClass.objective.quadratic.nonconvex = 0;
    localsolver = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
    if isempty(localsolver) | strcmpi(localsolver.tag,'bnb') | strcmpi(localsolver.tag,'kktqp')
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
        return
    end
    solver.lower = localsolver;
end

% *************************************************************************
%% DID WE SELECT THE LMIRANK?
% FIND SDP SOLVER FOR INITIAL SOLUTION
% *************************************************************************
if strcmpi(solver.tag,'lmirank')
    temp_options = options;
    temp_options.forceglobal = 0;
    temp_options.solver = options.lmirank.solver;
    tempProblemClass = ProblemClass;
    tempProblemClass.constraint.inequalities.rank = 0;
    tempProblemClass.constraint.inequalities.semidefinite.linear = 1;
    tempProblemClass.objective.linear = 1;
    initialsolver = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
    if isempty(initialsolver) | strcmpi(initialsolver.tag,'lmirank')
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
        return
    end
    solver.initialsolver = initialsolver;
end

% *************************************************************************
%% DID WE SELECT THE VSDP SOLVER? Define a solver for VSDP to use
% *************************************************************************
if strfind(solver.tag,'VSDP')
    temp_options = options;
    temp_options.forceglobal = 0;
    temp_options.solver = options.vsdp.solver;
    tempProblemClass = ProblemClass;
    tempProblemClass.interval = 0;
    tempProblemClass.constraint.inequalities.semidefinite.linear =  tempProblemClass.constraint.inequalities.semidefinite.linear | tempProblemClass.objective.quadratic.convex;
    tempProblemClass.constraint.inequalities.semidefinite.linear =  tempProblemClass.constraint.inequalities.semidefinite.linear | tempProblemClass.constraint.inequalities.secondordercone;
    tempProblemClass.constraint.inequalities.secondordercone = 0;
    tempProblemClass.objective.quadratic.convex = 0;
    initialsolver = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
    if isempty(initialsolver) | strcmpi(initialsolver.tag,'vsdp')
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
        return
    end
    solver.solver = initialsolver;
end

% *************************************************************************
%% DID WE SELECT THE INTERNAL BMIBNB SOLVER? SELECT UPPER/LOWER SOLVERs
% (UNLESS ALREADY SPECIFIED IN OPTIONS)
% *************************************************************************
if strcmpi(solver.tag,'bmibnb')
    [solver,diagnostic] = setupBMIBNB(solver,ProblemClass,options,solvers,socp_are_really_qc,F,h,logdetStruct,parametric,evaluation_based,F_vars,exponential_cone,allsolvers);
    if ~isempty(diagnostic)
        return
    end   
end

% *************************************************************************
%% DID WE SELECT THE INTERNAL SDPMILP SOLVER
% This solver solves MISDP problems by solving MILP problems and adding SDP
% cuts based on the infasible MILP solution.
% *************************************************************************
if strcmpi(solver.tag,'cutsdp')

    % Relax problem for lower solver
    tempProblemClass = ProblemClass;
    tempProblemClass.constraint.inequalities.elementwise.linear =  tempProblemClass.constraint.inequalities.elementwise.linear |     tempProblemClass.constraint.inequalities.semidefinite.linear | tempProblemClass.constraint.inequalities.secondordercone.linear;
    tempProblemClass.constraint.inequalities.semidefinite.linear = 0;
    tempProblemClass.constraint.inequalities.secondordercone.linear = 0;
    tempProblemClass.objective.quadratic.convex = 0;
    
    temp_options = options;
    temp_options.solver = options.cutsdp.solver;

    if strcmp(options.cutsdp.solver,'bnb')
        error('BNB can not be used in CUTSDP. Please install and use a better MILP solver');
    end
    
    [lowersolver,problem] = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);

    if ~isempty(lowersolver) & strcmpi(lowersolver.tag,'bnb')
        error('BNB can not be used in CUTSDP. Please install and use a better MILP solver');
    end
        
    if isempty(lowersolver) | strcmpi(lowersolver.tag,'cutsdp') |strcmpi(lowersolver.tag,'bmibnb') | strcmpi(lowersolver.tag,'bnb')
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
        return
    end
    solver.lower = lowersolver;
end

showprogress(['Solver chosen : ' solver.tag],options.showprogress);
 
% *************************************************************************
%% CONVERT SOS2 to binary constraints for solver not supporting sos2
% *************************************************************************
if  ProblemClass.constraint.sos2 & ~solver.constraint.sos2
    [F,binary_variables] = expandsos2(F,binary_variables);
end

% *************************************************************************
%% CONVERT MAXDET TO SDP USING GEOMEAN?
% *************************************************************************
% MAXDET using geometric mean construction
if ~isempty(logdetStruct)
    if isequal(solver.tag,'BNB')
        can_solve_maxdet = solver.lower.objective.maxdet.convex;
        can_solve_expcone = solver.lower.exponentialcone;
    else
        can_solve_maxdet = solver.objective.maxdet.convex;
        can_solve_expcone = solver.exponentialcone;
    end
    if ~can_solve_maxdet
        if isempty(h)
            h = 0;
        end
        if 0%can_solve_expcone
            for i = 1:length(logdetStruct.P)
                [vi,Modeli] = eigv(logdetStruct.P{i});
                F = [F, Modeli, logdetStruct.P{i} >= 0];
                log_vi = log(vi);
                h = h + logdetStruct.gain(i)*sum(log_vi);
                evalVariables = union(evalVariables,getvariables( log_vi));
            end
        else
        t = sdpvar(1,1);
        Ptemp = [];
        for i = 1:length(logdetStruct.P)
            Ptemp = blkdiag(Ptemp,logdetStruct.P{i});
        end
        P = {Ptemp};
        if length(F)>0
            if isequal(P,sdpvar(F(end)))
                F = F(1:end-1);
            end
        end
        F = F + detset(t,P{1});
        if isempty(h)
            h = -t;
            if length(logdetStruct.P) > 1 && options.verbose>0 && options.warning>0
                disp(' ')
                disp('Objective -sum logdet(P_i) has been changed to -sum det(P_i)^(1/(2^ceil(log2(length(P_i))))).')
                disp('This is not an equivalent transformation. You should use SDPT3 which supports MAXDET terms')
                disp('See the MAXDET section in the manual for details.')
                disp(' ')
            end
        else
            h = h-t;
            % Warn about logdet -> det^1/m
            if options.verbose>0 & options.warning>0
                disp(' ')
                disp('Objective c''x-sum logdet(P_i) has been changed to c''x-sum det(P_i)^(1/(2^ceil(log2(length(P_i))))).')
                disp('This is not an equivalent transformation. You should use SDPT3 which supports MAXDET terms')
                disp('See the MAXDET section in the manual for details.')
                disp(' ')
            end
        end
        end
        P = [];
        logdetStruct = [];
    end
end

% *************************************************************************
%% Change binary variables to integer?
% *************************************************************************
old_binary_variables = binary_variables;
if ~isempty(binary_variables) & (solver.constraint.binary==0)
    x_bin = recover(binary_variables(ismember(binary_variables,unique([getvariables(h) getvariables(F)]))));
    F = F + (x_bin<=1)+(x_bin>=0);
    integer_variables = union(binary_variables,integer_variables);
    binary_variables = [];
end

% *************************************************************************
%% Model quadratics using SOCP?
% Should not be done when using PENNLP or BMIBNB or FMINCON, or if we have relaxed the
% monmial terms or...Ouch, need to clean up all special cases, this sucks.
% *************************************************************************
convertQuadraticObjective = ~strcmpi(solver.tag,'pennlp-standard');
convertQuadraticObjective = convertQuadraticObjective & ~strcmpi(solver.tag,'bmibnb');
relaxed = (options.relax==1 | options.relax==3);
%| (~isempty(quad_info) & strcmp(solver.tag,'bnb') & localsolver.objective.quadratic.convex==0)
convertQuadraticObjective = convertQuadraticObjective & (~relaxed & (~isempty(quad_info) & solver.objective.quadratic.convex==0));
%convertQuadraticObjective = convertQuadraticObjective; % | strcmpi(solver.tag,'cutsdp');
if any(strcmpi(solver.tag,{'bnb','cutsdp'})) & ~isempty(quad_info)
    if solver.lower.objective.quadratic.convex==0
        convertQuadraticObjective = 1;
    end
end

if convertQuadraticObjective
    t = sdpvar(1,1);
    x = quad_info.x;
    R = quad_info.R;
    if ~isempty(R)
        c = quad_info.c;
        f = quad_info.f;
        F = F + lmi(cone([2*R*x;1-(t-f)],1+t-f));
        h = t+c'*x;
        if options.warmstart
            xx = value(x);
            ff = norm(quad_info.R*xx)^2+f;
            if ~isnan(ff)
                assign(t,ff);
            end
        end
    end
    quad_info = [];
end
if solver.constraint.inequalities.rotatedsecondordercone.linear == 0
    [F,changed] = convertlorentz(F);
    if changed
        options.saveduals = 0; % We cannot calculate duals since we change the problem
    end
end
% Whoa, horrible tests to find out when to convert SOCP to SDP
% This should not be done if :
%   1. Problem is actually a QCQP and solver supports this
%   2. Problem is integer, local solver supports SOCC
%   3. Solver supports SOCC
if ~((solver.constraint.inequalities.elementwise.quadratic.convex == 1) & socp_are_really_qc)
    if ~(strcmp(solver.tag,'bnb') & socp_are_really_qc & localsolver.constraint.inequalities.elementwise.quadratic.convex==1 )
        if ((solver.constraint.inequalities.secondordercone.linear == 0) | (strcmpi(solver.tag,'bnb') & localsolver.constraint.inequalities.secondordercone.linear==0))
            if solver.constraint.inequalities.semidefinite.linear
                [F,changed] = convertsocp(F);
            else
                [F,changed] = convertsocp2NONLINEAR(F);
            end
            if changed
                options.saveduals = 0; % We cannot calculate duals since we change the problem
            end
        end
    end
end

% *************************************************************************
%% Add logaritmic barrier cost/constraint for MAXDET and SDPT3-4. Note we
% have to add it her in order for a complex valued matrix to be converted.
% *************************************************************************
if ~isempty(logdetStruct) & solver.objective.maxdet.convex==1 & solver.constraint.inequalities.semidefinite.linear
    for i = 1:length(logdetStruct.P)
        F = F + (logdetStruct.P{i} >= 0);
        if ~isreal(logdetStruct.P{i})
            logdetStruct.gain(i) = logdetStruct.gain(i)/2;
            ProblemClass.complex = 1;
        end
    end
end

if ((solver.complex==0) & ProblemClass.complex) | ((strcmp(solver.tag,'bnb') & localsolver.complex==0) & ProblemClass.complex)
    showprogress('Converting to real constraints',options.showprogress)
    F = imag2reallmi(F);
    if ~isempty(logdetStruct) 
        for i = 1:length(logdetStruct.P)
            P{i} = sdpvar(F(end-length(logdetStruct.P)+i));
        end
    end
    options.saveduals = 0; % We cannot calculate duals since we change the problem
%else
%    complex_logdet = zeros(length(P),1);
end

% *************************************************************************
%% CREATE OBJECTIVE FUNCTION c'*x+x'*Q*x
% *************************************************************************
showprogress('Processing objective function',options.showprogress);
try
    % If these solvers, the Q term is placed in c, hence quadratic terms
    % are treated as any other nonlinear term
    geometric = strcmpi(solver.tag,'fmincon-geometric')| strcmpi(solver.tag,'gpposy') | strcmpi(solver.tag,'mosek-geometric') | strcmpi(solver.tag,'snopt-geometric') | strcmpi(solver.tag,'ipopt-geometric') | strcmpi(solver.tag,'pennon-geometric');
    if strcmpi(solver.tag,'bnb')
        lowersolver = lower([solver.lower.tag '-' solver.lower.version]);
        if strcmpi(lowersolver,'fmincon-geometric')| strcmpi(lowersolver,'gpposy-') |  strcmpi(lowersolver,'mosek-geometric')
            geometric = 1;
        end
    end
    if strcmpi(solver.tag,'bmibnb') | strcmpi(solver.tag,'sparsepop') | strcmpi(solver.tag,'pennlp-standard') | geometric | evaluation_based ;
        tempoptions = options;
        tempoptions.relax = 1;
        [c,Q,f]=createobjective(h,logdetStruct,tempoptions,quad_info);
    else
        [c,Q,f]=createobjective(h,logdetStruct,options,quad_info);
    end
catch
    error(lasterr)
end

% *************************************************************************
%% Convert {F(x),G(x)} to a numerical SeDuMi-like format
% *************************************************************************
showprogress('Processing constraints',options.showprogress);
F = lmi(F);
[F_struc,K,KCut,schur_funs,schur_data,schur_variables] = lmi2sedumistruct(F);
% We add a field to remember the dimension of the logarithmic cost.
% Actually, the actually value is not interesting, we know that the
% logarithmic cost corresponds to the last LP or SDP constraint anyway
if isempty(logdetStruct)
    K.m = 0;
else
    for i = 1:length(logdetStruct.P)
        K.m(i) = length(logdetStruct.P{i});        
    end
    K.maxdetgain = logdetStruct.gain;
end

if ~isempty(schur_funs)
    if length(schur_funs)<length(K.s)
        schur_funs{length(K.s)}=[];
        schur_data{length(K.s)}=[];
        schur_variables{length(K.s)}=[];
    end
end
K.schur_funs = schur_funs;
K.schur_data = schur_data;
K.schur_variables = schur_variables;

% *************************************************************************
%% SOME HORRIBLE CODE TO DETERMINE USED VARIABLES
% *************************************************************************
% Which sdpvar variables are actually in the problem
used_variables_LMI = find(any(F_struc(:,2:end),1));
used_variables_obj = find(any(c',1) | any(Q));
if isequal(used_variables_LMI,used_variables_obj)
    used_variables = used_variables_LMI;
else
    used_variables = uniquestripped([used_variables_LMI used_variables_obj]);
end
if ~isempty(K.sos)
    for i = 1:length(K.sos.type)
        used_variables = uniquestripped([used_variables K.sos.variables{i}(:)']);
    end
end
% The problem is that linear terms might be missing in problems with only
% nonlinear expressions
[monomtable,variabletype] = yalmip('monomtable');
if (options.relax==1)|(options.relax==3)
    monomtable = [];
    nonlinearvariables = [];
    linearvariables = used_variables;
else
    nonlinearvariables = find(variabletype);
    linearvariables = used_variables(find(variabletype(used_variables)==0));
end
needednonlinear = nonlinearvariables(ismembcYALMIP(nonlinearvariables,used_variables));
linearinnonlinear = find(sum(abs(monomtable(needednonlinear,:)),1));
missinglinear = setdiff(linearinnonlinear(:),linearvariables);
used_variables = uniquestripped([used_variables(:);missinglinear(:)]);

% *************************************************************************
%% So are we done now? No... What about variables hiding inside so called
% evaluation variables. We detect these, and at the same time set up the
% structures needed to support general functions such as exp, log, etc
% NOTE : This is experimental code
% FIX  : Clean up...
% *************************************************************************
[evalMap,evalVariables,used_variables,nonlinearvariables,linearvariables] = detectHiddenNonlinear(used_variables,options,nonlinearvariables,linearvariables,evalVariables);

% Attach information on the evaluation based variables that was generated
% when the model was expanded
if ~isempty(evalMap)
    for i = 1:length(operators)
        index = find(operators{i}.properties.models(1) == used_variables(evalVariables));
        if ~isempty(index)
            evalMap{index}.properties = operators{i}.properties;
        end
    end
    for i = 1:length(evalMap)
        for j = 1:length(evalMap{i}.computes)
            evalMap{i}.computes(j) = find(evalMap{i}.computes(j) == used_variables);
        end
    end
    
    % Add information about which argument is the variable
    for i = 1:length(evalMap)
        for j = 1:length(evalMap{i}.arg)-1
            if isa(evalMap{i}.arg{j},'sdpvar')
                evalMap{i}.argumentIndex = j;
                break
            end
        end
    end
end

% *************************************************************************
%% REMOVE UNNECESSARY VARIABLES FROM PROBLEM
% *************************************************************************
if length(used_variables)<yalmip('nvars')
    c = c(used_variables);
    if 0
        % very slow in some extreme cases
        Q = Q(:,used_variables);Q = Q(used_variables,:);
    else
        [i,j,s] = find(Q);
        keep = ismembcYALMIP(i,used_variables) & ismembcYALMIP(j,used_variables);
        i = i(keep);
        j = j(keep);
        s = s(keep);
        [ii,jj] = ismember(1:length(Q),used_variables);
        i = jj(i);
        j = jj(j);
        Q = sparse(i,j,s,length(used_variables),length(used_variables));
    end
                        
    if ~isempty(F_struc)
        F_struc = sparse(F_struc(:,[1 1+used_variables]));
    elseif size(F_struc,1) == 0
        F_struc = [];
    end
end

% *************************************************************************
%% Map variables and constraints in low-rank definition to local stuff
% *************************************************************************
if ~isempty(lowrankdetails)
    % Identifiers of the SDP constraints
    lmiid = getlmiid(F);
    for i = 1:length(lowrankdetails)
        lowrankdetails{i}.id = find(ismember(lmiid,lowrankdetails{i}.id));
        if ~isempty(lowrankdetails{i}.variables)
            index = ismember(used_variables,lowrankdetails{i}.variables);
            lowrankdetails{i}.variables = find(index);
        end
    end
end

% *************************************************************************
%% SPECIAL VARIABLES
% Relax = 1 : relax both integers and nonlinear stuff
% Relax = 2 : relax integers
% Relax = 3 : relax nonlinear stuff
% *************************************************************************
if (options.relax==1) | (options.relax==3)
    nonlins = [];
end
if (options.relax == 1) | (options.relax==2)
    integer_variables = [];
    binary_variables  = [];
    semicont_variables  = [];
    old_binary_variables  = find(ismember(used_variables,old_binary_variables));
else
    integer_variables = find(ismember(used_variables,integer_variables));
    binary_variables  = find(ismember(used_variables,binary_variables));
    semicont_variables = find(ismember(used_variables,semicont_variables));
    old_binary_variables  = find(ismember(used_variables,old_binary_variables));
end
parametric_variables  = find(ismember(used_variables,parametric_variables));
extended_variables =  find(ismember(used_variables,yalmip('extvariables')));
aux_variables =  find(ismember(used_variables,yalmip('auxvariables')));
if ~isempty(K.sos)
    for i = 1:length(K.sos.type)
        K.sos.variables{i} =  find(ismember(used_variables,K.sos.variables{i}));
        K.sos.variables{i} = K.sos.variables{i}(:); 
    end
end

% *************************************************************************
%% Equality constraints not supported or supposed to be removed
% *************************************************************************
% We may save some data in order to reconstruct
% dual variables related to equality constraints that
% have been removed.
oldF_struc = [];
oldQ = [];
oldc = [];
oldK = K;
Fremoved = [];
if (K.f>0)
    if (isequal(solver.tag,'BNB') && ~solver.lower.constraint.equalities.linear) ||  (isequal(solver.tag,'BMIBNB') && ~solver.lowersolver.constraint.equalities.linear)
        badLower = 1;
    else
        badLower = 0;
    end
    % reduce if user explicitely says remove, or user says nothing but
    % solverdefinitions does, and there are no nonlinear variables
    if ~badLower && ((options.removeequalities==1 | options.removeequalities==2) & isempty(intersect(used_variables,nonlinearvariables))) | ((options.removeequalities==0) & (solver.constraint.equalities.linear==-1))
        showprogress('Solving equalities',options.showprogress);
        [x_equ,H,A_equ,b_equ,factors] = solveequalities(F_struc,K,options.removeequalities==1);
        % Exit if no consistent solution exist
        if (norm(A_equ*x_equ-b_equ,'inf')>1e-5)%sqrt(eps)*size(A_equ,2))
            diagnostic.solvertime = 0;
            diagnostic.info = yalmiperror(1,'YALMIP');
            diagnostic.problem = 1;
            solution = diagnostic;
            solution.variables = used_variables(:);
            solution.optvar = x_equ;
            % And we are done! Save the result
            % sdpvar('setSolution',solution);
            return
        end
        % We dont need the rows for equalities anymore
        oldF_struc = F_struc;
        oldc = c;
        oldQ = Q;
        oldK = K;
        F_struc = F_struc(K.f+1:end,:);
        K.f = 0;
        Fold = F;
        [nlmi neq]=sizeOLD(F);
        iseq = is(Fold(1:(nlmi+neq)),'equality');
        F = Fold(find(~iseq));
        Fremoved = Fold(find(iseq));

        % No variables left. Problem solved!
        if size(H,2)==0
            diagnostic.solvertime = 0;
            diagnostic.info = yalmiperror(0,'YALMIP');
            diagnostic.problem = 0;
            solution = diagnostic;
            solution.variables = used_variables(:);
            solution.optvar = x_equ;
            % And we are done! Save the result
            % Note, no dual is saved
            yalmip('setSolution',solution);
            p = check(F);
            if any(p<1e-5)
                diagnostic.info = yalmiperror(1,'YALMIP');
                diagnostic.problem = 1;
            end
            return
        end
        showprogress('Converting problem to new basis',options.showprogress)

        % objective in new basis
        f = f + x_equ'*Q*x_equ;
        c = H'*c + 2*H'*Q*x_equ;
        Q = H'*Q*H;Q=((Q+Q')/2);
        % LMI in new basis
        F_struc = [F_struc*[1;x_equ] F_struc(:,2:end)*H];
    else      
        % Solver does not support equality constraints and user specifies
        % double-sided inequalitis to remove them, or solver is used from
        % lmilab or similiar solver
        if (solver.constraint.equalities.linear==0 | options.removeequalities==-1 | badLower)
            % Add equalities
            F_struc = [-F_struc(1:1:K.f,:);F_struc];
            K.l = K.l+K.f*2;
            % Keep this in mind...
            K.fold = K.f;
            K.f = 0;
        end
        % For simpliciy we introduce a dummy coordinate change
        x_equ   = 0;
        H       = 1;
        factors = [];
    end
else
    x_equ   = 0;
    H       = 1;
    factors = [];
end


% *************************************************************************
%% Setup the initial solution
% *************************************************************************
x0 = [];
if options.warmstart
    if solver.supportsinitial == 0
        error(['You have specified an initial point, but the selected solver (' solver.tag ') does not support warm-starts through YALMIP']);
    end
    if options.relax
        x0_used = relaxdouble(recover(used_variables));
    else
        %FIX : Do directly using yalmip('solution')
        %solution = yalmip('getsolution');
        x0_used = double(recover(used_variables));
    end
    x0 = zeros(sdpvar('nvars'),1);
    x0(used_variables)  = x0_used(:);
    if ~solver.supportsinitialNAN
        x0(isnan(x0))=0;
    end
end
if ~isempty(x0)
    % Get a coordinate in the reduced space
    x0 = H\(x0(used_variables)-x_equ);
end

% Monomial table for nonlinear variables
% FIX : Why here!!! mt handled above also
[mt,variabletype] = yalmip('monomtable');
if size(mt,1)>size(mt,2)
    mt(size(mt,1),size(mt,1)) = 0;
end
% In local variables
mt = mt(used_variables,used_variables);
variabletype = variabletype(used_variables);
if (options.relax == 1)|(options.relax==3)
    mt = speye(length(used_variables));
    variabletype = variabletype*0;
end

% FIX : Make sure these things work...
lub = yalmip('getbounds',used_variables);
lb = lub(:,1)-inf;
ub = lub(:,2)+inf;
lb(old_binary_variables) = max(lb(old_binary_variables),0);
ub(old_binary_variables) = min(ub(old_binary_variables),1);

% This does not work if we have used removeequalities, so we clear them for
% safety. note that bounds are not guaranteed to be used according to the
% manual, so this is allowed, although it might be a bit inconsistent to
% some users.
if ~isempty(oldc)
    lb = [];
    ub = [];
end

% Sanity check
if (~isempty(c) && any(isnan(c) )) || (~isempty(Q) && any(any(isnan(Q))))
    disp('You have NaNs in objective (<a href="yalmip.github.io/naninmodel">learn to debug</a>)')
    error('NaN in objective')
end
if (~isempty(ub) && any(isnan(ub))) || (~isempty(lb) && any(isnan(lb)))
    disp('You have NaNs in a bound (<a href="yalmip.github.io/naninmodel">learn to debug</a>)')
    error('NaN in bounds')
end
if ~isempty(F_struc)
    if any(any(isnan(F_struc)))
        disp('You have NaNs in your constraints (<a href="yalmip.github.io/naninmodel">learn to debug</a>)')
        error('NaN in model')
    end
end

% *************************************************************************
%% Does the solver support high precision input? If not, make sure the 
% input is only double precision.
% *************************************************************************

if solver.supportshighprec ~= 1
    if isa(F_struc, 'gem') || isa(F_struc, 'sgem')
        F_struc = double(F_struc);
    end
    if isa(c, 'gem') || isa(c, 'sgem')
        c = double(c);
    end
    if isa(Q, 'gem') || isa(Q, 'sgem')
        Q = double(Q);
    end
    if isa(f, 'gem') || isa(f, 'sgem')
        f = double(f);
    end
    if isa(lb, 'gem') || isa(lb, 'sgem')
        lb = double(lb);
    end
    if isa(ub, 'gem') || isa(ub, 'sgem')
        ub = double(ub);
    end
    if isa(x0, 'gem') || isa(x0, 'sgem')
        x0 = double(x0);
    end
    if isa(oldF_struc, 'gem') || isa(oldF_struc, 'sgem')
        oldF_struc = double(oldF_struc);
    end
    if isa(oldc, 'gem') || isa(oldc, 'sgem')
        oldc = double(oldc);
    end
end


% *************************************************************************
%% GENERAL DATA EXCHANGE WITH SOLVER
% *************************************************************************
interfacedata.F_struc = F_struc;
interfacedata.c = c;
interfacedata.Q = Q;
interfacedata.f = f;
interfacedata.K = K;
interfacedata.lb = lb;
interfacedata.ub = ub;
interfacedata.x0 = x0;
interfacedata.options = options;
interfacedata.solver  = solver;
interfacedata.monomtable = mt;
interfacedata.variabletype = variabletype;
interfacedata.integer_variables   = integer_variables;
interfacedata.binary_variables    = binary_variables;
interfacedata.semicont_variables    = semicont_variables;
interfacedata.implied_integer = [];
interfacedata.semibounds = [];
interfacedata.uncertain_variables = [];
interfacedata.parametric_variables= parametric_variables;
interfacedata.extended_variables  = extended_variables;
interfacedata.aux_variables  = aux_variables;
interfacedata.used_variables      = used_variables;
interfacedata.lowrankdetails = lowrankdetails;
interfacedata.problemclass = ProblemClass;
interfacedata.KCut = KCut;
interfacedata.getsolvertime = 1;
% Data to be able to recover duals when model is reduced
interfacedata.oldF_struc = oldF_struc;
interfacedata.oldc = oldc;
interfacedata.oldK = oldK;
interfacedata.factors = factors;
interfacedata.Fremoved = Fremoved;
interfacedata.evalMap = evalMap;
interfacedata.evalVariables = evalVariables;
interfacedata.evaluation_scheme = [];
interfacedata.lift = [];
if strcmpi(solver.tag,'bmibnb')
    interfacedata.equalitypresolved = 0;
    interfacedata.presolveequalities = 1;
else
    interfacedata.equalitypresolved = 1;
    interfacedata.presolveequalities = 1;
end
interfacedata.ProblemClass = ProblemClass;
interfacedata.dualized = is(F,'dualized');
interfacedata.presolved = 0;
if interfacedata.options.usex0==1
    interfacedata.options.warmstart=1;
end

% *************************************************************************
%% GENERAL DATA EXCANGE TO RECOVER SOLUTION AND UPDATE YALMIP VARIABLES
% *************************************************************************
recoverdata.H = H;
recoverdata.x_equ = x_equ;
recoverdata.used_variables = used_variables;

%%
function yesno = warningon

s = warning;
if isa(s,'char')
    yesno = isequal(s,'on');
else
    yesno = isequal(s(1).state,'on');
end

%%
function [evalMap,evalVariables,used_variables,nonlinearvariables,linearvariables] = detectHiddenNonlinear(used_variables,options,nonlinearvariables,linearvariables,eIN)

%evalVariables = yalmip('evalVariables');
evalVariables = eIN;
old_used_variables = used_variables;
goon = 1;
if ~isempty(evalVariables)
    while goon
        % Which used_variables are representing general functions
     %   evalVariables = yalmip('evalVariables');
     evalVariables = eIN;
        usedEvalVariables = find(ismember(used_variables,evalVariables));
        evalMap =  yalmip('extstruct',used_variables(usedEvalVariables));
        if ~isa(evalMap,'cell')
            evalMap = {evalMap};
        end
        % Find all variables used in the arguments of these functions
        hidden = [];
        for i = 1:length(evalMap)
            % Find main main argument (typically first argument, but this
            % could be different in a user-specified sdpfun object) 
            for aux = 1:length(evalMap{i}.arg)-1
                if isa(evalMap{i}.arg{aux},'sdpvar')
                    X = evalMap{i}.arg{aux};
                    break
                end
            end
            n = length(X);
            if isequal(getbase(X),[spalloc(n,1,0) speye(n)])% & is(evalMap{i}.arg{1},'linear')
                for j = 1:length(evalMap{i}.arg)-1
                    % The last argument is the help variable z in the
                    % transformation from f(ax+b) to f(z),z==ax+b. We should not
                    % use this transformation if the argument already is unitary
                    hidden = [hidden getvariables(evalMap{i}.arg{j})];
                end
            else
                for j = 1:length(evalMap{i}.arg)
                    % The last argument is the help variable z in the
                    % transformation from f(ax+b) to f(z),z==ax+b. We should not
                    % use this transformation if the argument already is unitary
                    hidden = [hidden getvariables(evalMap{i}.arg{j})];
                end
            end
        end
        used_variables = union(used_variables,hidden);

        % The problem is that linear terms might be missing in problems with only
        % nonlinear expressions
        [monomtable,variabletype] = yalmip('monomtable');
        if (options.relax==1)|(options.relax==3)
            monomtable = [];
            nonlinearvariables = [];
            linearvariables = used_variables;
        else
            nonlinearvariables = find(variabletype);
            linearvariables = used_variables(find(variabletype(used_variables)==0));
        end
        needednonlinear = nonlinearvariables(ismembcYALMIP(nonlinearvariables,used_variables));
        linearinnonlinear = find(sum(abs(monomtable(needednonlinear,:)),1));
        missinglinear = setdiff(linearinnonlinear(:),linearvariables);
        used_variables = uniquestripped([used_variables(:);missinglinear(:)]);


        usedEvalVariables = find(ismember(used_variables,evalVariables));
        evalMap =  yalmip('extstruct',used_variables(usedEvalVariables));
        if ~isa(evalMap,'cell')
            evalMap = {evalMap};
        end
        evalVariables = usedEvalVariables;

        for i = 1:length(evalMap)
            for aux = 1:length(evalMap{i}.arg)-1
                if isa(evalMap{i}.arg{aux},'sdpvar')
                    X = evalMap{i}.arg{aux};
                    break
                end
            end
            n = length(X);
            if isequal(getbase(X),[spalloc(n,1,0) speye(n)])
                index = ismember(used_variables,getvariables(X));
                evalMap{i}.variableIndex = find(index);
            else
                index = ismember(used_variables,getvariables(evalMap{i}.arg{end}));
                evalMap{i}.variableIndex = find(index);
            end
        end
        goon = ~isequal(used_variables,old_used_variables);
        old_used_variables = used_variables;
    end
else
    evalMap = [];
end


function evalVariables = determineEvaluationBased(operators)
evalVariables = [];
for i = 1:length(operators)
    if strcmpi(operators{i}.properties.model,'callback')
        evalVariables = [evalVariables operators{i}.properties.models];
    end
end

function [Fnew,changed] = convertsocp2NONLINEAR(F);
changed = 0;
socps = find(is(F,'socp'));
Fsocp = F(socps);
Fnew = F;
if length(socps) > 0
    changed = 1;
    Fnew(socps) = [];
    for i = 1:length(Fsocp)
        z = sdpvar(Fsocp(i));
        Fnew = [Fnew, z(1)>=0, z(1)^2 >= z(2:end)'*z(2:end)];
    end
end




