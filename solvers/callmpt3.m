function output = callmpt3(interfacedata)

% Speeds up solving LPs in mpmilp
global MPTOPTIONS
if ~isstruct(MPTOPTIONS)
    mpt_error
end

% Convert
Matrices = yalmip2mpt(interfacedata);

% Get some MPT options
options = interfacedata.options;
options.mpt.lpsolver = MPTOPTIONS.lpsolver;
options.mpt.milpsolver = MPTOPTIONS.milpsolver;
options.mpt.verbose = options.verbose;

if options.savedebug
    save mptdebug Matrices
end

if options.mp.unbounded
    Matrices = removeExplorationConstraints(Matrices);
end

[dummy,un] = unique([Matrices.G Matrices.E Matrices.W],'rows');
Matrices.G = Matrices.G(un,:);
Matrices.E = Matrices.E(un,:);
Matrices.W = Matrices.W(un,:);

if isempty(Matrices.binary_var_index)

    showprogress('Calling MPT',options.showprogress);
    solvertime = tic;
    if options.mp.presolve
        [Matrices.lb,Matrices.ub] = mpt_detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);
    end        
    
    if any(Matrices.lb(end-Matrices.nx+1:end) == Matrices.ub(end-Matrices.nx+1:end))
        model = [];
    else        
        model = mpt_solvenode(Matrices,Matrices.lb,Matrices.ub,Matrices,[],options);
    end
    solvertime = toc(solvertime);

else  
    % Pre-solve required on binary problems
    options.mp.presolve = 1;

    solvertime = tic;        
    switch options.mp.algorithm
        case 1
            showprogress('Calling MPT via enumeration',options.showprogress);
            model = mpt_enumeration_mpmilp(Matrices,options);
        case 2
            % Still experimental and just for fun. Not working!
            showprogress('Calling MPT via parametric B&B',options.showprogress);
            model = mpt_parbb(Matrices,options);         
            
       case 3
            showprogress('Calling MPT via delayed enumeration',options.showprogress);
            [Matrices.SOS,Matrices.SOSVariables] =  mpt_detect_sos(Matrices);
            [Matrices.lb,Matrices.ub] = mpt_detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);                                
            model = mpt_de_mpmilp(Matrices,options,[]);                                                
        otherwise
            
    end
    solvertime = toc(solvertime);
end

if isempty(model)
    model = {model};
end

if options.verbose
    if ~isempty(model{1})
        if length(model) == 1
            disp(['-> Generated 1 partition.'])            
        else
            disp(['-> Generated ' num2str(length(model)) ' partitions.'])
        end
    end
end

problem = 0;

% Save all data sent to solver?
if options.savesolverinput
    solverinput.Matrices = Matrices;
    solverinput.options  = [];
else
    solverinput = [];
end

% Save all data from the solver?
% This always done
if options.savesolveroutput
    solveroutput.model = model;
    solveroutput.U = interfacedata.used_variables(Matrices.free_var);%(Matrices.free_var <= length( interfacedata.used_variables)));
    solveroutput.x = interfacedata.used_variables(Matrices.param_var);
else
    solveroutput = [];
end

% Standard interface
Primal      = nan*ones(length(interfacedata.c),1);
Dual        = [];
output = createOutputStructure(Primal,Dual,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);

function Matrices = removeExplorationConstraints(Matrices);
candidates = find((~any(Matrices.G,2)) & (sum(Matrices.E | Matrices.E,2) == 1));
if ~isempty(candidates)
    Matrices.bndA = -Matrices.E(candidates,:);
    Matrices.bndb = Matrices.W(candidates,:);
    Matrices.G(candidates,:) = [];
    Matrices.E(candidates,:) = [];
    Matrices.W(candidates,:) = [];
end
