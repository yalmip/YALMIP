function output = callPOP(interfacedata)

% temporarily convert from YALMIP to MPT format
Matrices = yalmip2mpt(interfacedata);

% Some preprocessing to remove bad stuff
Matrices = removeExplorationConstraints(Matrices);
[dummy,un] = unique([Matrices.G Matrices.E Matrices.W],'rows');
Matrices.G = Matrices.G(un,:);
Matrices.E = Matrices.E(un,:);
Matrices.W = Matrices.W(un,:);

% Convert from MPT format to POP format
Matrices = mpt2pop(Matrices);

if interfacedata.options.savedebug
    save popdebug Matrices
end

solution = mpQP(Matrices);
solvertime=0;
problem = 0;
infostr = yalmiperror(problem,'POP');

% Save all data sent to solver?
if interfacedata.options.savesolverinput
    solverinput.Matrices = Matrices;
    solverinput.options  = [];
else
    solverinput = [];
end

% Save all data from the solver?
% This always done
if interfacedata.options.savesolveroutput
    solveroutput.model = {solution};
else
    solveroutput = [];
end

% Standard interface
Primal      = nan*ones(length(interfacedata.c),1);
Dual        = [];
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);