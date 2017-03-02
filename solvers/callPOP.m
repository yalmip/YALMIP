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

ops = interfacedata.options.pop;

if interfacedata.options.verbose
    if isequal(lower(ops.Progress),'on') || isequal(lower(ops.Progress),'off')
        % User hasn't specified fancy progress, so we simply turn on text
        ops.Progress = 'text';
    end
else
    ops.Progress = 'Off';
end

if interfacedata.options.savedebug
    save popdebug Matrices ops
end

if ~isempty(Matrices.E)
    solution = mpMIQP(Matrices,ops);
else
    solution = mpQP(Matrices,ops);
end

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