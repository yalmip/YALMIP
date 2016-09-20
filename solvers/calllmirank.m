function output = calllmirank(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

pars = options.lmirank;
pars.fid = double(options.verbose);

if ~options.usex0
    interfacedata.options.verbose =  max(0,interfacedata.options.verbose - 1);
    initialsolver = eval(['@' interfacedata.solver.initialsolver.call]);
    start = K.l;
    b = zeros(size(F_struc,2)-1,1);
    for j=1:length(K.s)
        if K.s(j)~=K.rank(j)
            ind = find(speye(K.s(j)));           
            for i=1:size(F_struc,2)-1
                b(i,1) = b(i,1) + sum(F_struc(start+ind,1+i));
            end
        end
        start=start+K.s(j)^2;
    end
    interfacedata.c = b;
    interfacedata.solver = interfacedata.solver.initialsolver;
    output = feval(initialsolver,interfacedata);
    if output.problem ~= 1
        options.usex0 = 1;       
        x0 = output.Primal;
    else
        return;
    end
end

if options.savedebug
    At = -F_struc(:,2:end);
    C  = F_struc(:,1);
    save lmirankdebug C At K pars x0
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
if options.usex0
    [y,info] = lmirank(-F_struc(:,2:end),F_struc(:,1),K,pars,x0);
else
    [y,info] = lmirank(-F_struc(:,2:end),F_struc(:,1),K,pars);
end
solvertime = toc(solvertime);
x = y;

switch info.solved
    case 1
        problem = 0;
    otherwise
        problem = 11;
end

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.At = -F_struc(:,2:end);
    solverinput.c = F_struc(:,1);
    solverinput.K = K;  
    solverinput.pars = pars;
    solverinput.x0 = x0;
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.y = y;
    solveroutput.info = info;    
else
    solveroutput = [];
end

% Standard interface 
infostr = yalmiperror(problem,interfacedata.solver.tag);	
output = createOutputStructure(x,[],[],problem,infostr,solverinput,solveroutput,solvertime);