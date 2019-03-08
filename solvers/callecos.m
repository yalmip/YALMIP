function output = callecos(yalmipmodel)
originalModel = yalmipmodel;
% Fix ECOS options
options = yalmipmodel.options;
options.ecos.verbose = options.verbose~=0;
if ~isempty(yalmipmodel.binary_variables)
    options.ecos.bool_vars_idx = yalmipmodel.binary_variables;
end
if ~isempty(yalmipmodel.integer_variables)
    options.ecos.int_vars_idx = yalmipmodel.integer_variables;
end
model.opts = options.ecos;

if ~isempty(yalmipmodel.ub)
    [yalmipmodel.F_struc,yalmipmodel.K] = addStructureBounds(yalmipmodel.F_struc,yalmipmodel.K,yalmipmodel.ub,yalmipmodel.lb);
end

% Write nonlinear functions using exponential cone operators, if possible
[yalmipmodel,output] = normalizeExponentialCone(yalmipmodel);
if output.problem
    return
end

% Map to ECOS syntax
cones = [];
cones.f = yalmipmodel.K.f;
cones.l = yalmipmodel.K.l;
cones.q = yalmipmodel.K.q;
cones.s = yalmipmodel.K.s;
cones.ep = yalmipmodel.K.e;
model.c = full(yalmipmodel.c);
model.A = -yalmipmodel.F_struc(1:cones.f,2:end);if isempty(model.A);model.A = [];end
model.b = full(yalmipmodel.F_struc(1:cones.f,1));if isempty(model.b);model.b = [];end
model.G = -yalmipmodel.F_struc(1+cones.f:end,2:end);if isempty(model.G);model.G = [];end
model.h = full(yalmipmodel.F_struc(1+cones.f:end,1));if isempty(model.h);model.h = [];end
if nnz(cones.l) > 0
    model.dims.l = cones.l;
end
if nnz(cones.q) > 0
    model.dims.q = cones.q;
end
if nnz(cones.ep) > 0
    model.dims.e = cones.ep;
    tempG = model.G(1:model.dims.l+sum(cones.q),:);
    temph = model.h(1:model.dims.l+sum(cones.q));
    top = 1+model.dims.l+sum(cones.q);
    for i = 1:model.dims.e
        tempG = [tempG;model.G(top + [0;2;1],:)];
        temph = [temph;model.h(top + [0;2;1],:)];
        top = top  + 3;
    end
    model.G = tempG;
    model.h = temph;
end

if options.savedebug
    save ecosdebug model
end


if options.showprogress;showprogress(['Calling ' yalmipmodel.solver.tag],options.showprogress);end
if isempty(model.A)   
    solvertime = tic;
    [x,y,info,s,z] = ecos(model.c,model.G,model.h,model.dims,model.opts);  
    solvertime = toc(solvertime);
else    
    solvertime = tic;
    [x,y,info,s,z] = ecos(model.c,model.G,model.h,model.dims,model.A,model.b,model.opts);
    solvertime = toc(solvertime);
end

% Internal format, only keep original variablesop
if ~isempty(x)
    Primal = x(1:length(originalModel.c));
else
    Primal = nan(length(originalModel.c),1);
end
if ~isempty(yalmipmodel.evalMap)
    % No support for duals when exponential cones yet
    Dual = [];
else
    Dual   = [y;z];
end

switch info.exitflag
    case 0
        problem = 0;
    case 1
        problem = 1;
    case 2
        problem = 2;
    case -1
        problem = 3;
    case {-2,-3}
        problem = 4;
    case -7
        problem = 9;
    case {10,11}
        problem = 3;
    otherwise
        problem = 9;
end

infostr = yalmiperror(problem,yalmipmodel.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.model = model;
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.y = y;
    solveroutput.info = info;
    solveroutput.s = s;
    solveroutput.z = z;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);