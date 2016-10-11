function output = callecos(yalmipmodel)

% Retrieve needed data
options = yalmipmodel.options;

options.ecos.verbose = options.verbose~=0;
if ~isempty(yalmipmodel.binary_variables)
    options.ecos.bool_vars_idx = yalmipmodel.binary_variables;
end
if ~isempty(yalmipmodel.integer_variables)
    options.ecos.int_vars_idx = yalmipmodel.integer_variables;
end
model.opts = options.ecos;


% first setup an scs model
K = yalmipmodel.K;
F_struc = yalmipmodel.F_struc;
ub      = yalmipmodel.ub;
lb      = yalmipmodel.lb;
c = yalmipmodel.c;
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

data.A = -F_struc(:,2:end);
data.b = full(F_struc(:,1));
data.c = full(c);
cones = [];
cones.f = K.f;
cones.l = K.l;
cones.q = K.q;
cones.s = K.s;
cones.ep = 0;
[data,cones,output] = addExponentialCone(data,cones,yalmipmodel);
[data,cones,lift] = compressModel(data,cones,yalmipmodel);
if output.problem == -4
    return
end

% Map to ECOS syntax
model.c = full(data.c);
model.A = data.A(1:cones.f,:);if isempty(model.A);model.A = [];end
model.b = data.b(1:cones.f,:);if isempty(model.b);model.b = [];end
model.G = data.A(1+cones.f:end,:);if isempty(model.G);model.G = [];end
model.h = data.b(1+cones.f:end,:);if isempty(model.h);model.h = [];end
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
    Primal(lift.linearIndex) = x;
    if ~isempty(lift.liftedIndex)
        Primal(lift.liftedIndex) = lift.d + lift.T*x;
    end
    Primal = Primal(1:length(yalmipmodel.c));
    Primal = Primal(:);
else
    Primal = nan(length(yalmipmodel.c),1);
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


% Compress model and remove YALMIPs intermediate auxilliary variables
function [data,cones,lift] = compressModel(data,cones,model)

lift.liftedIndex = [];
lift.linearIndex = 1:length(size(data.A,2));

if isempty(model.aux_variables) || isempty(cones.f) || (cones.f == 0)
    return
end
    
liftedIndex = setdiff(model.aux_variables,model.evalVariables);
if isempty(liftedIndex)
    return
end

linearIndex = setdiff(1:size(data.A,2),liftedIndex);

lift.liftedIndex = liftedIndex;
lift.linearIndex = linearIndex;

lift.d = [];
lift.T = [];
lift.S = [];
if cones.f > 0
    definingLift = find(any(data.A(1:cones.f,liftedIndex),2));
    lift.d = [lift.d;-data.b(definingLift)];
    lift.T = [lift.T;data.A(definingLift,linearIndex)];
    lift.S = [lift.S;-data.A(definingLift,liftedIndex)];
    keep = any(lift.S,1);
    if ~all(keep)
        % FIXME: Sort out this case
        return
    end   
    [i,j,k] = find(lift.S);
    lift.T = lift.T(i,:);
    lift.d = lift.d(i);   
    data.A(definingLift,:) = [];
    data.b(definingLift,:) = [];
    cones.f = cones.f - length(definingLift);
     
    Ax = data.A(:,linearIndex);
    Ay = data.A(:,liftedIndex);
    Ax = Ax + Ay*lift.T;
	b = data.b - Ay*lift.d;
	data.A = Ax;
	data.b = b;
    data.c = data.c(linearIndex) + ((data.c(liftedIndex))'*lift.T)';
end
