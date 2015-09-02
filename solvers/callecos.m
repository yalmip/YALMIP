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
    [F_struc,K] = addbounds(F_struc,K,ub,lb);
end

data.A = -F_struc(:,2:end);
data.b = full(F_struc(:,1));
data.c = full(c);
cones = [];
cones.f = K.f;
cones.l = K.l;
cones.q = K.q;
cones.s = K.s;
[data,cones] = addExponentialCone(data,cones,yalmipmodel);

% Map to ECOS syntax
model.A = data.A(1:cones.f,:);if isempty(model.A);model.A = [];end
model.b = data.b(1:cones.f,:);if isempty(model.b);model.b = [];end
model.G = data.A(1+cones.f:end,:);if isempty(model.G);model.G = [];end
model.h = data.b(1+cones.f:end,:);if isempty(model.h);model.h = [];end
model.dims.l = cones.l;if nnz(model.dims.l)==0;model.dims.l = [];end
model.dims.q = cones.q;if nnz(model.dims.q)==0;model.dims.q = [];end
model.dims.e = cones.ep;if nnz(model.dims.e)==0;model.dims.e = [];end

if options.savedebug
    save ecosdebug model
end


if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end
if isempty(model.A)   
    solvertime = tic;
    [x,y,info,s,z] = ecos(model.c,model.G,model.h,model.dims,model.opts);  
    solvertime = toc(solvertime);
else    
    solvertime = tic;
    [x,y,info,s,z] = ecos(model.c,model.G,model.h,model.dims,model.A,model.b,model.opts);
    solvertime = toc(solvertime);
end

% Internal format
Primal = x;
Dual   = [y;z];

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