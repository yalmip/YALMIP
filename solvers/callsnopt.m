function output = callsnopt(model)

% Build trees etc. This common format is used for fmincon, snopt and ipopt
% interfaces
model = yalmip2nonlinearsolver(model);

if ~model.derivative_available
    disp('Derivate-free call to snopt not yet implemented')
    error('Derivate-free call to snopt not yet implemented')
end

% Figure out which variables are artificially introduced to normalize
% arguments in callback operators to simplify chain rules etc. We can do
% our function evaluations and gradient computations in our lifted world,
% but only expose the model in the original variables to the nonlinear
% solver. 
% model = compressLifted(model);

if model.options.savedebug
    save snoptdebug model
end

showprogress('Calling SNOPT',model.options.showprogress);

x0 = model.x0;
solveopt = 1;
xlow = model.lb;
xupp = model.ub;
xmul   = zeros(length(xupp),1);
xstate = zeros(length(xupp),1);

Fupp = [ inf;
         repmat(0,length(model.bnonlinineq),1);
         repmat(0,nnz(model.K.q),1);
         repmat(0,length(model.bnonlineq),1);
         repmat(0,length(model.b),1);
         repmat(0,length(model.beq),1)];
     
Flow = [ -inf;
         repmat(-inf,length(model.bnonlinineq),1);
         repmat(-inf,nnz(model.K.q),1);
         repmat(0,length(model.bnonlineq),1);
         repmat(-inf,length(model.b),1);
         repmat(0,length(model.beq),1)];

Fmul   = zeros(length(Fupp),1);
Fstate = zeros(length(Fupp),1);
ObjAdd = 0;
ObjRow = 1;
A = [];
iAfun = [];
jAvar = [];
% Sparsity pattern of jacobian
if ~isempty(Fupp)
    G = jacobiansparsityfromnonlinear(model);   
else
    G = [];
end
% Add a row for objective. No sparsity declared
G = [ones(1,size(G,2));G];
[iGfun,jGvar] = find(G);
model.sparsityElements = find(G);

usrf = 'snopt_callback';
snopt_callback([],model);

global latest_xevaled
global latest_x_xevaled
latest_xevaled = [];
latest_x_xevaled = [];

solvertime = tic;
if strcmpi(model.solver.version,'cmex')
    % Some old interface? Keep for safety
    if model.options.verbose == 0
        snscreen('off')
    else
        snscreen('on');
    end
    snseti('Minimize',1)
    [xout,F,xmul,Fmul,inform, xstate, Fstate, ns, ninf, sinf, mincw, miniw, minrw] = snoptcmex( solveopt, x0, xlow, xupp, xmul, xstate, Flow, Fupp, Fmul, Fstate,ObjAdd, ObjRow, A, iAfun(:), jAvar(:),iGfun(:), jGvar(:), usrf );
else
    Astruct.A = A;
    Astruct.row = iAfun;
    Astruct.col = jAvar;
    if model.options.verbose == 0
       model.options.snopt.screen = 'off';
    end
    [xout,F,inform,xmul,Fmul] = snopt(x0, xlow, xupp, xmul, xstate,Flow, Fupp, Fmul, Fstate, usrf,ObjAdd, ObjRow,A,G,model.options.snopt);
end
  
solvertime = toc(solvertime);

lambda = Fmul(2:end);   

if ~isempty(xout) && ~isempty(model.lift);
    x = zeros(length(model.linearindicies),1);
    x(model.lift.linearIndex) = xout;
    x(model.lift.liftedIndex) = model.lift.T*xout + model.lift.d;
    x = RecoverNonlinearSolverSolution(model,x);
else
    x = RecoverNonlinearSolverSolution(model,xout);
end

problem = 0;

% Internal format for duals
D_struc = [];

% Check, currently not exhaustive...
switch inform
    case {1}
        problem = 0;
    case {11,12,13,14,40,43,91} % 1 is sent when I test
        problem = 1;
    case 21
        problem = 2;
    case {31,32}
        problem = 3;
    case {2,3,33,41,42,43,44}
        problem = 4;
    case {71,72,73,74}
        problem = 16;
    otherwise 
        problem = 11;
end

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.model = model; 
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = x;
    solveroutput.F = F;
    solveroutput.xmul = xmul;
    solveroutput.inform=inform;
    solveroutput.Fstate=Fstate;
    solveroutput.ns = ns;
    solveroutput.ninf = ninf;
    solveroutput.sinf = sinf;
    solveroutput.mincw = mincw;
    solveroutput.miniw = miniw;
    solveroutput.miniw = miniw;
    
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'SNOPT',solverinput,solveroutput,solvertime);