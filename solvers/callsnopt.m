function output = callsnopt(model)

% Build trees etc. This common format is used for fmincon, snopt and ipopt
% interfaces
model = yalmip2nonlinearsolver(model);

if ~model.derivative_available
    disp('Derivate-free call to snopt not yet implemented')
    error('Derivate-free call to snopt not yet implemented')
end
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
if model.options.verbose == 0
    snscreen('off')
else
    snscreen('on');
end
snseti('Minimize',1)  

tic
[xout,F,xmul,Fmul,inform, xstate, Fstate, ns, ninf, sinf, mincw, miniw, minrw] = snoptcmex( solveopt, x0, xlow, xupp, xmul, xstate, Flow, Fupp, Fmul, Fstate,ObjAdd, ObjRow, A, iAfun(:), jAvar(:),iGfun(:), jGvar(:), usrf );
solvertime = toc;
   
lambda = Fmul(2:end);   

x = RecoverNonlinearSolverSolution(model,xout);

problem = 0;

% Internal format for duals
D_struc = [];

% Check, currently not exhaustive...
switch inform
    case {1}
        problem = 0;
    case {1,11,12,13,14,40,43,91} % 1 is sent when I test
        problem = 1;
    case {2,33}
        problem = 4;
    otherwise        
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