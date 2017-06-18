function output = callpenlabbmi(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
Q       = model.Q;
K       = model.K;
x0      = model.x0;
monomtable = model.monomtable;
ub      = model.ub;
lb      = model.lb;

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end
bmimodel.penstruct = sedumi2penbmi(F_struc,full(c),2*Q,K,monomtable,options,x0);
penlabmodel=yalmip2bmi(bmimodel);
penlabmodel = bmi_define(penlabmodel);
prob = penlab(penlabmodel);
switch options.verbose
    case 0
        options.penlab.outlev = 0;
    otherwise
        options.penlab.outlev = 1+options.verbose;
end
prob.opts = options.penlab;
showprogress('Calling PENLAB',model.options.showprogress);

solvertime = tic;
solve(prob);
solvertime = toc(solvertime);

xout = prob.x;
x = zeros(length(model.c),1);
if ~isempty(xout)
    x(find(model.variabletype==0))=xout;
end

% Duals currently not supported
D_struc = [];

switch 0
    case {0,1}
        problem = 0;
    case {2}
        problem = 1;
    case {-1}
        problem = 3;
    case {3,4,-2,-3}
        problem = 4;
    case {-11,-12,-13}
        problem = 7;
    case {-10,-100,-101,-102,-199}
        problem = 11;
    otherwise
        problem = -1;
end

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.prob = prob;
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.prob = prob;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'PENLAB',solverinput,solveroutput,solvertime);


