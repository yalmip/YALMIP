function output = calllsqlin(interfacedata)

% Author Johan Löfberg
% $Id: calllsqlin.m,v 1.17 2007-08-02 11:39:36 joloef Exp $

K = interfacedata.K;
c = interfacedata.c;
CA = interfacedata.F_struc;
options = interfacedata.options;

% To begin with, with try to figure out if this is a simple non-negative
% least squares in disguise
if length(K.q)~=1
    output = error_output;
    return
end

% total number of variables, #x+1 since YALMIP has introduced an epi-graph
% variable to model norm(Cx-d)
n = size(CA,2)-1;

if nnz(c)~=1
    output = error_output;
    return
end
if c(end)~=1
    output = error_output;
    return
end
if c(end)~=1
    output = error_output;
    return
end

Cones = CA(K.f+K.l+1:end,:);
LPs = CA(1:K.l+K.f,:);

if any(LPs(:,end))
    % Epigraph involved in LPs
    output = error_output;
    return
end

if isempty(LPs)
    Aeq = [];
    beq = [];
    A = [];
    b = [];
else
    Aeq = -LPs(1:1:K.f,2:end-1);
    beq = LPs(1:1:K.f,1);        
    A =-LPs(K.f+1:end,2:end-1);
    b = LPs(K.f+1:end,1);   
end

% Is it really cone(C*x-d,t)
if ~all(Cones(1,:) == [zeros(1,n) 1])
    output = error_output;
    return
end
if ~all(Cones(:,end) == [1;zeros(K.q-1,1)])
    output = error_output;
    return
end

model.d = -Cones(2:end,1);
model.C = Cones(2:end,2:end-1);
model.Aineq = A;
model.Aeq = Aeq;
model.bineq = b;
model.beq = beq;
model.x0 = [];
model.options = interfacedata.options.lsqlin;
model.solver = 'lsqlin';

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock;
[X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = lsqlin(model);
solveroutput.X = X;
solveroutput.RESNORM = RESNORM;
solveroutput.RESIDUAL = RESIDUAL;
solveroutput.EXITFLAG = EXITFLAG;
solveroutput.OUTPUT = OUTPUT;
solveroutput.LAMBDA = LAMBDA;
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

solution.x = [X;RESNORM];
solution.D_struc = [];

switch EXITFLAG
    case 1
        solution.problem = 0;
    case 2
        solution.problem = 1;
    case 0
        solution.problem = 3;
    case {3,-4,-7}
        solution.problem = 4;
    otherwise
        solution.problem = 9;
end

% Save all data sent to solver?
if ~options.savesolverinput
    model = [];
end

% Save all data from the solver?
if ~options.savesolveroutput
    solveroutput = [];
end

% Standard interface
output.Primal      = solution.x(:);
output.Dual        = solution.D_struc;
output.Slack       = [];
output.problem     = solution.problem;
output.infostr     = yalmiperror(solution.problem,interfacedata.solver.tag);
output.solvertime  = solvertime;
if ~options.savesolverinput
    output.solverinput = [];
else
    output.solverinput = model;
end
if ~options.savesolveroutput
    output.solveroutput = [];
else
    output.solveroutput = solveroutput;
end


function output = error_output
output.Primal      = [];
output.Dual        = [];
output.Slack       = [];
output.problem     = 9;
output.infostr     = yalmiperror(-4,'LSQLIN');
output.solvertime  = 0;


