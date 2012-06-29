function output = calllsqnonneg(interfacedata)

% Author Johan Löfberg
% $Id: calllsqnonneg.m,v 1.17 2007-08-02 11:39:36 joloef Exp $

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

if K.l~=n-1
    output = error_output;
    return
end

Cones = CA(K.l+1:end,:);
LPs = CA(1:K.l,:);

if any(LPs(:,end))
    % Epigraph involved in LPs
    output = error_output;
    return
end

% Are we constraining all x>0
if any(LPs(:,1))
    output = error_output;
    return
end
[i,j,k] = find(LPs(:,2:end));
if ~all(k==1)
    output = error_output;
    return
end
if ~all(LPs*ones(n+1,1)==1)
    output = error_output;
    return
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
model.options = interfacedata.options.lsqnonneg;
if isempty(interfacedata.x0)
    model.x0 = [];
else
    model.x0 = interfacedata.x0(1:end-1);
end

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock;
try
    [X,RESNORM,RESIDUAL,EXITFLAG] = lsqnonneg(model.C,model.d,model.options);
catch
    [X,RESNORM,RESIDUAL,EXITFLAG] = lsqnonneg(model.C,model.d,model.x0,model.options);
end
solveroutput.X = X;
solveroutput.RESNORM = RESNORM;
solveroutput.RESIDUAL = RESIDUAL;
solveroutput.EXITFLAG = EXITFLAG;
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

solution.x = [X;RESNORM];
solution.D_struc = [];

switch EXITFLAG
    case 1
        solution.problem = 0;
    case 0
        solution.problem = 3;
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


function output = error_output(e)
if nargin == 0   
    e = -4;
end
output.Primal      = [];
output.Dual        = [];
output.Slack       = [];
output.problem     = e;
output.infostr     = yalmiperror(e,'LSQNONNEG');
output.solvertime  = 0;


