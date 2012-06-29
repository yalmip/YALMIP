function output = callclp(interfacedata)

% Author Johan Löfberg 
% $Id: callclp.m,v 1.16 2010-03-14 12:57:16 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
Q       = interfacedata.Q; 

n = length(c);

if options.showprogress;showprogress(['Calling CLP'],options.showprogress);end

if K.f>0
    Aeq = -F_struc(1:K.f,2:end);
    beq = F_struc(1:K.f,1);
else
    Aeq = [];
    beq = [];
end

if K.l > 0
    A = -F_struc(1+K.f:end,2:end);
    b = F_struc(1+K.f:end,1);
else
    A = [];
    b = [];
end

% Fix for bug in mexclp or clp
fixed = 0;
if size(A,1)==0 & size(Aeq,2)==0
    A = [ones(1,length(c))];
    if ~isempty(ub)
        A(isinf(ub)) = 0;        
        dummy = ub;
        dummy(isinf(ub)) = 0;
        b = sum(dummy);
    else
        b = 1e8;
    end
    fixed = 1;
    options.saveduals = 0;
end
% lb(lb==-inf)=-1e12
% ub(ub==inf)=1e12

ops = options.clp;
ops.verbose = options.verbose;

if options.savedebug
    save clpdebug
end

% Call mex-interfacec
solvertime = clock; 
[x,lambda,problem] = clp(2*Q,c,A,b,Aeq,beq,lb,ub,ops);%,interfacedata.integer_variables);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

if options.saveduals
    D_struc = -lambda;    
else
    D_struc = [];
end

infostr = yalmiperror(problem,interfacedata.solver.tag);

% Save all data sent to solver?
if options.savesolverinput
	solverinput = [];	
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
	solveroutput = [];	
else
	solveroutput = [];
end

% Standard interface 
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;