function output = callqsopt(interfacedata)

% Author Johan Löfberg 
% $Id: callqsopt.m,v 1.3 2005-05-07 13:53:20 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

n = length(c);
% Bounded variables converted to constraints
if ~isempty(ub)
    % TODO : MAKE SURE SIZES ARE CORRECT  
    lb(lb==-inf) = -1e30;
    ub(ub==inf) = 1e30;
else
    lb=repmat(-1e30,n,1);    % just for sure
    ub=repmat(1e30,n,1);         
end

if options.showprogress;showprogress('Calling QSOPT',options.showprogress);end

options.qsopt.verbose = options.verbose;

% Call mex-interface
if options.savedebug
    save qsoptdebug
end
solvertime = clock; 

[x,lambda,STATUS] = qsopt(c,-F_struc(1+K.f:end,2:end),F_struc(1+K.f:end,1),-F_struc(1:K.f,2:end),F_struc(1:K.f,1),lb,ub,options.qsopt);
solvertime = etime(clock,solvertime);
problem = 0;

if options.saveduals
    D_struc = -lambda;    
else
    D_struc = [];
end

% Check, currently not exhaustive...
switch STATUS
    case 1
        problem = 0;
    case 2
        problem = 1;
    case 3
        problem = 2;
    case 4        
        problem = 3;
    case {5,6}
        problem = 4;  
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'QSOPT');	

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
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;