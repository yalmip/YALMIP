function output = callcdd(interfacedata)

% Author Johan Löfberg 
% $Id: callcdd.m,v 1.3 2005-05-07 13:53:20 joloef Exp $

% Standard input interface
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
Q       = interfacedata.Q;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

if ~isempty(ub)
    [F_struc,K] = addbounds(F_struc,K,ub,lb);
end

IN = struct('obj',full(c(:))','A',full(-F_struc(:,2:end)),'B',full(F_struc(:,1)));
if K.f>0
    IN.lin = 1:K.f;
end

% % CDD does not check this..
% % Fix the duals...
% valid_constraints = find(any(IN.A'));
% IN.A = IN.A(valid_constraints,:);
% IN.B = IN.B(valid_constraints,:);

showprogress('Calling CDD',options.showprogress);
solvertime = clock; 

if options.savedebug
    save cdddebug IN
end

switch options.cdd.method
    case 'criss-cross'
        OUT = cddmex('solve_lp',IN); 
    case 'dual-simplex'
        OUT = cddmex('solve_lp_DS',IN);
    otherwise
end
solvertime = etime(clock,solvertime);
problem = 0;

% Internal format for duals
D_struc = [-OUT.lambda];

switch OUT.how
    case 1
        problem = 0;
    case {2,7}
        problem = 1;
    case {3,6} 
        problem = 2;
    case {0,4,5}
        problem = 11;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'CDD');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.IN = IN;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.OUT = OUT;
else
    solveroutput = [];
end

% Standard interface 
output.Primal      = OUT.xopt;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;