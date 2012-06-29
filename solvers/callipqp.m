function output = callipqp(varargin)
% Author Johan Löfberg 
% $Id: callipqp.m,v 1.2 2004-06-30 07:37:03 johanl Exp $


% Retrieve needed data
interfacedata = varargin{1};
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = 2*interfacedata.Q;
ub      = interfacedata.ub;
lb      = interfacedata.lb;


showprogress('Calling IPQP',options.showprogress);

if K.f>0
    A = -F_struc(1:K.f,2:end);
    b = F_struc(1:K.f,1);
else
    A = [];
    b = [];
end

if K.l>0
C = -F_struc(1+K.f:K.f+K.l,2:end);
d = F_struc(1+K.f:K.f+K.l,1);
else
    C = [];
    d = [];
end

if ~isempty(ub)
    C = [C;eye(length(c))];
    d = [d;ub];
end
if ~isempty(lb)
    C = [C;-eye(length(c))];
    d = [d;-lb];
end

solvertime = clock; 
[x,problem] = ipqp(2*Q,c,C,d,A,b);
solvertime = etime(clock,solvertime);

% Internal format for duals
D_struc = [];

switch problem
case 0
    problem = 0;
    infostr = yalmiperror(problem,'IPQP');       
case 1
    problem = 1;
    infostr = yalmiperror(problem,'IPQP');   
otherwise
    
    problem = 2;
    infostr = yalmiperror(problem,'IPQP');       
end    

% Save all data sent to solver?
if options.savesolverinput
%     solverinput.A = A;
%     solverinput.b = b;
%     solverinput.C = Aq;
%     solverinput.q = beq;
%     solverinput.c = c;
%     solverinput.H = Q;
%     solverinput.options = options.nag;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
%     solveroutput.x = x;
%     solveroutput.fmin = fmin;
%     solveroutput.flag = flag;
%     solveroutput.output=output;
%     solveroutput.lambda=lambda;  
else
    solveroutput = [];
end



% Standard interface 
output.x           = x;
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;
output.D_struc = D_struc;