function output = callnage04naf(varargin)
% Author Johan Löfberg 
% $Id: callnage04naf.m,v 1.12 2006-04-10 09:34:47 joloef Exp $

% Hack for NAG
persistent Q

if nargin>1
    output = Q*varargin{6};
    return
end

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

bigbnd = options.nag.bigbnd;
if isempty(ub)
    ub =  bigbnd*ones(length(c),1);
    lb = -bigbnd*ones(length(c),1);
end

% NAG fails on semidefinite
if nnz(Q)>0
    ii = find(diag(Q)==0);
    if ~isempty(ii)
        d = zeros(length(Q),1);
        d(ii) = sqrt(eps);
        Q = Q + diag(d);
    end
end

showprogress('Calling NAG',options.showprogress);

if ~isempty(F_struc)
    A =-F_struc(1:end,2:end);
    beq = F_struc(1:1:K.f,1);
    b = F_struc(K.f+1:end,1);   
    
    ub = full([ub;beq;b]);
    lb = full([lb;beq;-bigbnd*ones(length(b),1)]);
elseif ~isempty(lb)% NAG stinks on problems with no A<b constraints
    A  = eye(length(c));
    ub = full([ub;ub]);
    lb = full([lb;lb]);
else
    A  = zeros(1,length(c));
    ub = full([ub;1]);
    lb = full([lb;0]);    
end

% Bug prevents lp=full(nnz(...
if nnz(Q)==0
    lp = 1;
else
    lp = 0;
end

cold = 1;
istate = zeros(length(ub),1);
featol = options.nag.featol*ones(length(ub),1);
switch options.verbose
case 0
    msglev = -1;
    ifail = 1;
case 1
    msglev = 2;
    ifail = -1;
otherwise
    msglev = 5*options.verbose;
    ifail = -1;
end
solvertime = clock; 

[x,iter,obj,clambda,istate,ifail] = e04naf(full(lb),full(ub),'callnage04naf',zeros(length(c),1),full(c),full(A),0,lp, cold,istate,featol,msglev,options.nag.itmax,options.nag.bigbnd,options.nag.orthog,ifail);

solvertime = etime(clock,solvertime);
problem = 0;

% Internal format for duals
D_struc = -clambda(length(c)+1:end);

switch ifail
case {-1,0,1,3}
    problem = 0;
case 2
    problem = 2;
case {4,7}
    problem = 4;
case {5,8}
    problem = 3;
case 6
    problem = 1;
otherwise
    problem = 9;
end    
infostr = yalmiperror(problem,'NAG');       

% Save all data sent to solver?
if options.savesolverinput
    solverinput.bl = full(lb);
    solverinput.bu = full(ub);
    solverinput.X = zeros(length(c),1);
    solverinput.cvec = full(c);
    solverinput.A = full(A);
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.iter = iter;
    solveroutput.obj = obj;
    solveroutput.clambda=clambda;
    solveroutput.istate=istate;  
    solveroutput.ifail=ifail;  
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


% NAG-hack
clear Q