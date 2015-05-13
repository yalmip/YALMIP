function output = callnage04mbf(varargin)

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
sent_bounds = 1;
if isempty(ub)
    sent_bounds = 0;
    ub =  bigbnd*ones(length(c),1);
    lb = -bigbnd*ones(length(c),1);
end   

showprogress('Calling NAG',options.showprogress);

if ~isempty(F_struc)
    A =-F_struc(1:end,2:end);
    beq = F_struc(1:1:K.f,1);
    b = F_struc(K.f+1:end,1);   
    
    ub = full([ub;beq;b]);
    lb = full([lb;beq;-bigbnd*ones(length(b),1)]);
else
    A=zeros(1,length(c));
end

lp   = 0;
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
solvertime = tic;
[x,istate,objlp,clambda,ifail] = e04mbf(full(lb),full(ub),zeros(length(c),1),full(c),full(A),msglev,options.nag.itmax,ifail);%options.nag.bigbnd,options.nag.orthog);
solvertime = toc(solvertime);
problem = 0;

% Internal format for duals
D_struc = -clambda(length(c)+1:end);

switch ifail
case {-1,0}
    problem = 0;
case 1
    problem = 1;
case 2
    problem = 2;
case 3
    problem = 4;
case 4
    problem = 3;
case 5
    problem = 7;
otherwise
    problem = -1;
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
    solveroutput.objlp = objlp;
    solveroutput.clambda=clambda;
    solveroutput.istate=istate;  
    solveroutput.ifail=ifail;  
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);