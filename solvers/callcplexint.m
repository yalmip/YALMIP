function output = callcplexint(interfacedata)

% Author Johan Löfberg 
% $Id: callcplexint.m,v 1.21 2009-11-03 11:08:47 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
UB      = interfacedata.ub;
LB      = interfacedata.lb;

showprogress('Calling CPLEXINT',options.showprogress);

SENSE = 1;     


n_original = length(c);
if K.q(1)>0
    % To simplyfy code, we currently normalize everything to z'*z<x0^2
    nNEW = sum(K.q);
    if ~isempty(x0)
        x0 = [x0;full(F_struc(1+K.f+K.l:end,:))*[1;x0]];
    end
        
    F_strucSOCP = [F_struc(1+K.f+K.l:end,:) -speye(nNEW)];
    F_struc = [F_struc(1:K.f+K.l,:) spalloc(K.f+K.l,nNEW,0)];
    UB = [UB;inf(nNEW,1)];
    c = [c;spalloc(nNEW,1,0)];
    Q = blkdiag(Q,spalloc(nNEW,nNEW,0));
    
    iCone = n_original+1;
  %  ri = zeros(1,length(K.q));
  %  Li = spalloc(n_original+nNEW,length(K.q),0);
    for i = 1:length(K.q);
        QC(i).Q = (sparse(iCone:iCone+K.q(i)-1,iCone:iCone+K.q(i)-1,[-1 ones(1,K.q(i)-1)],n_original+nNEW,n_original+nNEW));
        QC(i).r = 0;
        QC(i).L = zeros(n_original+nNEW,1)';
        LB = [LB;0;-inf(K.q(i)-1,1)];
        iCone = iCone + K.q(i);
    end
    F_struc = [F_strucSOCP;F_struc];
    K.f = K.f + nNEW;
else
  QC = [];
end

C = full(c);   
if K.l+K.f+K.q == 0
    A = zeros(1,length(c));A(1)=1;
    B = 1e6;
else
    A =-(F_struc(1:K.f+K.l,2:end)); 
    B = full(F_struc(1:K.f+K.l,1));            
end

INDEQ = [];
if K.f>0
    INDEQ(1:K.f) = 1:K.f;
end

VARTYPE = repmat('C',size(A,2),1);
VARTYPE(integer_variables)='I'; % Integer variables
VARTYPE(binary_variables) ='B'; % Binary variables

if nnz(Q)==0
    H = [];
else
    H = full(2*Q);
end

PARAM = options.cplex.param;
OPTIONS.verbose = options.verbose;
OPTIONS.logfile = options.cplex.logfile;
if ~isempty(x0)
    OPTIONS.x0 = [(1:length(x0))' x0(:)];
end

if options.savedebug
    save cplexintdebug H C A B LB UB QC VARTYPE INDEQ PARAM OPTIONS
end

if ~isempty(PARAM.double)
     i = find(PARAM.double(:,1)==1025);
     if ~isempty(i)
         PARAM.double(i,2) = PARAM.double(i,2)-interfacedata.f;
     end
     i = find(PARAM.double(:,1)==1026);
     if ~isempty(i)
         PARAM.double(i,2) = PARAM.double(i,2)-interfacedata.f;
     end
end


% Call mex-interface
solvertime = clock; 
[x,FMIN,SOLSTAT,DETAILS] = cplexint(H, C, A, B, INDEQ, QC, LB, UB,VARTYPE,PARAM,OPTIONS);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end
problem = 0;

D_struc = -DETAILS.dual;    
if K.q(1)>0
    D_struc=[];
end

if isempty(x)
    x = zeros(n_original,1);
else
    x = x(1:n_original);
end

switch SOLSTAT
    case {1,101,102}
        problem = 0;
    case {3,22,103}
        problem = 1;
    case 108
        problem = 3;
    case {2}
        problem = 2;
    case {4,119}
        problem = 12;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'CPLEXINT');	

% Save all data sent to solver?
if options.savesolverinput
    solverinput.H = H;
    solverinput.A = A;
    solverinput.C = C;
    solverinput.INDEQ = INDEQ;
    solverinput.QC = QC;
    solverinput.B = B;
    solverinput.VARTYPE = VARTYPE;
    solverinput.LB = LB;
    solverinput.UB = UB;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.FMIN = FMIN;
    solveroutput.SOLSTAT = SOLSTAT;
    solveroutput.DETAILS=DETAILS;
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