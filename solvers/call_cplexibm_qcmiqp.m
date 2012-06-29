function output = call_cplexibm_miqcp(interfacedata)

% Author Johan Löfberg
% $Id: call_cplexibm_miqcp.m,v 1.21 2009-11-03 11:08:47 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;

UB      = interfacedata.ub;
LB      = interfacedata.lb;

showprogress('Calling CPLEX',options.showprogress);

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [LB,UB,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,LB,UB);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
    
    redundant = find(LB<=0 & UB>=0);
    semicont_variables = setdiff(semicont_variables,redundant);
end

n_original = length(c);
if K.q(1)>0
    % To simplify code, we currently normalize everything to z'*z<x0^2
    if 1
        nNEW = sum(K.q);
        if ~isempty(x0)
            x0 = [x0;full(F_struc(1+K.f+K.l:end,:))*[1;x0]];
        end
        
        Ftemp = F_struc(1+K.f+K.l:end,:);
        F_strucSOCP = [Ftemp -speye(nNEW)];
        F_struc = [F_struc(1:K.f+K.l,:) spalloc(K.f+K.l,nNEW,0)];
        UB = [UB;inf(nNEW,1)];
        c = [c;spalloc(nNEW,1,0)];
        Q = blkdiag(Q,spalloc(nNEW,nNEW,0));
        
        iCone = n_original+1;
        ri = zeros(1,length(K.q));
        Li = spalloc(n_original+nNEW,length(K.q),0);
        for i = 1:length(K.q);
            Qi{i} = sparse(iCone:iCone+K.q(i)-1,iCone:iCone+K.q(i)-1,[-1 ones(1,K.q(i)-1)],n_original+nNEW,n_original+nNEW);
            LB = [LB;0;-inf(K.q(i)-1,1)];
            iCone = iCone + K.q(i);
        end
        F_struc = [F_strucSOCP;F_struc];
        K.f = K.f + nNEW;
    else
        rhsROWS = K.f + K.l + [1 1+cumsum([K.q(1:end-1)])];
        nNEW = length(K.q);
        rhsDATA = F_struc(rhsROWS,:);
        F_struc(rhsROWS,:) = 0;
        F_struc = [rhsDATA -speye(nNEW);
        F_struc [spalloc(K.f+K.l,nNEW,0);spalloc(sum(K.q),nNEW,0)]];
        K.f = K.f + nNEW;
        UB = [UB;inf(nNEW,1)];
        LB = [LB;zeros(nNEW,1)];
        c = [c;spalloc(nNEW,1,0)];
        Q = blkdiag(Q,spalloc(nNEW,nNEW,0));
        
        iCone = n_original+1;
        ri = zeros(1,length(K.q));
        top = K.l+K.f + 1;
        Li = [];
        for i = 1:length(K.q);
            data = F_struc(top:top + K.q(i)-1,:);
            b = data(:,1);
            A = data(:,2:end);
            ri(i) = -b'*b;
            Li = [Li 2*A'*b];
            Qtemp = A'*A;Qtemp(n_original + i,n_original+i)=-1;
            Qi{i} = Qtemp;
            iCone = iCone + K.q(i);
            top = top + K.q(i);
        end
        
    end
else
    Qi = [];
    ri = [];
    Li = [];
end

if K.l+K.f>0
    A =-(F_struc(1:K.f+K.l,2:end));
    B = full(F_struc(1:K.f+K.l,1));
end

if K.f > 0
    Aeq = A(1:K.f,:);
    beq = B(1:K.f);
else
    Aeq = [];
    beq = [];
end
if K.l > 0
    Aineq = A(1+K.f:K.f+K.l,:);
    bineq = B(1+K.f:K.f+K.l);
else
    Aineq = [];
    bineq = [];
end


if nnz(Q)==0
    H = [];
else
    H = full(2*Q);
end

NegativeSemiVar = [];
if ~isempty(semicont_variables)
    NegativeSemiVar = find(LB(semicont_variables) < 0);
    if ~isempty(NegativeSemiVar)
        temp = UB(semicont_variables(NegativeSemiVar));
        UB(semicont_variables(NegativeSemiVar)) = -LB(semicont_variables(NegativeSemiVar));
        LB(semicont_variables(NegativeSemiVar)) = -temp;
        if ~isempty(Aineq)
            Aineq(:,semicont_variables(NegativeSemiVar)) = -Aineq(:,semicont_variables(NegativeSemiVar));
        end
        if ~isempty(Aeq)
            Aeq(:,semicont_variables(NegativeSemiVar)) = -Aeq(:,semicont_variables(NegativeSemiVar));
        end
        H(:,semicont_variables(NegativeSemiVar)) = -H(:,semicont_variables(NegativeSemiVar));
        H(semicont_variables(NegativeSemiVar),:) = -H(semicont_variables(NegativeSemiVar),:);
        for i = 1:length(Qi)
            Qi{i}(:,semicont_variables(NegativeSemiVar)) = -Qi{i}(:,semicont_variables(NegativeSemiVar));
            Qi{i}(semicont_variables(NegativeSemiVar),:) = -Qi{i}(semicont_variables(NegativeSemiVar),:);
        end
        c(semicont_variables(NegativeSemiVar)) = -c(semicont_variables(NegativeSemiVar));
    end
end

% Bug in cplex 12.1
if length(K.q) == 1
    if isa(Qi,'cell')    
        Qi = Qi{1};
    end
end    

ctype = char(ones(length(c),1)*67);
ctype(setdiff(integer_variables,semicont_variables)) = 'I';
ctype(binary_variables)  = 'B';  % Should not happen except from bmibnb
ctype(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
ctype(intersect(semicont_variables,integer_variables)) = 'N';

options.cplex.simplex.display = options.verbose;
options.cplex.mip.display = options.verbose;
options.cplex.barrier.display = options.verbose;
if str2num(interfacedata.solver.subversion)>=12.3
    if options.verbose
        options.cplex.diagnostics = 'on';
    else
        options.cplex.diagnostics = 'off';
    end
end

if ~isempty(K.sos.type)
    for i = 1:length(K.sos.weight)
        K.sos.weight{i} = full(K.sos.weight{i}(:));
    end
    if length(K.sos.weight)==1
        K.sos.weight = K.sos.weight{1};
        K.sos.variables = K.sos.variables{1};
    end
end

if options.savedebug
    save cplexintdebug H c Aineq bineq Aeq beq Li Qi ri LB UB ctype
end

% Call mex-interface
solvertime = clock;
if isempty(integer_variables) & isempty(binary_variables) & isempty(semicont_variables) & isempty(K.sos.type)
    [x,fval,exitflag,output] = cplexqcp(H, c(:), Aineq,bineq,Aeq,beq,Li,Qi,ri,LB,UB,x0,options.cplex);    
else   
    [x,fval,exitflag,output] = cplexmiqcp(H, c(:), Aineq,bineq,Aeq,beq,Li,Qi,ri,K.sos.type,K.sos.variables,K.sos.weight,LB,UB,ctype',x0,options.cplex);   
end
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

if length(x) == length(c)
    if ~isempty(NegativeSemiVar)
       x(NegativeSemiVar) = -x(NegativeSemiVar);
    end
end
 
if isempty(x)
    x = zeros(n_original,1);
else
    x = x(1:n_original);
end

problem = 0;
D_struc = [];

% Check, currently not exhaustive...
switch output.cplexstatus
    case {1,101,102}
        problem = 0;
    case {3,103,106}
        problem = 1; % Infeasible
    case {2,20,21,118}
        problem = 2; % Unbounded
    case 4
        problem = 1;
    case {10,11,104,105,107,108,111,112}
        problem = 3; % Iteration/time
    case {5,6,109,110}
        problem = 4; % Numerics
    case 119
        problem = 15;        
    otherwise
        problem = -1;
end

infostr = yalmiperror(problem,'CPLEX-IBM');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.H = H;
    solverinput.f = c(:);
    solverinput.Aineq = Aineq;
    solverinput.bineq = bineq;
    solverinput.Aeq = Aeq;
    solverinput.beq = beq;
    solverinput.ctype = ctype;
    solverinput.LB = LB;
    solverinput.UB = UB;
    solverinput.x0 = [];
    solverinput.options = options.cplex;
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

function nSOCP=countSOCP(F_struc,K)
nSOCP = 0;
top = K.f+K.l + 1;
ri = zeros(1,length(K.q));
Li = [];
for i = 1:length(K.q)
    % [cx+d;Ax+b]   |Ax+b|<cx+d, originally a QCQP
    m = K.q(i);
    ci = F_struc(top,2:end)';
    di = F_struc(top,1);
    Ai = F_struc(top+1:top+m-1,2:end);
    if(min(eig(full(Ai'*Ai - ci*ci')))<0)
        nSOCP = nSOCP + 1;
    end
    top = top+m;
end
