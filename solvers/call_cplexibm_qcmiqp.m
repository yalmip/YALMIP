function output = call_cplexibm_qcmiqp(interfacedata)

% Author Johan Löfberg

% This is a gateway to all CPLEX interfaces
% Call LP/QP solver if sufficient
if isempty(interfacedata.K.q) | interfacedata.K.q(1)==0
    output = call_cplexibm_miqp(interfacedata);
    return
end

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
Q       = interfacedata.f;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;

ub      = interfacedata.ub;
lb      = interfacedata.lb;

n_original = length(c);
[F_struc,K,lb,ub,semicont_variables] = extractSemiContBounds(F_struc,K,lb,ub,semicont_variables);
[F_struc,K,c,Q,ub,lb,Qi,Li,ri] = append_normalized_socp(F_struc,K,c,Q,ub,lb);

if nnz(Q)==0
    H = [];
else
    H = 2*Q;
end

if K.l+K.f>0
    A =-(F_struc(1:K.f+K.l,2:end));
    B = full(F_struc(1:K.f+K.l,1));
end

[NegativeSemiVar,H,c,A,lb,ub,semicont_variables] = negateNegativeSemiContVariables(H,c,A,lb,ub,semicont_variables,Qi);

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
    save cplexdebug H c Aineq bineq Aeq beq Li Qi ri lb ub ctype
end

% Call mex-interfac
showprogress('Calling CPLEX',options.showprogress);
solvertime = clock;
if isempty(integer_variables) & isempty(binary_variables) & isempty(semicont_variables) & isempty(K.sos.type)
    [x,fval,exitflag,output] = cplexqcp(H, c(:), Aineq,bineq,Aeq,beq,Li,Qi,ri,lb,ub,x0,options.cplex);    
else   
    [x,fval,exitflag,output] = cplexmiqcp(H, c(:), Aineq,bineq,Aeq,beq,Li,Qi,ri,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype',x0,options.cplex);   
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
    solverinput.LB = lb;
    solverinput.UB = ub;
    solverinput.x0 = [];
    solverinput.Qi = Qi;
    solverinput.Li = Li;
    solverinput.ri = ri;
    solverinput.options = options.cplex;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.output=output;
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
