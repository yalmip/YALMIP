function output = callclp(interfacedata)

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

if any(K.l)
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
solvertime = tic;
[x,lambda,problem] = clp(2*Q,c,A,b,Aeq,beq,lb,ub,ops);%,interfacedata.integer_variables);
solvertime = toc(solvertime);

if options.saveduals
    D_struc = -lambda;    
else
    D_struc = [];
end

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
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);