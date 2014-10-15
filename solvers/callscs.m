function output = callscs(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
ub      = model.ub;
lb      = model.lb;

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addbounds(F_struc,K,ub,lb);
end

data.A = -F_struc(:,2:end);
data.b = full(F_struc(:,1));
data.c = c;
cones = K;

% Now add exponential cone information
if ~isempty(model.evalMap)
    % First check that we have only exponentials
    exponentials = [];
    logarithms = [];
    isexponential = zeros(1,length(model.evalMap));
    for i = 1:length(model.evalMap)
        if isequal(model.evalMap{i}.fcn,'exp') 
            exponentials = [exponentials model.evalMap{i}.computes];
            isexponential(i) = 1;           
        elseif isequal(model.evalMap{i}.fcn,'log') 
            logarithms = [logarithms model.evalMap{i}.computes];          
        else
            % Standard interface, return solver not applicable
            output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
            return
        end       
    end
    % Check that all exp/log enter in a convex fashion
    if model.K.f > 0
       if nnz(data.A(1:K.f,exponentials))>0
          output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
             return
        end 
    end
    if model.K.l > 0
        if nnz(data.A(1+model.K.f:model.K.f+model.K.l,exponentials)<0) || nnz(data.A(1+model.K.f:model.K.f+model.K.l,logarithms)>0)
             output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
             return
        end
    end
    if sum(model.K.q) + sum(model.K.s) > 0
         if nnz(data.A(1+model.K.f+K.l:end,exponentials))>0 || nnz(data.A(1+model.K.f+K.l:end,logarithms))>0
             output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
             return
        end
    end
    % We have to create new "y"-variables in yexp(x/y)<=z
    % Add constraint y = 1
    m = length(exponentials) + length(logarithms);
    cones.f = cones.f + m;
    data.b = [ones(m,1);data.b];
    data.A = [spalloc(m,size(data.A,2),0) speye(m);data.A spalloc(size(data.A,1),m,0)];
    data.c = [c;zeros(m,1)];
    
    % Now describe all exponential cones
    cones.ep = m;
    for i = 1:m
        if isexponential(i)
            % exp(xv) <= xc
            % y*exp(x/y)<= z.  y new variable, xv the original variable, and xc the "computed"
            y = length(model.c)+i; %New variable appended after the old
            x = model.evalMap{i}.variableIndex;
            z = model.evalMap{i}.computes;
            data.A = [data.A;sparse([1;2;3],[x y z],-1,3,size(data.A,2))];
            data.b = [data.b;zeros(3,1)];          
        else
            % log(xv) >= xc i.e. xv >= exp(xc)
            y = length(model.c)+i; %New variable appended after the old
            z = model.evalMap{i}.variableIndex;
            x = model.evalMap{i}.computes;
            data.A = [data.A;sparse([1;2;3],[x y z],-1,3,size(data.A,2))];
            data.b = [data.b;zeros(3,1)];       
        end
    end
end

if options.savedebug
    save scsdebug data cones
end

if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end

solvertime = clock; 
problem = 0;  
params = options.scs;
params.verbose = options.verbose;

switch  model.solver.tag
    case 'scs-direct'
         [x_s,y_s,s,info] = scs_direct(data,cones,params);
    otherwise
        [x_s,y_s,s,info] = scs_indirect(data,cones,params);
end

% solvertime = cputime - solvertime;%etime(clock,solvertime
if model.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% Internal format. Only recover the original variables
Primal = x_s(1:length(model.c)); 
Dual   = y_s;

switch info.status
    case 'Solved'
        problem = 0;
    case 'Infeasible'
        problem = 1;
     case 'Unbounded'
        problem = 2;    
    otherwise
        status = 9;
end

infostr = yalmiperror(problem,model.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.data = data
    solverinput.cones = cones;
    solverinput.param = param;  
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x_s;
    solveroutput.y = y_s;
    solveroutput.s = s;
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface 
output.Primal      = Primal;
output.Dual        = Dual;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;