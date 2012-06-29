function output = callfmincon(model)

% Author Johan Löfberg
% $Id: callfminsearch.m,v 1.3 2007-08-07 06:55:57 joloef Exp $

% Retrieve needed data
options = model.options;
c       = model.c;
x0      = model.x0;
Q       = model.Q;
lb      = model.lb;
ub      = model.ub;
K       = model.K;
monomtable = model.monomtable;

nobounds = isempty(lb) | (~isempty(lb) & all(isinf(lb)));
nobounds = nobounds & (isempty(ub) | (~isempty(ub) & all(isinf(ub))));
if any(K.l) | any(K.q) | any(K.s) | ~nobounds
    output = createoutput(zeros(length(c),1),[],[],-4,'FMINSEARCH',[],[],0);
    return
end

if isempty(model.evaluation_scheme)
    model = build_recursive_scheme(model);
end

switch options.verbose
    case 0
        options.fmincon.Display = 'off';
    case 1
        options.fmincon.Display = 'final';
    otherwise
        options.fmincon.Display = 'iter';
end

% Do some pre-calc to be used in calls from fmincon
nonlinearindicies = union(find(model.variabletype~=0),model.evalVariables);
linearindicies    = setdiff(find(model.variabletype==0),nonlinearindicies);
model.nonlinearindicies = nonlinearindicies;
model.linearindicies    = linearindicies;

model.Anonlinineq = [];
model.bnonlinineq = [];
model.Anonlineq = [];
model.bnonlineq = [];

% Extract linear and nonlinear equality constraints
Aeq = [];
beq = [];
A = [];
b = [];

% This helps with robustness in bnb in some cases
x0candidate = zeros(length(c),1);

if isempty(x0)
    x0 = x0candidate(linearindicies);
else
    x0 = x0(linearindicies);
end

if options.savedebug
    ops = options.fminsearch;
    save fminsearch model x0 ops
end

showprogress('Calling FMINSEARCH',options.showprogress);

% Precalc for the callbacks
model = setup_fmincon_params(model);
if (model.SimpleQuadraticObjective | model.SimpleNonlinearObjective) & isempty(model.evalMap)
    options.fmincon.GradObj = 'on';    
end

solvertime = clock;
[xout,fmin,flag,output] = fminsearch('fmincon_fun',x0,options.fminsearch,model);
solvertime = etime(clock,solvertime);

if isempty(nonlinearindicies)
    x = xout(:);
else
    x = zeros(length(c),1);
    for i = 1:length(linearindicies)
        x(linearindicies(i)) = xout(i);
    end
    x = x(1:length(c));
end

problem = 0;

% Internal format for duals
D_struc = [];

% Check, currently not exhaustive...
if flag==0
    problem = 3;
else
    problem = 0;
end

% Save all data sent to solver?
if options.savesolverinput   
    solverinput.fun = 'fmincon_fun';
    solverinput.x0 = x0;
    solverinput.ops = options.fminsearch;
    solverinput.param = model;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'FMINSEARCH',solverinput,solveroutput,solvertime);





% % function output = callfminsearch(interfacedata)
% % 
% % % Author Johan Löfberg
% % % $Id: callfminsearch.m,v 1.3 2007-08-07 06:55:57 joloef Exp $
% % 
% % % Retrieve needed data
% % options = interfacedata.options;
% % F_struc = interfacedata.F_struc;
% % c       = interfacedata.c;
% % K       = interfacedata.K;
% % x0      = interfacedata.x0;
% % Q       = interfacedata.Q;
% % lb      = interfacedata.lb;
% % ub      = interfacedata.ub;
% % monomtable = interfacedata.monomtable;
% % 
% % switch options.verbose
% %     case 0
% %         options.fmincon.Display = 'off';
% %     case 1
% %         options.fmincon.Display = 'final';
% %     otherwise
% %         options.fmincon.Display = 'iter';
% % end
% % 
% % % Do some pre-calc to be used in calls from fmincon
% % temp = sum(interfacedata.monomtable,2)>1;
% % nonlinearindicies = find(temp);
% % linearindicies = find(~temp);
% % interfacedata.nonlinearindicies = nonlinearindicies;
% % interfacedata.linearindicies = linearindicies;
% % 
% % Aeq = [];
% % beq = [];
% % A = [];
% % b = [];
% % 
% % if isfield(options.fminsearch,'LargeScale')
% %     if isequal(options.fminsearch.LargeScale,'off')
% %         A = full(A);
% %         b = full(b);
% %         Aeq = full(Aeq);
% %         beq = full(beq);
% %     end
% % end
% % 
% % if isempty(x0)
% %     x0 = zeros(length(linearindicies),1);
% % else
% %     x0 = x0(linearindicies);
% % end
% % 
% % if options.savedebug
% %     ops = options.fminsearch;
% %     save fmincondebug interfacedata A b Aeq beq x0 ops
% % end
% % 
% % showprogress('Calling FMINSEARCH',options.showprogress);
% % fmincon_fun([],[],[],interfacedata); % initialize persistent
% % fmincon_con([],[],[],interfacedata); % initialize persistent
% % solvertime = clock;
% % [xout,fmin,flag,output] = fminsearch('fmincon_fun',x0,options.fminsearch);    
% % solvertime = etime(clock,solvertime);
% % 
% % if isempty(nonlinearindicies)
% %     x = xout(:);
% % else
% %     x = zeros(length(interfacedata.c),1);
% %     for i = 1:length(linearindicies)
% %         x(linearindicies(i)) = xout(i);
% %     end
% % end
% % 
% % problem = 0;
% % 
% % % Internal format for duals
% % D_struc = [];
% % 
% % % Check, currently not exhaustive...
% % if flag==0
% %     problem = 3;
% % else
% %     if flag>0
% %         problem = 0;
% %     else
% %         if isempty(x)
% %             x = repmat(nan,length(c),1);
% %         end
% %         if 0%any((A*x-b)>sqrt(eps)) | any( abs(Aeq*x-beq)>sqrt(eps))
% %             problem = 1; % Likely to be infeasible
% %         else
% %             if c'*x<-1e10 % Likely unbounded
% %                 problem = 2;
% %             else          % Probably convergence issues
% %                 problem = 5;
% %             end
% %         end
% %     end
% % end
% % infostr = yalmiperror(problem,'FMINSEARCH');
% % 
% % % Save all data sent to solver?
% % if options.savesolverinput
% %     solverinput.A = A;
% %     solverinput.b = b;
% %     solverinput.Aeq = Aq;
% %     solverinput.beq = beq;
% %     solverinput.c = c;
% %     solverinput.H = Q;
% %     solverinput.options = options.quadprog;
% % else
% %     solverinput = [];
% % end
% % 
% % % Save all data from the solver?
% % if options.savesolveroutput
% %     solveroutput.x = x;
% %     solveroutput.fmin = fmin;
% %     solveroutput.flag = flag;
% %     solveroutput.output=output;
% %     solveroutput.lambda=lambda;
% % else
% %     solveroutput = [];
% % end
% % 
% % 
% % 
% % % Standard interface
% % output.Primal      = x;
% % output.Dual        = D_struc;
% % output.Slack       = [];
% % output.problem     = problem;
% % output.infostr     = infostr;
% % output.solverinput = solverinput;
% % output.solveroutput= solveroutput;
% % output.solvertime  = solvertime;