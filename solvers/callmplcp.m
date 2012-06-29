function output = callmpt(interfacedata)

% Author Johan Löfberg
% $Id: callmplcp.m,v 1.1 2006-03-30 13:28:06 joloef Exp $

% Speeds up solving LPs in mpmilp
global mptOptions
if ~isstruct(mptOptions)
    mpt_error
end

% Convert to MPT
Matrices = yalmip2mpt(interfacedata);

% % Convert to colin
% parametric = find(~any(Matrices.G,2) & any(Matrices.E,2));
% nonparametric = find(~(~any(Matrices.G,2) & any(Matrices.E,2)));
% [M,Q,q,T,Tth,t] = lcp_mpqp(Matrices.H,Matrices.F',Matrices.Cf',Matrices.G(nonparametric,:),Matrices.E(nonparametric,:),Matrices.W(nonparametric,:));
% Ht = [Matrices.E(parametric,:) Matrices.W(parametric)];
% [BB,dPiv,dkPiv] = mplcp(M,q,Q,[Matrices.E(parametric,:) Matrices.W(parametric)]);

% 
% m = size(M,1);
% A  = [eye(m) -M];
% [At,bt] = a2s(Ht);
% bases = BB.bases;
% for i=[1:size(bases,1)]
%   iB = inv(A(:,bases(i,:)));
%   P(i) = polytope([At;-iB*Q],[bt;iB*q]);
% end;
% 
% model{1}.Pn = P;

%model = mpt_solvenode(Matrices,Matrices.lb,Matrices.ub,Matrices,[],options);


% Get some MPT options
options = interfacedata.options;
options.mpt.lpsolver = mptOptions.lpsolver;
options.mpt.milpsolver = mptOptions.milpsolver;
options.mpt.verbose = options.verbose;

if options.savedebug
    save mptdebug Matrices
end

if isempty(Matrices.binary_var_index)

    showprogress('Calling MPT',options.showprogress);
    solvertime = clock;
    if options.mp.presolve        
        [Matrices.lb,Matrices.ub] = mpt_detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);
    end        

    model = mpt_solvenode(Matrices,Matrices.lb,Matrices.ub,Matrices,[],options);
    solvertime = etime(clock,solvertime);

else  
    % Pre-solve required on binary problems
    options.mp.presolve = 1;

    solvertime = clock;    
    switch options.mp.algorithm
        case 1
            showprogress('Calling MPT via enumeration',options.showprogress);
            model = mpt_enumeration_mpmilp(Matrices,options);
        case 2
            showprogress('Calling MPT via parametric B&B',options.showprogress);
            model = bb_mpmilp(Matrices,Matrices,[],options)
        case 3
            showprogress('Calling MPT via delayed enumeration',options.showprogress);
            model = mpt_de_mpmilp(Matrices,options,[]);
            
        otherwise
    end
    solvertime = etime(clock,solvertime);
end

if isempty(model)
    model = {model};
end

if options.verbose
    if ~isempty(model{1})
    disp(['-> Generated ' num2str(length(model)) ' partitions.'])
    end
end


problem = 0;
infostr = yalmiperror(problem,'MPT');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.Matrices = Matrices;
    solverinput.options  = [];
else
    solverinput = [];
end

% Save all data from the solver?
% This always done
if options.savesolveroutput
    solveroutput.model = model;
    solveroutput.U = interfacedata.used_variables(Matrices.free_var);
    solveroutput.x = interfacedata.used_variables(Matrices.param_var);
else
    solveroutput = [];
end

% Standard interface
output.Primal      = nan*ones(length(interfacedata.c),1);
output.Dual        = [];
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;




