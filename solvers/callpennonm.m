function output = callpennonm(model)

% Author Johan Löfberg
% $Id: callpennonm.m,v 1.7 2008-05-05 14:51:54 joloef Exp $


% Pull out SDP-specific data
norig = length(model.c);
if ~isequal(model.K.s,0)    
    sdp_data = model.F_struc(1+model.K.f+model.K.l:sum(model.K.s.^2)+model.K.f+model.K.l,:);    
    model.F_struc(1+model.K.f+model.K.l:sum(model.K.s.^2)+model.K.f+model.K.l,:) = [];
    % Loop through the SDP constraints and add the constraint S_i = F(y)
    start = 0;
    model.sdpvars = 0;
    for i = 1:length(model.K.s)
        ni = model.K.s(i);
        idx = find(triu(ones(ni)));
        ns = length(idx);
        triuFy = sdp_data(start+idx,:);start = start + ni^2;
        model.F_struc = [triuFy -eye(ns);model.F_struc zeros(size(model.F_struc,1),ns)];        
        model.sdpbasis{i} = getbase(sdpvar(ni));
        sdp_data = [sdp_data zeros(size(sdp_data,1),ns)];
        model.c = [model.c;zeros(ns,1)];
        model.Q = blkdiag(model.Q,zeros(ns));
        model.lb = [model.lb;-inf*ones(ns,1)];
        model.ub = [model.ub;inf*ones(ns,1)];
        model.monomtable = blkdiag(model.monomtable,eye(ns));
        model.variabletype = [model.variabletype zeros(1,ns)] ;       
    end
    model.K.f = sum((model.K.s).*(model.K.s+1)/2)+model.K.f;
else
    sdp_data = [];
end

% Build trees etc. This common format is used for fmincon, snopt and ipopt
% interfaces
model = yalmip2nonlinearsolver(model);

if ~model.derivative_available
    disp('Derivate-free call to PENNONM not yet implemented')
    error('Derivate-free call to PENNONM not yet implemented')
end
if model.options.savedebug
    save pennonmdebug model
end
showprogress('Calling PENNONM',model.options.showprogress);

n = length(model.linearindicies);
Infinity = 1.0E38;
pen.nvars = n;
pen.nlin = size(model.A,1)+size(model.Aeq,1);
pen.nconstr =1+size(model.Anonlineq,1)+size(model.Anonlinineq,1)+size(model.A,1)+size(model.Aeq,1);

pen.nnz_gradient = n;
pen.nnz_hessian = n^2;
pen.lbv = model.lb;pen.lbv(pen.lbv<-1e12) = -Infinity;
pen.ubv = model.ub;pen.ubv(pen.ubv> 1e12) =  Infinity;
pen.ubc = [repmat(0,length(model.bnonlinineq),1);
           repmat(0,length(model.bnonlineq),1);
           repmat(0,length(model.b),1);
           repmat(0,length(model.beq),1);0];
       
pen.lbc = [repmat(-Infinity,length(model.bnonlinineq),1);
         repmat(0,length(model.bnonlineq),1);
         repmat(-Infinity,length(model.b),1);
         repmat(0,length(model.beq),1);0];

pen.lbmv = [];
pen.ubmv = [];
pen.mnzs = [];
pen.mrow = [];
pen.mcol = [];

if ~isempty(sdp_data)
    pen.nsdp = length(model.K.s);
    pen.mtype = 0*model.K.s;
    pen.blks = model.K.s;
    pen.mnzs = (model.K.s).*(model.K.s+1)/2;
    pen.lbmv = zeros(length(model.K.s),1);
    pen.ubmv = zeros(length(model.K.s),1)+Infinity;
else
    pen.nsdp = [0];
    pen.mtype = [];
    pen.blks = [];
    pen.lbmv = [];
    pen.ubmv = [];
end

pen.xinit=model.x0;
pen.my_f = 'pennonm_callback_f';
pen.my_f_gradient = 'pennonm_callback_df';
pen.my_g = 'pennonm_callback_g';
pen.my_g_gradient = 'pennonm_callback_dg';
pen.ioptions = [100 100 2 0 0 0 1 0 0 1 0 0 0 -1 0 1 0 0]; 
pen.doptions = [1.0E-2 1.0E0 1.0E-0 1.0E-2 5.0E-1 5.0E-1 1.0E-6 1.0E-12 1.0e-7 0.05 1.0 1.0 1.0];
switch model.options.verbose
    case 0
        pen.ioptions(3) = 0;
    case 1 
        pen.ioptions(3) = 2;
    otherwise
        pen.ioptions(3) = 3;
end
if model.options.usex0
    pen.ioptions(8) = 0;
else
    pen.ioptions(8) = 1;
end

% YALMIP does not supply a Hessian
pen.ioptions(12) = 2;
pen.ioptions(13) = 2;
pen.my_f_hessian = '';
pen.my_g_hessian = '';

% Initialize the call-backs with a persistent model structure
pennonm_callback_f([],model);
pennonm_callback_df([],model);
pennonm_callback_g([],[],model);
pennonm_callback_dg([],[],model);

% These are needed to avoid recomputation due to pennons scalarized
% computation scheme in all-backs
global latest_x
global latest_x_g
global latest_df
global latest_G
global latest_g
latest_G= [];
latest_g = [];
latest_df = [];
latest_f = [];
latest_x = [];
latest_x_g = [];

% Solve
solvertime = clock;
[f,xout,u,inform,iresults,dresults] = pennonm(pen);
if ~isempty(sdp_data)
    % remove the S-variables
    xout(end-sum((model.K.s).*(model.K.s+1)/2)+1:end) = [];
    model.linearindicies(model.linearindicies>norig)=[];
    model.c = model.c(1:norig);
end
solvertime = etime(clock,solvertime);

x = RecoverNonlinearSolverSolution(model,xout);

% We currently don't extract dual solutions
D_struc = [];

% Check, currently not exhaustive...
switch inform
    case {0}
        problem = 0;
    case {2,3}
        problem = 1;
    case 4
        problem = 2;
    case 5
        problem = 3;
    case {1,7,6}
        problem = 4;
    otherwise       
        problem = 11;
end

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.model = model; 
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = x;
    solveroutput.F = F;
    solveroutput.xmul = xmul;
    solveroutput.inform=inform;
    solveroutput.Fstate=Fstate;
    solveroutput.ns = ns;
    solveroutput.ninf = ninf;
    solveroutput.sinf = sinf;
    solveroutput.mincw = mincw;
    solveroutput.miniw = miniw;
    solveroutput.miniw = miniw;
    
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'PENNON',solverinput,solveroutput,solvertime);