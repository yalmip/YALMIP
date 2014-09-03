function output = mpcvx(p)
%MPCVX          Approximate multi-parametric programming
%
% MPCVX is never called by the user directly, but is called by
% YALMIP from SOLVESDP, by choosing the solver tag 'mpcvx' in sdpsettings
%
% The behaviour of MPCVX can be altered using the fields
% in the field 'mpcvx' in SDPSETTINGS
%
% Implementation of
% A. Bemporad and C. Filippi,
% An Algorithm for Approximate Multiparametric Convex Programming
% Computational Optimization and Applications Volume 35, Number 1 /
% September, 2006

% *************************************************************************
% INITIALIZE DIAGNOSTICS IN YALMIP
% *************************************************************************
mpsolvertime = clock;
showprogress('mpcvx started',p.options.showprogress);

% *************************************************************************
% Display-logics
% *************************************************************************
switch max(min(p.options.verbose,3),0)
    case 0
        p.options.bmibnb.verbose = 0;
    case 1
        p.options.bmibnb.verbose = 1;
        p.options.verbose = 0;
    case 2
        p.options.bmibnb.verbose = 2;
        p.options.verbose = 0;
    case 3
        p.options.bmibnb.verbose = 2;
        p.options.verbose = 1;
    otherwise
        p.options.bmibnb.verbose = 0;
        p.options.verbose = 0;
end

% *************************************************************************
% No reason to save
% *************************************************************************
p.options.saveduals = 0;

% *************************************************************************
% The solver for approximation bounds
% *************************************************************************
solver    = eval(['@' p.solver.lower.call]); 

% *************************************************************************
% Generate an exploration set by finding an inner approximation of feasible
% parametric space
% *************************************************************************
[THETA,problem] = ray_shoot(p,solver);

% *************************************************************************
% Calculate optimal solution (x) in each initial vertex of exploration set
% *************************************************************************
X   = [];
OBJ = [];
for i = 1:size(THETA,2)
    [x_i,obj_i] = solve_node(p,solver,THETA(:,i));
    X   = [X x_i];
    OBJ = [OBJ obj_i];
end

% *************************************************************************
% Partition initial set to simplicies
% *************************************************************************
node_keeper;
node_keeper(THETA,X,OBJ);
T = delaunayn(THETA',{'Qz','Qt'});

if p.options.mpcvx.plot
    figure
    hold on
    plot(polytope(THETA'));
end

% *************************************************************************
% Run approximate algorithm recursively on all simplcies
% *************************************************************************
optimal_simplicies = [];
for i = 1:size(T,1)
    optimal_simplicies = [optimal_simplicies mp_simplex(p,solver,T(i,:)')];
end
mpsolvertime = etime(clock,mpsolvertime);

% *************************************************************************
% Create output compatible with MPT
% *************************************************************************
Pn = polytope;
j = 1;
for i = 1:size(optimal_simplicies,2)
    [theta,x,obj] = node_keeper(optimal_simplicies(:,i));
    m = size(theta,1);
    Minv = inv([ones(1,m+1);theta]);
    Vbar = obj'*Minv;
    Pn = [Pn polytope(theta')];
    Gi{j} = x(find(~ismember(1:length(p.c),p.parametric_variables)),:)*Minv(:,1);
    Fi{j} = x(find(~ismember(1:length(p.c),p.parametric_variables)),:)*Minv(:,2:end);
    Bi{i} = Vbar(2:end);
    Ci{i} = Vbar(1);
    j = j + 1;
end

% Prepare for generating costs
n = length(p.c) - length(p.parametric_variables);
m = length(p.parametric_variables);
Matrices = yalmip2mpt(p);
Matrices.getback.S1 = eye(n);
Matrices.getback.S2 = zeros(n,m);
Matrices.getback.S3 = zeros(n,1);
[Fi,Gi,details] = mpt_project_back_equality(Matrices,Fi,Gi,[],Matrices);
[Fi,Gi] = mpt_select_rows(Fi,Gi,Matrices.requested_variables);
[Fi,Gi] = mpt_clean_optmizer(Fi,Gi);
 
details.Ai = [];
details.Bi = Bi;
details.Ci = Ci;
output.problem = problem;
output.Primal      = nan*ones(length(p.c),1);
output.Dual        = [];
output.Slack       = [];
output.infostr     = yalmiperror(output.problem,'MPCVX');
output.solverinput = 0;
output.solveroutput.model = mpt_appendmodel([],polytope(THETA'),Pn,Fi,Gi,details);
output.solveroutput.U = p.used_variables(Matrices.free_var);
output.solveroutput.x = p.used_variables(Matrices.param_var);

output.solvertime   = mpsolvertime;

function simplex_solution = mp_simplex(p,solver,theta_indicies)

% Get the vertices, the optimal solution in the vertices and the optimal
% cost in the vertices
[theta,x,obj] = node_keeper(theta_indicies);

% Parametric dimension
m = size(theta,1);

% Define the polytope defined by the vertices, and add this polytope as a
% constraint to the YALMIP numerical model
Minv = inv([ones(1,m+1);theta]);
M1 = Minv(:,1);
M2 = zeros(size(M1,1),length(p.c));
M2(:,p.parametric_variables) = Minv(:,2:end);
pbound = p;
pbound.F_struc = [p.F_struc(1:p.K.f,:);M1 M2;p.F_struc(p.K.f + 1:end,:)];
pbound.K.l = pbound.K.l + size(M1,1);

% Save original linear cost
c = pbound.c;

% Linear approximation of optimal cost (which is an upper bound)
Vbar = obj'*Minv;
c2 = zeros(length(pbound.c),1);
c2(pbound.parametric_variables) = -Vbar(2:end);

% This should be the difference between bound and true optima
pbound.c = pbound.c + c2;
pbound.f = pbound.f - Vbar(1);

% Find minimum (we have negated error so this is the maximum, i.e. the
% point where the linear approximation of the optimal cost is worst)
output = feval(solver,pbound);

if p.options.mpcvx.plot
    plot(polytope(theta'),struct('color',rand(3,1)))
end

% Dig deeper?
upper = Vbar*[1;output.Primal(pbound.parametric_variables)];
lower = p.c'*output.Primal+output.Primal'*p.Q*output.Primal + p.f;
%eps_CP_S = -(pbound.c'*output.Primal+output.Primal'*pbound.Q*output.Primal + pbound.f);
eps_CP_S_abs = upper - lower;
eps_CP_S_rel = (upper - lower)/(1 + abs(upper) + abs(lower));
if eps_CP_S_abs > p.options.mpcvx.absgaptol & eps_CP_S_rel > p.options.mpcvx.relgaptol
    
    % Put back original cost to the model
    pbound.c = p.c;
    pbound.f = p.f;

    % Worst-case approximation point
    thetac = output.Primal(pbound.parametric_variables);

    % Optimizer in this point
    [x_i,obj_i] = solve_node(pbound,solver,thetac);
    
    % Add to list of vertices
    new_index = node_keeper(thetac(:),x_i(:),obj_i);

    % Partition current simplex and solve recursively in each new simplex
    % THe simplex is split such that the worst/case point is a vertex in
    % every new simplex
    simplex_solution = [];
    test = [];
  %  [xc,R] = facetcircle(polytope(theta'),1:3);   
    if 0%min(svd([ones(1,size(theta,2));theta])) < 0.2
        % Simplicies are getting too elongated. Add a random cut
        
        P1 = polytope(theta');
      %  xc = get(P1,'xCheb');
      %  xc = (0*xc + (1/3)*sum(theta,2))/1;
        [xd,R] = facetcircle(polytope(theta'),1:3);
        [i,j] = max(R);
        [H,K] = double(P1);
        x = sdpvar(2,1);
        Pnew1 =  polytope((H*x <= K) + (null(H(j,:))'*(x-xd(:,j)) >= 0));
        Pnew2 =  polytope((H*x <= K) + (null(H(j,:))'*(x-xd(:,j)) <= 0));
        plot(Pnew1);
        plot(Pnew2);
        theta1 = extreme(Pnew1)';
        theta2 = extreme(Pnew2)';
        X1 = [];X2 = [];
        OBJ1 = [];OBJ2 = [];
        indicies1 = []
        indicies2 = []
        for i = 1:size(theta1,2)
            [x_i,obj_i] = solve_node(p,solver,theta1(:,i));
            indicies1 = [indicies1  node_keeper(theta1(:,i),x_i(:),obj_i)];
            X1   = [X1 x_i];
            OBJ1 = [OBJ1 obj_i];
        end

        for i = 1:size(theta2,2)
            [x_i,obj_i] = solve_node(p,solver,theta2(:,i));
            indicies2 = [indicies2  node_keeper(theta2(:,i),x_i(:),obj_i)];
            X2   = [X2 x_i];
            OBJ2 = [OBJ2 obj_i];
        end
        T1 = (delaunayn(theta1',{'Qz','Qt'}));
        T2 = (delaunayn(theta2',{'Qz','Qt'}));
        for i = 1:size(T1,1)
            index = indicies1(T1);
            simplex_solution =  [simplex_solution mp_simplex(p,solver,index(i,:))];
        end
        for i = 1:size(T2,1)
            index = indicies2(T2);
            simplex_solution =  [simplex_solution mp_simplex(p,solver,index(i,:))];
        end

    else
        
        for i = 1:(size(theta,1)+1)
            j = 1:(size(theta,1)+1);
            j(i)=[];
            theta_test = [theta(:,j) thetac];
            if min(svd([ones(1,size(theta_test,2));theta_test]))>1e-4            
                simplex_solution =  [simplex_solution mp_simplex(p,solver,[theta_indicies(j);new_index])];
            end
        end
    end
    
% 	if length(simplex_solution) > 0
%         [theta,x,obj] = node_keeper(simplex_solution);
%         if max(abs(x(p.requested_variables,:)-mean(x(p.requested_variables,:))) < 1e-5)
%              simplex_solution = theta_indicies(:);
%         end
%     else
%         simplex_solution = theta_indicies(:);        
%     end
else
    % This simplex constitutes a node, report back
    simplex_solution = theta_indicies(:);
end

function varargout = node_keeper(varargin)

persistent savedTHETA
persistent savedX
persistent savedOBJ

switch nargin
    case 0 % CLEAR
        savedTHETA = [];
        savedX     = [];
        savedOBJ   = [];
    case 3 % Append
        savedTHETA = [savedTHETA varargin{1}];
        savedX     = [savedX varargin{2}];
        savedOBJ   = [savedOBJ varargin{3}];
        varargout{1} = size(savedTHETA,2);
    case 1 % Get data
        varargout{1} = savedTHETA(:,varargin{1});
        varargout{2} = savedX(:,varargin{1});
        varargout{3} = savedOBJ(:,varargin{1});varargout{3} = varargout{3}(:);
    otherwise
        error('!')
end


function [THETA,problem] = ray_shoot(p,solver)
THETA = [];

p_temp = p;
p_temp.c = p_temp.c*0;
p_temp.Q = 0*p_temp.Q;
n = length(p.parametric_variables);
for i = 1:eval(p.options.mpcvx.rays)
    p_temp.c(p.parametric_variables) = randn(length(p.parametric_variables),1);
    output = feval(solver,p_temp);
    THETA = [THETA output.Primal(p.parametric_variables)];
end
% Select unique and generate center
THETA = unique(fix(THETA'*1e4)/1e4,'rows')';
center = sum(THETA,2)/size(THETA,2);
THETA = [THETA*0.999+repmat(0.001*center,1,size(THETA,2))];
problem = 0;


function [x_i,obj_i] = solve_node(p,solver,theta);
p_temp = p;
p_temp.F_struc(:,1) = p_temp.F_struc(:,1) + p_temp.F_struc(:,1+p.parametric_variables)*theta;
p_temp.F_struc(:,1+p.parametric_variables) = [];

empty_rows = find(~any(p_temp.F_struc(p.K.f+1:p.K.f+p.K.l,2:end),2));
if ~isempty(empty_rows)
    if all(p_temp.F_struc(p.K.f+empty_rows,1)>=-1e-7)
        p_temp.F_struc(p.K.f+empty_rows,:)=[];
        p_temp.K.l = p_temp.K.l - length(empty_rows);
    else
        feasible = 0;
    end
end

x_var = find(~ismember(1:length(p.c),p.parametric_variables));
theta_var = p.parametric_variables;
Q11 = p.Q(x_var,x_var);
Q12 = p.Q(x_var,theta_var);
Q22 = p.Q(theta_var,theta_var);
c1 = p.c(x_var);
c2 = p.c(theta_var);

p_temp.Q = Q11;
p_temp.c = c1+2*Q12*theta;
p_temp.lb(theta_var) = [];
p_temp.ub(theta_var) = [];

%p_temp.c(p.parametric_variables) = [];
%p_temp.Q(:,p.parametric_variables) = [];
%p_temp.Q(p.parametric_variables,:) = [];
output = feval(solver,p_temp);

% Recover complete [x theta]
x_i = zeros(length(p.c),1);
x_i(find(~ismember(1:length(p.c),p.parametric_variables))) = output.Primal;
x_i(p.parametric_variables) = theta;
obj_i = x_i'*p.Q*x_i+p.c'*x_i;
