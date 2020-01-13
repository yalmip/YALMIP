function tests = test_gp_lakshminarayanan
tests = functiontests(localfunctions);

function test1(dummy)
% Architectural parameters
beta = 4.03 * 10 ^ -4;
alpha = 2.6 * 10 ^ -7;
tau = 2.02 * 10 ^ -7;

% Number of processors
Procs = 5;

% -------------  End of Parameters (constants) -----------------------

obj3d = 1:5;
obj2d = 1:5;
objstripMP = 1:5;
objstrip1P = 1:5;
N = 1:5;

MP_tis = 1:5;

v2d_tis = 1:5;
v2d_tjs = 1:5;

v3d_tis = 1:5;
v3d_tjs = 1:5;
v3d_tks = 1:5;

N_i = 1000;
N_j = 1000;
j = 1;
for n = N_i/10: N_i/10 : N_i
    Nk(j)  = n;
    j = j + 1;
end;
for n = 2*N_i : N_i : 10*N_i
    Nk(j) = n;
    j = j + 1;
end;

fprintf('Iteration ');
for i = 1:5

    N_k = Nk(i);

    % Call the 3D model
    [obj3d(i), v3d_tis(i), v3d_tjs(i), v3d_tks(i)] = f3D(alpha,beta,tau,Procs,N_i,N_j,N_k);
    
    % Call the 2D model
    [obj2d(i), v2d_tis(i), v2d_tjs(i)] = f2D_SemiOblique_PerPlane(alpha,beta,tau,Procs,N_i,N_j,N_k);
    
    % Call the Strip Baseline Multi Pass model
    [objstripMP(i), MP_tis(i)] = fStrip_MP(alpha,beta,tau,Procs,N_i,N_j,N_k);
    
    % Call the Strip Baseline One Pass model
    objstrip1P(i) = fStrip_1P(alpha,beta,tau,Procs,N_i,N_j,N_k);

    % Save N for plots
    N(i) = Nk(i);
    % print current iteration.
    fprintf('%g ',i);
    if  mod(i,25) == 0
        fprintf('\n');
    end;
end   
assert(abs(norm(obj3d) + norm(obj2d) + norm(objstripMP) + norm(objstrip1P)-1.742617616580068e+002) <= 1e-2);

function [obj_val, ti_val, tj_val, tk_val] = f3D( alpha, beta, tau, Procs, Ni, Nj, Nk)
% Computes the optimal tile sizes given architectural parameters and domain sizes.
%
% INPUT ARGUMENTS: (alpha, beta, tau, Procs, Ni, Nj, Nk)
% alpha -- time to compute an iteration
% beta -- time to transfer a word
% tau -- startup cost
% Procs -- number of processors
% Ni -- size of domain along dimension i
% Nj -- size of domain along dimension j
% Nk -- size of domain along dimension k
%
% RETURN VALUES: [obj_val, ti_val, tj_val, tk_val]
% obj_val -- Value of objective function at the optimizer
% ti_val, tj_val, tk_val -- the optimal tile sizes.


% Notes
% -----
% Linear processor array allocation
% Tile and skew along k
% projection along j for every plane
% do the planes sequentially


% -------------  Tile Variables ----------------------------

% ti, tj, and tk are tile sizes
ti = sdpvar(1,1);
tj = sdpvar(1,1);
tk = sdpvar(1,1);

% ------------- Execution Time  Model ----------------------------


% Latency count
Lambda = (1+ ti/tj) * (Procs-1);

% Time for computing one tile
TilePeriod = (alpha * ti * tj * tk) + (2 * tau * tj * tk) + ( 2 * beta) ;

% Number of planes
no_planes = Nk / tk ;

% Number of passes per plane
no_passes_per_plane = (1/Procs) * ( (Ni + tk) / ti );

% Number of tiles per pass per plane -- macro column
no_tiles_per_pass_per_plane = (Nj + ti + 2 * tk) / tj;

% Total running time
T = ( Lambda + ( no_planes * no_passes_per_plane * no_tiles_per_pass_per_plane ) ) * TilePeriod ;

% ------------- End of Execution Time  Model ----------------------------


% -------------  Constraints ----------------------------

% Lower bounds on the tile vars
F = (ti >= 1) ;
F = F + (tj >= 1) ;
F = F + (tk >= 1) ;

% Upper bounds on the tile vars
F = F + (ti <= ( Ni  / Procs)) ;
F = F + (tj <= Nj) ;
% Number of tiles per macro column is atleast 2
% F = F + ( (Nj+ti+2tk) < 2*tj); 
F = F + (tk <= Nk)  ;

% No idle time between passes
%F = F + set ( ((Procs-1)/Nj) * (ti+tj) <= 1)
%F = F + set ( ((Procs-1)/ (2*Nj)) * (ti+tj) < 1)
%F = F + set ( (Procs-1) * ( ti + tj ) * (1/(2*Nj)) < 1 );
F = F + (  (1/(Nj + 2*Nk) ) * ( (Procs-2) * ti + (Procs-1) * tj ) <= 1);
% Not sure this is correct
%F = F + ( tj <= ( Nj / (Procs-1) ) ) ; 

% solve an integer version of the problem
intConstraints = (integer(ti)) + (integer(tj)) + (integer(tk));
F = F + intConstraints;

% + ( (Procs-1) * tj <= Nj + 2*tk) ; 
% No idle time between passes

% -------------  End of Constraints -------------------------

%%%%%%% Obective function
obj = T;

%%%%%  Solve the optimization problem 

optimize(F,obj, sdpsettings('verbose',1,'solver','bnb','debug',1));

%%%%% Print the solutions 


% Return values
obj_val = value(obj);
ti_val = value(ti);
tj_val = value(tj);
tk_val = value(tk);


function [obj_val, ti_val, tj_val] = f2D_SemiOblique_PerPlane( alpha, beta, tau, Procs, Ni, Nj, Nk)
% Computes the optimal tile sizes given architectural parameters and domain sizes.
%
% INPUT ARGUMENTS: (alpha, beta, tau, Procs, Ni, Nj, Nk)
% alpha -- time to compute an iteration
% beta -- time to transfer a word
% tau -- startup cost
% Procs -- number of processors
% Ni -- size of domain along dimension i
% Nj -- size of domain along dimension j
% Nk -- size of domain along dimension k
%
% RETURN VALUES: [obj_val, ti_val, tj_val]
% obj_val -- Value of objective function at the optimizer
% ti_val, tj_val -- the optimal tile sizes.


% Notes
% -----
% In this model we skew the ij plane and do not skew or tile the 
% k dimension. We solve for ti and tj.

% -------------  Tile Variables ----------------------------

% ti and tj are tile sizes
ti = sdpvar(1,1);
tj = sdpvar(1,1);

% ------------- Execution Time  Model ----------------------------


% Latency count
Lambda = (1 + (ti/tj)) * (Procs-1);

% Time for computing one tile
TilePeriod = (alpha * ti * tj) + (2 * tau * tj) + (2 * beta) ;

% Number of planes
no_planes = Nk;

% Number of passes per plane
no_passes_per_plane = (1/Procs) * ( Ni / ti );

% Number of tiles per pass per plane -- macro column
no_tiles_per_pass_per_plane = (Nj + ti) / tj;

% Total running time
T = ( Lambda + ( no_planes * no_passes_per_plane * no_tiles_per_pass_per_plane ) ) * TilePeriod ;

% ------------- End of Execution Time  Model ----------------------------


% -------------  Constraints ----------------------------

% Lower bounds on the tile vars
F = (ti >= 1) ;

% To avoid idle time between execution of successive planes
% we need tj >= 2.
F = F + (tj >= 2 ) ;

% Upper bounds on the tile vars
F = F + (ti <= ( Ni  / Procs)) ;
F = F + (tj <= Nj) ;

% No idle time between passes
%F = F + (  ( (1/Nj) * tj * (Procs-1) * (ti+tj) ) <= 1 );
%F = F + (  ( (1/Nj) * (Procs-1) * (ti+tj) ) <= 1 );
F = F + (  (1/Nj) * ( (Procs-2) * ti + (Procs-1) * tj ) <= 1);

% solve an integer version of the problem
intConstraints = (integer(ti)) + (integer(tj));
F = F + intConstraints;

% -------------  End of Constraints -------------------------

%%%%%%% Obective function
obj = T;

%%%%%  Solve the optimization problem 

optimize(F,obj, sdpsettings('verbose',1,'solver','bnb'));

%%%%% Print the solutions 


% Return values
obj_val = value(obj);
ti_val = value(ti);
tj_val = value(tj);


function [obj_val, ti_val] = fStrip_MP(alpha, beta, tau, Procs, Ni, Nj, Nk)
% Computes the optimal tile sizes given architectural parameters and domain sizes.
%
% INPUT ARGUMENTS: (alpha, beta, tau, Procs, Ni, Nj, Nk)
% alpha -- time to compute an iteration
% beta -- time to transfer a word
% tau -- startup cost
% Procs -- number of processors
% Ni -- size of domain along dimension i
% Nj -- size of domain along dimension j
% Nk -- size of domain along dimension k
%
% RETURN VALUES: [obj_val, ti_val, tj_val]
% obj_val -- Value of objective function at the optimizer
% ti_val -- the optimal tile size.


% Notes
% -----
% Divide the ij plane into strips.  
% In this model we allow the strip width to vary, so that
% there may be multiple passes. 
% We solve to ti.

% -------------  Tile Variables ----------------------------

% ti and tj are tile sizes
ti = sdpvar(1,1);

% ------------- Execution Time  Model ----------------------------


% Latency count
Lambda = (Procs-1);

% Time for computing one tile
TilePeriod = (alpha * ti * Nj) + (2 * tau * Nj) + ( 2 * beta) ;

% Number of planes
no_planes = Nk;

% Number of passes per plane
no_passes_per_plane = (1/Procs) * ( Ni / ti );

% Total running time
T = ( Lambda + ( no_planes * no_passes_per_plane ) ) * TilePeriod ;

% ------------- End of Execution Time  Model ----------------------------


% -------------  Constraints ----------------------------

% Lower bounds on the tile vars
%F = (ti > 1) ;
% To avoid idle time between execution of successive planes
% we need ti >= 2.
F = (ti >= 2) ;

% Upper bounds on the tile vars
F = F + (ti <= ( Ni  / Procs)) ;

% No idle time between passes
%F = F + ( (Procs-1) <= Nk );

% solve an integer version of the problem
intConstraints = (integer(ti)) ;
F = F + intConstraints;

% -------------  End of Constraints -------------------------

%%%%%%% Obective function
obj = T;

%%%%%  Solve the optimization problem 

optimize(F,obj, sdpsettings('verbose',1,'solver','bnb'));

%%%%% Print the solutions 

%fprintf('Objective : %g \n', value(obj));
%fprintf('ti : %g \n', value(ti));

% Return values
obj_val = value(obj);
ti_val = value(ti);



function [obj_val] = fStrip_1P(alpha, beta, tau, Procs, Ni, Nj, Nk)
% Computes the running time given architectural parameters and domain sizes.
%
% INPUT ARGUMENTS: (alpha, beta, tau, Procs, Ni, Nj, Nk)
% alpha -- time to compute an iteration
% beta -- time to transfer a word
% tau -- startup cost
% Procs -- number of processors
% Ni -- size of domain along dimension i
% Nj -- size of domain along dimension j
% Nk -- size of domain along dimension k
%
% RETURN VALUES: [obj_val, ti_val, tj_val]
% obj_val -- the running time.


% Notes
% -----
% Divide the ij plane into strips and do the plane in one pass.
% In this model we DO NOT allow the strip width to vary, so that
% there is only one pass.


% ------------- Execution Time  Model ----------------------------


% Latency count
Lambda = (Procs-1);

% Time for computing one tile
TilePeriod = (alpha * Ni * Nj / Procs )  + (2* tau * Nj) + (2*beta) ;

% Number of planes
no_planes = Nk;

% Total running time
T = ( Lambda +  no_planes ) * TilePeriod ;

% ------------- End of Execution Time  Model ----------------------------

if (Ni / Procs) < 2 
    fprintf('>>>>>  fStrip_1P: Constraint Ni/Procs >= 2 violated\n');
end;

% Return values
obj_val = T;


