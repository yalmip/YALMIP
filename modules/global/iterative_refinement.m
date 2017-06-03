function output = iterative_refinement(model)
%ITERATIVE REFINEMENT metasolver
%
% ITERATVE_REFINEMENT is never called by the user directly, but is called by YALMIP from
% SOLVESDP, by choosing the solver tag 'refiner' in sdpsettings
%
% The behaviour of ITERATIVE_REFINEMENT can be altered using the fields in the field
% 'refiner' in SDPSETTINGS
%
% refiner.nbDigits           - Solver target precision in digits [real (35)]
% refiner.maxiter            - Maximum number of iterations [int (20)];
% refiner.verbose            - Verbosity level [ 0|1|2 (1)];
% refiner.internalSolver     - Internal solver [solver tag ('')]
% refiner.refinePrimal       - Whether the primal should be refined [true|false (true)]
% refiner.refineDual         - Whether the dual should be refined [true|false (true)]
% refiner.solvePrimalFirst   - Whether the first iteration is done in primal or dual form [true|false (false)]
% refiner.primalInPrimalForm - Whether the primal refinements should be solved in primal or dual form [true|false (false)]
% refiner.dualInPrimalForm   - Whether the primal refinements should be solved in primal or dual form [true|false (false)]
%
% Note : In order to take achieve precisions of more than 15 digits, this
% solver requires the high precision library GEM to be in matlab's path.
% This library can be downloaded from https://github.com/jdbancal/gem/releases
% Best results will be obtained by specifying the problem to be solved in
% terms of gem or sgem variables.
%
% Example:

% Author Jean-Daniel Bancal


showprogress('Iterative refinement started',model.options.showprogress);

% *************************************************************************
% INITIALIZE DIAGNOSTICS AND DISPLAY LOGICS (LOWER VERBOSE IN SOLVERS)
% *************************************************************************

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
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

if options.savedebug
    save refinerdebug F_struc c K options
end


% *********************************************
% Some general checks
% *********************************************

if (options.verbose >= 2) && (options.refiner.nbDigits > 15) && ((~isa(F_struc, 'gem') && ~isa(F_struc, 'sgem')) || (~isa(c, 'gem') && ~isa(c, 'sgem')))
    disp(' ');
    disp('Warning: Trying to solve a problem in high precision but all coefficients don''t have a high precision.');
    disp('         Additional zeros will be added at the end of these coefficients.');
    disp(' ');
end

% *********************************************
% Call the Refiner
% *********************************************
if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end

problem = 0;
warnState = warning;
solvertime = tic;

[fval x_s y_s z_s info] = refiner(-F_struc(:,2:end), -c, F_struc(:,1), K, options);
y_s = -y_s;

warning(warnState);
solvertime = toc(solvertime);

% Internal format
Primal = y_s;
Dual   = x_s;

infostr = yalmiperror(info.problem,model.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.A = -F_struc(:,2:end);
    solverinput.c = F_struc(:,1);
    solverinput.b = -c;
    solverinput.K = K;
    solverinput.options = options;
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x_s;
    solveroutput.y = y_s;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);


return;


end



function [fval x y z info] = refiner(Aeq, beq, f, K, options, l, fd, ld, beqd, Ld, xEst, yEst, zEst, beq0, l0, beqd0, Ld0, zoom)
% [fval x y z info] = refiner(Aeq, beq, f, K, options, l, fd, ld, beqd, Ld, xEst, yEst, zEst, beq0, l0, beqd0, Ld0, zoom)
%
% Tries to solve the semi-definite program
%     maximize   f'*x
%   subject to   Aeq*x = beq
%          and   x belongs to some self-dual homogeneous cones
%
% x can be a combination of scalars and matrices of various sizes, as 
% specified by K.f, K.l and K.s (a la sedumi).
%
% For better results, the above parameters should be provided as gem/sgem
% objects.
%
% l does not need to be specified. If specified, it sets a lower bound on
% the bounded variables (as specified by K.l).
% 
% fd, ld, beqd, Ld are specifications of the dual (they should NOT be
% provided by the user; they will be computed automatically).
%
% xEst, yEst, zEst are estimations of the primal and dual solutions. They 
% should NOT be provided by the user.
%
% beq0, l0, beqd0, Ld0 are the initial values of beq, l, beqd and Ld. They
% should NOT be provided manually.
%
% zoom is also an internal parameter.
%
% The output fval is the result of the optimization
% x is the optimal primal variable
% y and z are the optimal dual variables
% dimacs reports the standard errors associated with the solution

if false
    %As an example, doubleLP_gem(Aeq, beq, f, K) with
    Aeq = [1 0 1
          -1 0 0
           0 -1 0
           0 -1 0
           0 0 -1]
    beq = [0 -1 0];
    f = [-1 0 0 0 0]';
    %  l = [0]'; % <--- Is this constraint really needed?
    K.f = 1;     % <--- Take this into account (!!!)
    K.l = 0;
    K.s = 2;
    %solves the following problem:
    sedumi(Aeq, beq, -f, K)
    %or similarly :
    M=sdpvar(2)
    F = [trace(M) == 1, M >= 0]
    obj = M(1,2)
    optimize(F, obj)

    %Another example, which is a linear program:
    Aeq = [1 0 1
          -1 0 0
           0 -1 0
           0 -1 0
           0 0 -1]
    beq = [0 -1 0]';
    f = [-1 0 0 0 0]';
    %  l = [0 0 0 0]; % <--- Is this constraint really needed?
    K.f = 1;     % <--- Take this into account (!!!)
    K.l = 4;
    K.s = 0;
    % solves the following problem:
    sedumi(Aeq, beq, -f, K)
end

% Written by Jean-Daniel Bancal on 28 January 2016

%% Parameters setting

persistent nbIter ops highPrecisionSupported allDimacsTot


%% Arguments analysis and setup

% First of all, we check if high precision numbers are supported
if exist('gem','class') == 8
    highPrecisionSupported = true;
    
    % We setup shortcuts to the following functions
    gemify_f = @(x) gemify(x);
    toStrings_f = @(x,y) toStrings(x,y);
else
    highPrecisionSupported = false;

    % We emulate the missing functions
    gemify_f = @(x) x;
    toStrings_f = @(x,y) num2str(x,y);
end

% Arguments analysis
if nargin < 18
    % first call
    disp('Refiner 1.0 - Iterative meta-solver');
    tic;
    
    % We make sure the input is of high precision (if possible)
    if highPrecisionSupported
        if ~isa(Aeq, 'gem') && ~isa(Aeq, 'sgem')
            Aeq = gemify_f(Aeq);
        end
        if ~isa(beq, 'gem') && ~isa(beq, 'sgem')
            beq = gemify_f(beq);
        end
        if ~isa(f, 'gem') && ~isa(f, 'sgem')
            f = gemify_f(f);
        end
    end
    
    % Initialization of the iteration process
    nbIter = 1;
    if highPrecisionSupported
        zoom = gem(1);
    else
        zoom = 1;
    end

    % We set the initial values of the arguments
    if size(Aeq,2)==length(beq)
        Aeq = Aeq.';
    end
    
    if highPrecisionSupported
        l = gem(zeros(K.l,1));
    else
        l = zeros(K.l,1);
    end

    % We also prepare the specification of the dual program
    fd = beq;
    ld = l;
    beqd = -f;
    if highPrecisionSupported
        Ld = gem(zeros(size(l)));
    else
        Ld = zeros(size(l));
    end
    
    % We keep a copy of the initial arguments
    beq0 = beq;
    l0 = l;
    beqd0 = beqd;
    Ld0 = Ld;
    
    % We will monitor the precision improvement
    allDimacsTot = zeros(1,6);
else
    nbIter = nbIter + 1;
end


if (nargin >= 18) && (zoom < 0)
    % This is going to be the last iteration
    zoom = -zoom;
    onlyOnce = 1;
else
    onlyOnce = 0;
end


% We set the solver's parameters
if nargin < 5
    error('Not enough parameters');
else
    % We recover the parameters from the provided options
    maxNbIter = options.refiner.maxiter;
    verbose = options.verbose;
    refinePrimal = options.refiner.refinePrimal;
    refineDual = options.refiner.refineDual;
    solvePrimalFirst = options.refiner.solvePrimalFirst;
    primalInPrimalForm = options.refiner.primalInPrimalForm;
    dualInPrimalForm = options.refiner.dualInPrimalForm;
    if highPrecisionSupported
        precision = 10^gem(-options.refiner.nbDigits);
    else
        precision = 10^(-options.refiner.nbDigits);
        if (precision < 1e-15) 
            precision = 1e-15;
            if (nbIter == 1) && (verbose >= 1)
                warning('No high precision library found. Precision will be limitted to 1e-15.');
            end
        end
    end
    
    % Here are the options for the internal solver
    if nbIter == 1
        ops = options;
        ops.solver = options.refiner.internalSolver;
        ops.verbose = ops.verbose - 1;
%        ops.solver = 'sedumi';
%        ops.sedumi.eps = 1e-4;
        ops.dimacs = 1; % We ask for dimacs values
    end
end


% We make sure the gem default precision is larger than the required
% precision...
if highPrecisionSupported
    if gemWorkingPrecision < -log10(precision)+20
        if verbose >= 1
            warning(['Precision of the GEM library is low (', num2str(gemWorkingPrecision), ' digits), increasing it to ', num2str(-log10(precision)+20), ' digits.']);
        end
        gemWorkingPrecision(-log10(precision)+20);
    end
end




%% We solve the primal problem
if ((nbIter > 1) && refinePrimal) || ((nbIter == 1) && solvePrimalFirst)
    if primalInPrimalForm
        if verbose >= 2
            disp(' ');
            disp('-- Solving the primal refinement in primal form ');
            disp(' ');
        end
        % To solve the primal in primal form
        x = sdpvar(size(f,1), 1);
        obj = double(f)'*x;
        F = [double(Aeq)*(x+[zeros(K.f,1); double(l)]) == double(beq), x(K.f+[1:K.l]) >= 0];
        sol = solvesdp(F,obj,ops);

        if (sol.problem ~= 0) && (sol.problem ~= 3) && (sol.problem ~= 4) && (sol.problem ~= 5)
            % Problem is badly formulated, we come back to the previous
            % solution...
            nbIter = nbIter - 1;
            fval = [];
            x = value(x)+[zeros(1,K.f);l];
            y = dual(F(1));
            z = dual(F(2));

            % Invert unbounded and infeasible here
            if (sol.problem == 1) || (sol.problem == 2)
                sol.problem = 3-sol.problem;
            end
            
            info = sol;
            return;
        end
        
        % We save the result
        xi = gemify_f(double(x))+[zeros(1,K.f);l];
        yi = gemify_f(dual(F(1)));
        zi = gemify_f(dual(F(2)));
        
        dimacsp = sol.dimacs;
    else
        if verbose >= 2
            disp(' ');
            disp('-- Solving the primal refinement in dual form ');
            disp(' ');
        end
        % To solve the primal in dual form
        ypd = sdpvar(size(beq,1),1);
        zpd = sdpvar(size(l,1),1);
        objpd = double(beq)'*ypd - double(l)'*zpd;
        Fpd = [double(Aeq)'*ypd - [zeros(K.f,1); zpd] == -double(f), zpd >= 0];
        sol = solvesdp(Fpd,objpd,ops);
        
        if (sol.problem ~= 0) && (sol.problem ~= 3) && (sol.problem ~= 4) && (sol.problem ~= 5)
            % Problem is badly formulated, we come back to the previous
            % solution...
            nbIter = nbIter - 1;
            fval = [];
            x = -dual(Fpd(1));
            y = double(ypd);
            z = double(zpd);
            info = sol;
            return;
        end
        
        % We save the result
        xi = gemify_f(-dual(Fpd(1)));
        yi = gemify_f(value(ypd));
        zi = gemify_f(value(zpd));
        
        dimacsp = sol.dimacs;
    end
else
    % We didn't solve the primal, so we keep neutral solutions
    xi = zeros(size(f));
    yi = zeros(size(beq));
    zi = zeros(size(f));
    dimacsp = NaN*ones(1,6);
end


%% Now we solve the dual independently
if ((nbIter > 1) && refineDual) || ((nbIter == 1) && ~solvePrimalFirst)
    if dualInPrimalForm
        if verbose >= 2
            disp(' ');
            disp('-- Solving the dual refinement in primal form ');
            disp(' ');
        end
        % We solve the dual in primal form (so that it coincides with the
        % primal for the first iteration)

        xdd = sdpvar(size(beqd,1), 1);
        objdd = double(beqd)'*xdd + double(Ld)'*xdd(K.f+[1:K.l]);
        Fdd = [double(Aeq)*(xdd+[zeros(K.f,1);double(ld)]) == double(fd), xdd(K.f+[1:K.l]) >= 0];
        sold = solvesdp(Fdd,-objdd,ops);

        if (sold.problem ~= 0) && (sold.problem ~= 3) && (sold.problem ~= 4) && (sold.problem ~= 5)
            % Problem is badly formulated, we come back to the previous
            % solution...
            nbIter = nbIter - 1;
            fval = [];
            x = value(xdd)+[zeros(K.f,1);ld];
            y = dual(Fdd(1));
            z = dual(Fdd(2))+Ld;

            % Invert unbounded and infeasible here
            if (sold.problem == 1) || (sold.problem == 2)
                sold.problem = 3-sold.problem;
            end
            
            info = sold;
            return;
        end
        
        xdi = gemify_f(double(xdd))+[zeros(K.f,1);ld];
        ydi = gemify_f(dual(Fdd(1)));
        zdi = gemify_f(dual(Fdd(2))+Ld);
        
        dimacsd = sold.dimacs;
    else
        if verbose >= 2
            disp(' ');
            disp('-- Solving the dual refinement in dual form ');
            disp(' ');
        end
        % Solve the dual as a dual
        yd = sdpvar(size(fd,1),1);
        zd = sdpvar(size(ld,1),1);
        objd = double(fd)'*yd - double(ld)'*zd;
        Fd = [double(Aeq)'*yd - [zeros(K.f,1); (zd + double(Ld))] == double(beqd), zd >= 0];
        sold = solvesdp(Fd,objd,ops);

        if (sold.problem ~= 0) && (sold.problem ~= 3) && (sold.problem ~= 4) && (sold.problem ~= 5)
            % Problem is badly formulated, we come back to the previous
            % solution...
            nbIter = nbIter - 1;
            fval = [];
            x = -dual(Fd(1));
            y = value(yd);
            z = value(zd)+Ld;
            info = sold;
            return;
        end
        
        % We save the result
        ydi = gemify_f(double(yd));
        zdi = gemify_f(double(zd))+Ld;
        xdi = -gemify_f(dual(Fd(1)));
        %tdi = gemify_f(dual(Fd(2))); % This is basically the same as xdi...

        dimacsd = sold.dimacs;
    end
else
    % We didn't solve the dual, so we keep neutral solutions
    xdi = zeros(size(f));
    ydi = zeros(size(beq));
    zdi = zeros(size(f));
    dimacsd = NaN*ones(1,6);
end


%% We put both solutions together

% We keep the worst dimacs in memory:
dimacs = max(abs([dimacsp; dimacsd]));

if nargin < 13
    % Then we don't have a previous estimation of the solution. The best we
    % have is what we just obtained
    if solvePrimalFirst
        xTot = xi;
        yTot = yi;
        zTot = zi;

        % We need to assign the following variables, because the dual program
        % was not solved
        ydi = yi;
        zdi = zi;
    else
        xTot = xdi;
        yTot = ydi;
        zTot = zdi;

        % We need to assign the following variables, because the primal program
        % was not solved
        xi = xdi;
        yi = ydi;
        zi = zdi;
    end
else
    % Then we are not in the first iteration, so we take into account
    % all previous estimations...
    xTot = xEst + xi/zoom;
    yTot = yEst + ydi/zoom; % The dual variables are taken from the solution of the dual program
    zTot = zEst + zdi/zoom;
end

% Now we estimate the error for all rounds until now:
dimacsTot = computedimacs(beq0, f, Aeq, xTot, -yTot, [], K);
allDimacsTot(nbIter,:) = dimacsTot;

if highPrecisionSupported
    dimacsTotStr = toStrings_f(dimacsTot,5);
    for i=1:length(dimacsTotStr)
        dimacsTotStr{i} = [dimacsTotStr{i} '  '];
    end
    dimacsTotStr = [dimacsTotStr{:}];
else
    dimacsTotStr = num2str(dimacsTot);
end

% norm1 = @(x) max(abs(x));
% norm2 = @(x) sqrt(sum(x.^2));
% err1 = norm2(Aeq*xTot - beq0)/(1+norm1(beq0));
% err2 = max([0,-min(xTot(K.f+[1:K.l])-l0)/(1+norm1(beq0))]);
% err3 = norm2((yTot'*Aeq)' - [zeros(K.f,1);zTot] - beqd0)/(1+norm1(beqd0));
% err4 = max([0,-min(zTot-Ld0)/(1+norm1(beqd0))]);
% err5 = (-f'*xTot - (fd'*yTot - ld'*zTot))/(1 + abs(f'*xTot) + abs(fd'*yTot - ld'*zTot));
% err6 = (xTot'*[zeros(K.f,1);zTot])/(1 + abs(f'*xTot) + abs(fd'*yTot - ld'*zTot));
% dimacsTot = [err1 err2 err3 err4 err5 err6];


%% Diagnosis

% If we have enough precision, if we didn't gain enough precision at this
% run, or if we performed too many optimizations, we stop
if (nbIter >= maxNbIter) ...
    || (max(abs(dimacsTot)) < precision) ...
    || ((max(abs(dimacsTot(1:4))) < precision) && (max(abs(dimacsTot(5:6)))/max(abs(dimacsTot(1:4))) > 1e5)) ...
    || (max(abs(dimacs)) > 1e-1) ...
    || (onlyOnce == 1) ...
    || (~refineDual && (max(abs(dimacsTot(1:2))) < precision)) ...
    || (~refinePrimal && (max(abs(dimacsTot(3:4))) < precision)) ...
    || ((nbIter >= 6) && (max(abs(allDimacsTot(end,:))./mean(abs(allDimacsTot(end-5:end-1,:)))) > 1e-1))

    if verbose >= 1
        disp(['Iteration number ', num2str(nbIter), ', ', num2str(toc), 's. Current errors : ', dimacsTotStr, '  (', num2str(max(abs(dimacsp))), ', ', num2str(max(abs(dimacsd))), ') ']);
        disp(' ');
    end

    % We prepare the variables for the output
    if refinePrimal && refineDual
        % If both sides were optimized, we return both values
        fval = [-f'*xTot (fd'*yTot - ld'*zTot)];
    elseif refineDual
        % If only the dual was iterated, we return the dual objective
        % function
        fval = (fd'*yTot - ld'*zTot);
    else
        % By default we return the primal objective function
        fval = -f'*xTot;
    end
    x = xTot;
    y = yTot;
    z = zTot;

    % Let's give a short summary message
    if max(abs(dimacsTot)) <= precision
        info.problem = 0;
        if verbose >= 1
            disp(['The solver achieved a precision of ', toStrings_f(max(abs(dimacsTot)),5), ' in ', num2str(nbIter), ' iterations.']);
            disp(['Value of the objective function : ', toStrings_f(fval(1), -log10(precision))])
        end
    else
        if nbIter >= maxNbIter
            info.problem = 3;
            if verbose >= 1
                disp(['Iterative solver stopped after having reached the maximum number of rounds (', num2str(maxNbIter), ' iterations).']);
                disp(['The precision of the current solution is ', toStrings_f(max(abs(dimacsTot)),5)]);
                disp(['Current value of the objective function : ', toStrings_f(fval(1), -log10(precision))])
            end
        elseif (max(abs(dimacs)) > 1e-1)
            info.problem = 5;
            if verbose >= 1
                disp(['Iterative solver stopped because the precision of the ', num2str(nbIter), 'th iteration is only ', num2str(max(abs(dimacs)))]);
                disp(['The precision of the current solution is ', toStrings_f(max(abs(dimacsTot)),5)]);
                disp(['Current value of the objective function : ', toStrings_f(fval(1) -log10(precision))])
            end
        elseif ((nbIter >= 6) && (max(abs(allDimacsTot(end,:))./mean(abs(allDimacsTot(end-5:end-1,:)))) > 1e-1))
            info.problem = 5;
            if verbose >= 1
                disp(['Iterative solver stopped because the precision improvement over the last iterations is too small']);
                disp(['The precision of the current solution is ', toStrings_f(max(abs(dimacsTot)),5)]);
                disp(['Current value of the objective function : ', toStrings_f(fval(1), -log10(precision))])
            end
        elseif (refinePrimal && ~refineDual && max(abs(dimacsTot(1:2))) < precision)
            info.problem = 0;
            if verbose >= 1
                disp(['Primal precision of ', num2str(max(abs(dimacsTot(1:2)))), ' achieved in ', num2str(nbIter), ' iterations.']);
                disp(['Value of the objective function : ', toStrings_f(fval(1), -log10(precision))])
            end
        elseif (~refinePrimal && refineDual && max(abs(dimacsTot(3:4))) < precision)
            info.problem = 0;
            if verbose >= 1
                disp(['Dual precision of ', num2str(max(abs(dimacsTot(3:4)))), ' achieved in ', num2str(nbIter), ' iterations.']);
                disp(['Value of the objective function : ', toStrings_f(fval(1), -log10(precision))])
            end
        elseif (refinePrimal && refineDual && max(abs(dimacsTot(1:4))) < precision)
            info.problem = 4;
            if verbose >= 1
                disp(['The primal and dual solutions each have a precision of ', toStrings_f(max(abs(dimacsTot(1:4))),5), ' after ', num2str(nbIter), ' iterations.']);
                disp(['However, these solutions are only compatible up to ', toStrings_f(max(abs(dimacsTot(5:6))),5), ': the dual and primal values are respectively']);
                if highPrecisionSupported
                    display([fval(2); fval(1)], -log10(precision));
                else
                    num2str([fval(2); fval(1)], 15);
                end
            end
        else
            info.problem = 9;
        end
    end
    
    if verbose >= 1
        disp(' ');
    end
    
    return;
end


%% We prepare the next call

% Here is our renormalization factor. It defines how many digits from the
% primal we are ready to trust... There is at least an uncertainty of order
% epsilon.
if highPrecisionSupported
    zoom2 = gem(10^(floor(log10(1/max([eps abs(dimacs)])))-1));
else
    zoom2 = 10^(floor(log10(1/max([eps abs(dimacs)])))-1);
end

if verbose >= 1
    disp(['Iteration number ', num2str(nbIter), ', ', num2str(toc), 's. Current errors : ', dimacsTotStr, '  (', num2str(max(abs(dimacsp))), ', ', num2str(max(abs(dimacsd))), ') ']);
end

% If no refinement was requested, we prepare the output variables and exit.
if (refinePrimal ~= 1) && (refineDual ~= 1)
    fval = -f'*xTot;
    x = xTot;
    y = yTot;
    z = zTot;
    info.problem = 0;
    if max(abs(dimacsTot)) > precision
        info.problem = 4;
    end
    return;
end



factor = 1;
nbRescaling = 0;
zoom2 = zoom2*1;
fval = [];
while (nbRescaling < 10) && (zoom2 > 1) && isempty(fval)
    nbRescaling = nbRescaling + 1;
    if zoom2 >= 100
        zoom2 = zoom2/10;
    else
        zoom2 = zoom2/2;
    end
    
    % The new primal :
    f2 = f;
    beq2 = zoom2*(beq-Aeq*xi);
    l2 = max([-factor*ones(size(l)), zoom2*(l-xi(K.f+[1:K.l]))],[],2);

    % The new dual :
    fd2 = fd;
    ld2 = ld;
    beqd2 = zoom2*(beqd-(ydi'*Aeq)'+[zeros(K.f,1);zdi]);
    Ld2 = max([-factor*ones(size(Ld)), zoom2*(Ld-zdi)],[],2);

    % We solve these new problems
    [fval, x, y, z, info] = refiner(Aeq, beq2, f2, K, options, l2, fd2, ld2, beqd2, Ld2, xTot, yTot, zTot, beq0, l0, beqd0, Ld0, zoom*zoom2);

    % If the result is empty, this means that there was a problem in the 
    % following iteration. So we try to rescale by a smaller factor...
    % If this doesn't work for a few times, we abort
    if isempty(fval) && (verbose >= 1)
        disp(['Rescaled problem appears to be badly conditioned for zoom = 10^', num2str(double(log10(zoom2))), ', trying again with a smaller zoom factor.']);
    end
end


if isempty(fval)
    % If we arrive here, we encountered numerical problems
    x = xTot;
    y = yTot;
    z = zTot;
    info = [];
    info.problem = 4;
end

return;

end
