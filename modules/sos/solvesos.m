function [sol,m,Q,residuals,everything] = solvesos(F,obj,options,params,candidateMonomials)
%SOLVESOS Sum of squares decomposition
%
%    [sol,m,B,residual] = solvesos(F,h,options,params,monomials) is used
%    for finding SOS decompositions of polynomials.
%
%    The coefficients of the polynomials are assumed linear w.r.t a set of
%    decision variables params and polynomial with respect to a variable x.
%    
%    An extension with a nonlinear parameterization in params is possible.
%    Note though that this gives BMIs or PMIs, solvable (locally) only if
%    PENBMI is installed, or by specifying 'moment' as solver to try to
%    solve the nonconvex semidefinite programming problem using a
%    semidefinite relaxation based on moments.
%
%    The SOS problem can be formulated as
%
%              min h(params)
%
%       subject to  F(i) >(=) 0 or F(i) is SOS w.r.t x
%
%   INPUT
%    F         : SET object with SOS constrained polynomials and constraints on variables params
%    h         : scalar SDPVAR object (can be [])
%    options   : options structure obtained from SDPSETTINGS (can be [])
%    params    : SDPVAR object defining parametric variables (can be [])
%    monomials : SDPVAR object with user-specified monomials for decomposition (can be [])
%
%   OUTPUT
%    sol       : Solution diagnostic from SDP problem
%    v         : Cell with monomials used in decompositions
%    Q         : Cell with Gram matrices, p = v{i}'*Q{i}*v{i}, where p is the ith SOS polynomial in your model.
%    residuals : Mismatch between p and decompositions. Same values (modulo numerical issue) as checkset(find(is(F,'sos')))
%                Warning, these residuals are not computed on matrix sos-of-squares
%
%   EXAMPLE
%    x = sdpvar(1);solvesos((sos(x^4+x^3+1)));                    % Simple decompositions
%    x = sdpvar(1);t = sdpvar(1);solvesos((sos(x^4+x^3+1-t)),-t); % Lower bound by maximizing t
%
%   NOTES
%
%    Variables not part of params, but part of non-SOS constraints in F
%    or objective h will automatically be appended to the params list.
%
%    To extract SOS decomposition, use the command SOSD (or compute from use v and Q)
%
%    If the 5th input argument is used, no additional monomial reduction is
%    performed (Newton, inconstency, congruence). It is thus assumed that
%    the supplied candidate monomials constitute a sufficient basis.
%
%    The field options.sos can be used to tune the SOS-calculations. See HTML help for details
%
%     sos.model          - Kernel (1) or image (2) representation of SOS problem [0|1|2 (0, YALMIP decides)]
%     sos.newton         - Use Newton polytope to reduce size [0|1 (1)]
%     sos.congruence     - Block-diagonalize using congruence classes [0|1|2 (2)]
%     sos.scale          - Scale polynomial [0|1 (1)]
%     sos.numblkdg       - Try to perform a-posteriori block-diagonalization [real  (0)]
%     sos.inconsistent   - Remove diagonal-inconsistent monomials [0|1|2 (0)]
%     sos.clean          - Remove monomials with coefficients < clean [real > 0 (1e-4)]
%     sos.traceobj       - Minimize trace of Gram matrix in problems without objective function [0|1 (0)]
%     sos.extlp          - Extract simple translated LP cones when performing dualization [0|1 (1)]
%
% See also SOS, SOSD, SDPSETTINGS, SOLVEMOMENT, SDPVAR, SDISPLAY

%% Time YALMIP
yalmip_time = clock;

% ************************************************
%% Check #inputs
% ************************************************
if nargin<5
    candidateMonomials = [];
    if nargin<4
        params = [];
        if nargin<3
            options = sdpsettings;
            if nargin<2
                obj = [];
                if nargin<1
                    help solvesos
                    return
                end
            end
        end
    end
end

if isa(F,'constraint')
    F = lmi(F);
end

if isempty(options)
    options = sdpsettings;
else
    if ~isa(options,'struct')
        error('The third argument should be an options structure');
    end
end

% Lazy syntax (not official...)
if nargin==1 & isa(F,'sdpvar')
    F = (sos(F));
end

if ~isempty(options)
    if options.sos.numblkdg
        [sol,m,Q,residuals,everything] = solvesos_find_blocks(F,obj,options,params,candidateMonomials);
        return
    end
end

[F,obj,m,everything,modeltype] = compilesos(F,obj,options,params,candidateMonomials);

if isempty(everything)
    % compilesos has detected trivially infeasible
    sol.yalmipversion = yalmip('version');
    sol.matlabversion = version;
    sol.yalmiptime = etime(clock,yalmip_time);
    sol.solvertime = 0;         
    sol.info = yalmiperror(1,'solvesos compilation');
    sol.problem = 1;  
    return
end

p = everything.p;
normp = everything.normp;
sizep = everything.sizep;
BlockedQ = everything.BlockedQ;
BlockedA = everything.BlockedA;
BlockedN = everything.BlockedN;
Blockedx = everything.Blockedx;
Blockedvarchange = everything.Blockedvarchange;
Blockedb = everything.Blockedb;
ranks  = everything.ranks;
UncertainData = everything.UncertainData;
ParametricVariables  = everything.ParametricVariables;
sol = everything.sol;

% % **********************************************
% %% SOLVE SDP
% % **********************************************
if sol.problem == 0
    if options.verbose > 0
        disp(' ');
    end
    if all(isinf(ranks))
        % Standard case
        sol =  solvesdp(F,obj,options);
    else
        % We have to alter the problem slightly if there are rank
        % constraints on the decompositions
        sol =  solveranksos(F,obj,options,ranks,BlockedQ);      
    end
end

% **********************************************
%% Recover solution (and rescale model+solution)
% **********************************************
for constraint = 1:length(p)
    for i = 1:length(BlockedQ{constraint})
        doubleQ{constraint}{i} = normp(constraint)*double(BlockedQ{constraint}{i});
    end
    doubleb{constraint}=normp(constraint)*double(Blockedb{constraint});
end

% **********************************************
%% Minor post-process
% **********************************************
if all(sizep<=1)
    [doubleQ,residuals] = postprocesssos(BlockedA,doubleb,doubleQ,[],options);
else
    residuals = 0;
end

% **********************************************
%% Safety check for bad solvers...
% **********************************************
if any(residuals > 1e-3) & (sol.problem == 0) & options.verbose>0
    disp(' ');
    disp('-> Although the solver indicates no problems,')
    disp('-> the residuals in the problem are really bad.')
    disp('-> My guess: the problem is probably infeasible.')
    disp('-> Make sure to check how well your decomposition')
    disp('-> matches your polynomial (see manual)')
    disp('-> You can also try to change the option sos.model')
    disp('-> or use another SDP solver.')
    disp(' ');
end

% **********************************************
%% Confused users. Primal dual kernel image?...
% **********************************************
if options.verbose > 0
    if sol.problem == 1
        if modeltype == 1
            disp(' ')
            disp('-> Solver reported infeasible dual problem.')
            disp('-> Your SOS problem is probably unbounded (SOS is dualized).')
            sol.problem = 2;
            sol.info = strrep(sol.info,'Infeasible problem','Unbounded objective function');
        elseif modeltype == 2
            disp(' ')
            disp('-> Solver reported infeasible primal problem.')
            disp('-> Your SOS problem is probably infeasible.')
        end
    elseif sol.problem == 2
        if modeltype == 1
            disp(' ')
            disp('-> Solver reported unboundness of the dual problem.')
            disp('-> Your SOS problem is probably infeasible (SOS is dualized).')
            sol.problem = 1;
            sol.info = strrep(sol.info,'Unbounded objective function','Infeasible problem');
        elseif modeltype == 2
            disp(' ')
            disp('-> Solver reported unboundness of the primal problem.')
            disp('-> Your SOS problem is probably unbounded.')            
        end
    elseif sol.problem == 12
            disp(' ')
            disp('-> Solver reported unboundness or infeasibility of the primal problem.')
            disp('-> Your SOS problem is probably unbounded but to clarify you should')   
            disp('-> solve the problem again without an objective')            
    end
end

% **********************************************
%% De-block
% **********************************************
for constraint = 1:length(p)
    Qtemp = [];
    for i = 1:length(BlockedQ{constraint})
        Qtemp = blkdiag(Qtemp,doubleQ{constraint}{i});
    end
    Q{constraint} = full(Qtemp);
end

% **********************************************
%% Experimental code for declaring sparsity
% **********************************************
if options.sos.impsparse == 1
    somesmall = 0;
    for i = 1:length(BlockedQ)
        for j =  1:length(BlockedQ{i})
            small = find(abs(double(BlockedQ{i}{j}))<options.sos.sparsetol);
            nullThese{i}{j} = small;
            if ~isempty(small)
                somesmall = 1;
            end
        end
    end
    if somesmall
        [F,obj,BlockedQ,Primal_matrices,Free_variables] = create_kernelmodel(BlockedA,Blockedb,F_parametric,parobj,options,nullThese);
        sol =  solvesdp(F,obj,options);
        for constraint = 1:length(p)
            for i = 1:length(BlockedQ{constraint})
                doubleQ{constraint}{i} = normp(constraint)*double(BlockedQ{constraint}{i});
            end
            doubleb{constraint}=normp(constraint)*double(Blockedb{constraint});
        end

        % **********************************************
        %% Post-process
        % **********************************************
        [doubleQ,residuals] = postprocesssos(BlockedA,doubleb,doubleQ,nullThese,options);
        for constraint = 1:length(p)
            Qtemp = [];
            for i = 1:length(BlockedQ{constraint})
                Qtemp = blkdiag(Qtemp,doubleQ{constraint}{i});
            end
            Q{constraint} = Qtemp;
        end
    end
end

% *********************************************
%% EXTRACT DECOMPOSITION
% *********************************************
switch sol.problem
    case {0,1,2,3,4,5} % Well, it didn't f**k up completely at-least

        % *********************************************
        %% GENERATE MONOMIALS IN SOS-DECOMPOSITION
        % *********************************************
        for constraint = 1:length(p)

            if constraint > 1 && isequal(BlockedN{constraint},BlockedN{constraint-1}) && isequal(Blockedx{constraint},Blockedx{constraint-1}) && isequal(Blockedvarchange{constraint},Blockedvarchange{constraint-1}) && isequal(sizep(constraint),sizep(constraint-1))
                monoms{constraint} = monoms{constraint-1};
            else
                monoms{constraint} = [];
                totalN{constraint} = [];
                N = BlockedN{constraint};
                x = Blockedx{constraint};
                for i = 1:length(N)
                    % Original variable
                    for j = 1:size(N{i},1)
                        N{i}(j,:)=N{i}(j,:).*Blockedvarchange{constraint};
                    end
                    if isempty(N{i})
                        monoms{constraint} = [monoms{constraint};[]];
                    else
                        mi = kron(eye(sizep(constraint)),recovermonoms(N{i},x));
                        monoms{constraint} = [monoms{constraint};mi];
                    end
                end
                if isempty(monoms{constraint})
                    monoms{constraint}=1;
                end
            end

            % For small negative eigenvalues
            % this is a good quick'n'dirty approximation
            % Improve...use shifted eigenvalues and chol or what ever...
            if ~any(any(isnan(Q{constraint})))
                if isempty(Q{constraint})
                    Q{constraint}=0;
                    h{constraint}=0;
                else
                    usedVariables = find(any(Q{constraint},2));
                    if length(usedVariables)<length(Q{constraint})
                        Qpart = Q{constraint}(usedVariables,usedVariables);
                        [U,S,V]=svd(Qpart);
                        R = sqrt(S)*V';
                        h0 = R*monoms{constraint}(usedVariables);
                        if isa(h0,'sdpvar')
                            h{constraint} = clean(R*monoms{constraint}(usedVariables),options.sos.clean);
                            h{constraint} = h{constraint}(findelements(h{constraint}));
                        else
                            h{constraint} = h0;
                        end
                    else
                        [U,S,V]=svd(mid(Q{constraint}));
                        R = sqrt(S)*V';
                        h0 = R*monoms{constraint};

                        if isa(h0,'sdpvar')
                            h{constraint} = clean(R*monoms{constraint},options.sos.clean);
                            if isa(h{constraint},'sdpvar')
                               h{constraint} = h{constraint}(findelements(sum(h{constraint},2)),:);                                
                            end
                        else
                            h{constraint} = h0;
                        end
                    end
                end
                if isempty(ParametricVariables)
                    ParametricVariables = [];
                end
                setsos(p{constraint},h{constraint},ParametricVariables,Q{constraint},monoms{constraint});
            else
                if options.verbose>0;
                    if UncertainData
                        disp(' ');
                        disp('-> Only partial decomposition is returned (since you have uncertain data).');
                        disp('');
                    else
                        disp(' ');
                        disp('-> FAILURE : SOS decomposition not available.');
                        disp('-> The reason is probably that you are using a solver that does not deliver a dual (LMILAB)');
                        disp('-> Use sdsettings(''sos.model'',2) to circumvent this, or use another solver (SDPT3, SEDUMI,...)');
                        disp('');
                        disp('-> An alternative reason is that YALMIP detected infeasibility during the compilation phase.');
                    end
                end
            end
        end

        m = monoms;

    otherwise
        Q = [];
        m = [];
end

% Don't need these outside
yalmip('cleardual')

% Done with YALMIP, this is the time it took, minus solver
if ~isfield(sol,'solvertime')
    sol.solvertime = 0;
end

sol.yalmiptime = etime(clock,yalmip_time)-sol.solvertime;
