function Fchance = modelchance(F,options,rec)

% Goes through all probability constraints and checks for cases where we
% can use analytic expressions.

% Find chance constraints
if ~isempty(F)
    chanceDeclarations = find(is(F,'chance'));
else
    chanceDeclarations = [];
end
if isempty(chanceDeclarations)
    Fchance = F;
    return
end

% Find variables with attached distributions
randomDeclarations = find(is(F,'random'));
if isempty(randomDeclarations)
    error('Chance constraints without any distributions');
end

if nargin < 2
    options = sdpsettings;
end
if nargin < 3
    % Keep track of recursive calls for print-outs
    rec = 0;
end

keep = ones(length(F),1);
keep(chanceDeclarations)=0;
randomVariables = extractRandomDefinitions(F(randomDeclarations));
if options.verbose
    if ~rec
        disp('***** Starting YALMIP chance constraint module. *******************')
    else
        disp(' - (recursive application of chance constraints)')
    end
    disp([' - Detected ' num2str(length(randomVariables)) ' distribution models.'])
end

[randomVariables,map] = mergeDistributions(randomVariables);
if options.verbose && length(map)>max(map)
    disp([' - Merged to ' num2str(length(randomVariables)) ' distribution models.'])
end

groupedChanceConstraints = groupchanceconstraints(F);

if options.verbose
    disp([' - Detected ' num2str(length(groupedChanceConstraints)) ' chance constraints.'])
end

[Fchance,eliminatedConstraints,recursive] = deriveChanceModel(groupedChanceConstraints,randomVariables,options);
Fchance = Fchance + F(find(keep)) + F(find(keep(~eliminatedConstraints)));
if recursive
    Fchance = modelchance(Fchance,options,1);
end
if ~rec && options.verbose
    disp('***** Modeling of chance constraints done. ************************')
end


function [Fchance,eliminatedConstraints,recursive] = deriveChanceModel(groupedChanceConstraints,randomVariables,options);

recursive = 0;
Fchance = [];
eliminatedConstraints = zeros(length(groupedChanceConstraints),1);

allwVars = [];
for i = 1:length(randomVariables)
    allwVars = [allwVars;getvariables(randomVariables{i}.variables)];
end

for uncertaintyGroup = 1:length(randomVariables)
    
    wVars = getvariables(randomVariables{uncertaintyGroup}.variables);
    
    savedParameters = randomVariables{uncertaintyGroup}.distribution.parameters;
    nMixtures = length(randomVariables{uncertaintyGroup}.distribution.mixture);
       
    for ic = 1:length(groupedChanceConstraints)
        if length(groupedChanceConstraints{ic})>1
            error('Joint chance constraint not supported');
        end
        if ~is(groupedChanceConstraints{ic},'elementwise')
            error('Only elementwise chance constraints supported')
        end
        
        confidencelevel = struct(groupedChanceConstraints{ic}).clauses{1}.confidencelevel;
        gamma = 1-confidencelevel;              
        confidencelevelAllMixtures = confidencelevel;
        gammaAllMixtures = gamma;
        confidencelevelIndividualMixtures = sdpvar(1,nMixtures);
                                
        Xvec = sdpvar(groupedChanceConstraints{ic});
        
        for ix = 1:length(Xvec)
            X = Xvec(ix);
            
            if nMixtures > 1 && options.verbose
                disp([' - Expanding a mixture ']);
            end
                                    
            for mixtureComponent = 1:max(1,nMixtures)
                if nMixtures > 1
                    % Set distribution to current component
                    for i = 2:length(randomVariables{uncertaintyGroup}.distribution.parameters)
                        Param = savedParameters{i};
                        randomVariables{uncertaintyGroup}.distribution.parameters{i} = Param{mixtureComponent};
                    end
                    confidencelevel = confidencelevelIndividualMixtures(mixtureComponent);
                    gamma = 1-confidencelevel;                     
                end
                
                
                % Extract quadratic part, X = fX + X, where fx is other stuff
                [fX,X] = functionSeparation(X);
                
                if isa(fX,'sdpvar') && ~isempty(intersect(deepdepends(fX),wVars))
                    error('Stochastic uncertainty in nonlinear operator not supported yet.');
                end
                
                allVars = depends(X);
                if ~isempty(intersect(wVars,allVars))
                    xVars = setdiff(allVars,wVars);
                    x = recover(xVars);
                    w = recover(wVars);
                    
                    fail = 0;
                    [A,cx,b,cw,fail] = quadraticDecomposition(X,x,w);
                    
                    % Remap to original ordering on variables in distribution
                    % Base = full(getbase(randomVariables{uncertaintyGroup}.variables));
                    % Base = Base(:,2:end);
                    % [ii,jj,kk] = find(Base)
                    %  cw = cw*Base;
                    %  A = A*Base;
                    
                    % b(x) + c(x)'*w >= 0
                    if isempty(b)
                        b = 0;
                    end
                    b = b + fX;
                    if ~isempty(cx)
                        b = b + cx'*x;
                    end
                    c = cw';
                    if ~isempty(A)
                        c = c + A'*x;
                    end
                    
                    newConstraint = [];
                    if ~fail
                    %    confidencelevel = struct(groupedChanceConstraints{ic}).clauses{1}.confidencelevel;
                    %    gamma = 1-confidencelevel;
                        if strcmp(func2str(randomVariables{uncertaintyGroup}.distribution.generator),'random')
                            distName = randomVariables{uncertaintyGroup}.distribution.parameters{1};
                            switch distName
                                case 'dro'
                                    newConstraint = droChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                    printout(options.verbose,'dro',randomVariables{uncertaintyGroup}.distribution);
                                    eliminatedConstraints(ic)=1;
                                case 'moment'
                                    if isequal(options.chance.method,'momentchebyshev')
                                        newConstraint = momentChebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                    else
                                        newConstraint = momentChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                    end
                                    printout(options.verbose,'moment',randomVariables{uncertaintyGroup}.distribution);
                                    eliminatedConstraints(ic)=1;
                                case 'momentf'
                                    newConstraint = momentfactorizedChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                    printout(options.verbose,'factorized moment',randomVariables{uncertaintyGroup}.distribution);
                                    eliminatedConstraints(ic)=1;
                                case {'normal','normalm'}
                                    newConstraint = normalChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                    printout(options.verbose,'exact normal',randomVariables{uncertaintyGroup}.distribution);
                                    eliminatedConstraints(ic)=1;
                                case 'normalf'
                                    newConstraint = normalfactorizedChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                    printout(options.verbose,'exact normalf',randomVariables{uncertaintyGroup}.distribution);
                                    eliminatedConstraints(ic)=1;
                                case {'logistic', 'laplace','uniform'}
                                    newConstraint = symmetricUnivariateChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                    printout(options.verbose,['exact symmetric univariate'],randomVariables{uncertaintyGroup}.distribution);
                                    eliminatedConstraints(ic)=1;
                                case {'gamma','exponential','weibull','gamma','uniform'}
                                    newConstraint = nonsymmetricUnivariateChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                    printout(options.verbose,['exact nonsymmetric univariate'],randomVariables{uncertaintyGroup}.distribution);
                                    eliminatedConstraints(ic)=1;
                                otherwise
                                    switch options.chance.method
                                        case 'dro'
                                            newConstraint = droChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                        case {'chebyshev','chebychev'}
                                            newConstraint = sampledchebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                        case {'momentchebyshev','momentchebychev'}
                                            newConstraint = sampledmomentChebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                        case {'moment'}
                                            newConstraint = sampledmomentChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                        case 'markov'
                                            newConstraint =  sampledmarkovChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                        case 'chernoff'
                                            newConstraint =  sampledchernoffChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                        case 'integer'
                                            newConstraint =  sampledMIChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                        case 'scenario'
                                            newConstraint =  sampledScenarioChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                        otherwise
                                            error('Chance modeling approach not recognized');
                                    end
                                    printout(options.verbose,options.chance.method,randomVariables{uncertaintyGroup}.distribution);
                                    eliminatedConstraints(ic)=1;
                            end
                        else
                            switch options.chance.method
                                case 'chebyshev'
                                    newConstraint = sampledchebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                case 'moment'
                                    newConstraint = sampledmomentChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                case 'momentchebyshev'
                                    newConstraint = sampledmomentChebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                case 'markov'
                                    newConstraint =  sampledmarkovChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                case 'chernoff'
                                    newConstraint =  sampledchernoffChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                case 'integer'
                                    newConstraint =  sampledMIChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                otherwise
                                    error('Chance modeling approach not recognized');
                            end
                            printout(options.verbose,options.chance.method,randomVariables{uncertaintyGroup}.distribution);
                            eliminatedConstraints(ic)=1;
                        end
                    end
                    if ~isempty(newConstraint)
                        if ~isempty(intersect(depends(newConstraint),allwVars))
                            % New uncertainties popped up,i.e. parameters in a
                            % distribution, are distributions them selves
                            Fchance = [Fchance, probability(newConstraint)>=confidencelevel];
                            recursive = 1;
                        else
                            Fchance = [Fchance, newConstraint];
                        end
                    end
                end
            end
            if nMixtures > 1
                weights = randomVariables{uncertaintyGroup}.distribution.mixture;
                Fchance = [Fchance, sum(weights(:).*confidencelevelIndividualMixtures(:)) == confidencelevelAllMixtures];
            end
        end
    end
end
if any(eliminatedConstraints == 0)
    for weirdconstraints = find(eliminatedConstraints == 0)
        % This was listed as a probabilistic chance constraint, but there
        % appears to have been no random variables in the definition
        if options.verbose
            disp(' - Chance constraint with no random variables detected...')
        end
        confidencelevel = struct(groupedChanceConstraints{weirdconstraints}).clauses{1}.confidencelevel;
        gamma = 1-confidencelevel;
        Xvec = sdpvar(groupedChanceConstraints{weirdconstraints});
        % So basically X>=0, but if confidencelevel is <= 0 it can be removed
        if isa(confidencelevel,'double')
            % Simple case, something like X >= 0.5 hence it
            % must be satisfied, otherwise removed
            if confidencelevel > 0
                Fchance = [Fchance, Xvec >= 0];
            end
        elseif isa(confidencelevel,'sdpvar')
            % This is nasty. Probability is a decision variable, so this is
            % basically a combinatorial case. If probability > 0, then it must
            % be true
            binvar satisfied
            Fchance = [Fchance,  0 <= confidencelevel <= 1,
                implies(satisfied, [confidencelevel >= 1e-4, Xvec >=0])
                implies(1-satisfied, [confidencelevel <= 0])];
        end
        eliminatedConstraints(weirdconstraints) = 1;
    end
end

function [AAA,ccc,b,c_wTbase,fail] = quadraticDecomposition(X,x,w)
b = [];
A = [];
% Some pre-calc
xw = [x;w];
fail = 0;
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
[Qs,cs,fs,dummy,nonquadratic] = vecquaddecomp(X,xw);
c_wTbase = [];
AAA = [];
ccc = [];
for i = 1:length(X)
    Q = Qs{i};
    c = cs{i};
    f = fs{i};
    if nonquadratic
        error('Constraints can be at most quadratic, with the linear term uncertain');
    end
    Q_ww = Q(wind,wind);
    if nnz(Q_ww)>0
        fail = 1;
        return
    end
    Q_xw = Q(xind,wind);
    Q_xx = Q(xind,xind);
    c_x = c(xind);
    c_w = c(wind);
    
    %b = [b;f + c_w'*w];
    %A = [A;-c_x'-w'*2*Q_xw'];
    % A = [A -c_x-2*Q_xw*w];
    AAA = [AAA;sparse(2*Q_xw)];
    ccc = [ccc;sparse(c_x)];
    b = [b;f+x'*Q_xx*x];
    c_wTbase = [c_wTbase;c_w'];
end

function printout(verbose,method,distribution)

if verbose
    if strcmpi(func2str(distribution.generator),'random')
        disp([' - Using ''' method '''-filter on constraint with ''' distribution.parameters{1} ''' distribution']);
    else
        disp([' - Using ''' method '''-filter on constraint with data created by @' distribution.name']);
    end
end