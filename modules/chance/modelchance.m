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
        disp('***** Starting YALMIP chance constraint module. ***************************')
    else
        disp(' - (recursive application of chance constraints)')
    end
    nw = length(randomVariables);
    disp([' - Detected ' num2str(nw) ' distribution declaration' pluralS(nw)])
end
% Merging distribution
% 1. Multiple with same distribution are placed as one vector-valued
% 2. The different versions of Gaussian are combined and normalized, and
% will have three properties, mean, covariance, and factorized covariance,
% and all always called 'normal' from now on
[randomVariables,map] = mergeDistributions(randomVariables);
randomVariables = setupCharacterstics(randomVariables);

if options.verbose && length(map)>max(map)
    nw = length(randomVariables);
    disp([' - Merged to ' num2str(nw) ' distribution model' pluralS(nw)])
end

groupedChanceConstraints = groupchanceconstraints(F);

if options.verbose
    nc = length(groupedChanceConstraints);
    disp([' - Detected ' num2str(nc) ' chance constraint' pluralS(nc)])
end

% Some strategies exploit simplex structure
Simplicies = F(find(is(F,'simplex')));
SimplexInfo = sparse([]);
if length(Simplicies)>0
    for i = 1:length(Simplicies)
        SimplexInfo(end+1,getvariables(Simplicies(i))) = 1;
    end
else
    Candidates = F(find(is(F,'equality')));
    for i = 1:length(Candidates)
        C = sdpvar(Candidates(i));
        B = getbase(C);
        if size(B,1)==1 && all((B(1) == -B(2:end))) && (abs(B(1)==1))
            SimplexInfo(end+1,getvariables(C)) = sparse(1);
        end
    end
end
binaryVariables = union(getvariables(F(find(is(F,'binary')))),yalmip('binvariables'));

[Fchance,eliminatedConstraints,recursive] = deriveChanceModel(groupedChanceConstraints,randomVariables,options,SimplexInfo,binaryVariables);
Fchance = Fchance + F(find(keep)) + F(find(keep(~eliminatedConstraints)));
if recursive
    Fchance = modelchance(Fchance,options,1);
end
if ~rec && options.verbose
    disp('***** Modeling of chance constraints done. ********************************')
end


function [Fchance,eliminatedConstraints,recursive] = deriveChanceModel(groupedChanceConstraints,randomVariables,options,SimplexInfo,binaryVariables)

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
    
    for ic = 1:length(groupedChanceConstraints)
        if length(groupedChanceConstraints{ic})>1
            error('Joint chance constraint not supported');
        end
        if ~is(groupedChanceConstraints{ic},'elementwise')
            error('Only elementwise chance constraints supported')
        end
        
        confidencelevel = struct(groupedChanceConstraints{ic}).clauses{1}.confidencelevel;
        gamma = 1-confidencelevel;
        
        X = sdpvar(groupedChanceConstraints{ic});
                
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
            
            % Derive model b(x) + c(x)'*w >= 0
            fail = 0;
            b = [];
            c = [];
            for i = 1:length(X)
                [Qxwi,cxi,bi_,cwi,fail] = quadraticDecomposition(X(i),x,w);
                if isempty(bi_)
                    bi_ = 0;
                end
                bi = bi_ + fX;
                if ~isempty(cxi)
                    bi = bi + cxi'*x;
                end
                ci = cwi';
                if ~isempty(Qxwi)
                    ci = ci + Qxwi'*x;
                end
                c = [c;ci'];
                b = [b;bi];
            end
            c = c';
            
            % Used in characteristics stuff but there notation is h(x)+g(x)^Tw <= 0
            % Assumes this is an individual constraint
            funcs.h = @(x)(-bi_-cxi(:)'*x);
            funcs.dh =@(x)(-cxi(:));
            funcs.g = @(x)(-cwi(:)-Qxwi'*x);
            funcs.dg =@(x)(-Qxwi');
            
            
            % Are all variables in c constrained to a simple (this
            % is exploited in the normalChancefilter
            DisjointWeight = 0;
            if ~isempty(SimplexInfo) && isa(c,'sdpvar')
                x_in_c = getvariables(c);
                for i = 1:size(SimplexInfo,1)
                    if isequal(find(SimplexInfo(i,:)),x_in_c)
                        if all(ismember(x_in_c,binaryVariables))
                            DisjointWeight = 1;
                        else
                            DisjointWeight = -1;
                        end
                    end
                end
            end
            
            newConstraint = [];
            if ~fail
                if strcmp(func2str(randomVariables{uncertaintyGroup}.distribution.generator),'random')
                    distName = randomVariables{uncertaintyGroup}.distribution.parameters{1};
                    switch distName
                        case 'dro'
                            newConstraint = droChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                            printout(options.verbose,'dro',randomVariables{uncertaintyGroup}.distribution,ic,length(groupedChanceConstraints));
                            eliminatedConstraints(ic)=1;
                        case 'moment'
                            if isequal(options.chance.method,'momentchebyshev')
                                newConstraint = momentChebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                            else
                                newConstraint = momentChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                            end
                            printout(options.verbose,'moment',randomVariables{uncertaintyGroup}.distribution);
                            eliminatedConstraints(ic)=1;
                        case 'momentf'
                            newConstraint = momentfactorizedChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                            printout(options.verbose,'factorized moment',randomVariables{uncertaintyGroup}.distribution,ic,length(groupedChanceConstraints));
                            eliminatedConstraints(ic)=1;
                        case {'normal'}
                            newConstraint = normalChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options,funcs,x,DisjointWeight);
                            printout(options.verbose,'exact normal',randomVariables{uncertaintyGroup}.distribution,ic,length(groupedChanceConstraints));
                            eliminatedConstraints(ic)=1;
                        case {'logistic', 'laplace','uniform','t','tlocationScale','cauchy'}
                            newConstraint = symmetricUnivariateChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options,funcs,x);
                            printout(options.verbose,['exact symmetric univariate'],randomVariables{uncertaintyGroup}.distribution,ic,length(groupedChanceConstraints));
                            eliminatedConstraints(ic)=1;
                        case {'stable'}
                            newConstraint = conditionallysymmetricUnivariateChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options,funcs,x);
                            printout(options.verbose,['exact conditionally symmetric univariate'],randomVariables{uncertaintyGroup}.distribution,ic,length(groupedChanceConstraints));
                            eliminatedConstraints(ic)=1;
                        case {'gamma','weibull','exponential'}
                            newConstraint = nonsymmetricUnivariateChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options,funcs,x);
                            printout(options.verbose,['exact nonsymmetric univariate'],randomVariables{uncertaintyGroup}.distribution,ic,length(groupedChanceConstraints));
                            eliminatedConstraints(ic)=1;
                        otherwise
                            switch options.chance.method
                                case 'dro'
                                    newConstraint = droChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                                case {'chebyshev','chebychev'}
                                    newConstraint = sampledchebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                                case {'momentchebyshev','momentchebychev'}
                                    newConstraint = sampledmomentChebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                                case {'moment'}
                                    newConstraint = sampledmomentChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                                case 'markov'
                                    newConstraint =  sampledmarkovChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                                case 'chernoff'
                                    newConstraint =  sampledchernoffChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                                case 'integer'
                                    newConstraint =  sampledMIChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                                case 'scenario'
                                    newConstraint =  sampledScenarioChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,w,options);
                                otherwise
                                    error('Chance modeling approach not recognized');
                            end
                            printout(options.verbose,options.chance.method,randomVariables{uncertaintyGroup}.distribution,ic,length(groupedChanceConstraints));
                            eliminatedConstraints(ic)=1;
                    end
                else
                    switch options.chance.method
                        case 'chebyshev'
                            newConstraint = sampledchebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                        case 'moment'
                            newConstraint = sampledmomentChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                        case 'momentchebyshev'
                            newConstraint = sampledmomentChebyshevChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                        case 'markov'
                            newConstraint =  sampledmarkovChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                        case 'chernoff'
                            newConstraint =  sampledchernoffChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
                        case 'integer'
                            newConstraint =  sampledMIChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,gamma,options);
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
end
if any(eliminatedConstraints == 0)
    for weirdconstraints = find(eliminatedConstraints(:)' == 0)
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
    if isempty(x) || (nnz(Q_xx)==0)
        b = [b;f];
    else
        b = [b;f+x'*Q_xx*x];
    end
    c_wTbase = [c_wTbase;c_w'];
end

function printout(verbose,method,distribution,count,amountSame)

if verbose && count == 1
    if strcmpi(func2str(distribution.generator),'random')
        disp([' - Using ''' method '''-filter on ' num2str(amountSame) ' constraint' pluralS(amountSame) ' with ''' distribution.parameters{1} ''' distribution']);
    else
        disp([' - Using ''' method '''-filter on ' num2str(amountsame) ' constraint' pluralS(amountSame) ' with data created by @' distribution.name']);
    end
end

function s = pluralS(n)
if n == 1;s = '';else s = 's';end


function randomVariables = setupCharacterstics(randomVariables)

for k = 1:length(randomVariables)
    if strcmpi(randomVariables{k}.distribution.type,'stochastic')
        if isequal(randomVariables{k}.distribution.generator,@random)
            switch randomVariables{k}.distribution.parameters{1}
                
                case 'normal'
                    % Warning. When we normalize all Gaussian versions, we
                    % finally use covariance matrix as second parameter
                    % (defining  'normal' uses vector std. dev. while
                    % 'mvnrnd' uses covariance. However, the characterstic
                    % is setup with std.dev as parameter
                    
                    phi = @(t,mu,sigma) exp(1i*t.*mu(:) - 0.5*sigma(:).^2.*t.^2);
                    dphi = @(t,mu,sigma) (1i*mu(:) - sigma(:).^2.*t).*exp(1i*t.*mu(:) - 0.5*sigma(:).^2.*t.^2);
                    reldphi = @(t,mu,sigma) dphi(t,mu,sigma)./phi(t,mu,sigma);
                    
                    randomVariables{k}.distribution.characteristicfunction = phi;
                    randomVariables{k}.distribution.characteristicfunction_derivative = dphi;
                    randomVariables{k}.distribution.characteristicfunction_relativederivative = reldphi;
                    
                case 'exponential'
                    phi = @(t,mu)1./(1-t*1i.*mu(:));
                    dphi = @(t,mu)1i*(mu(:))./(1-1i*t.*mu(:)).^2;
                    reldphi = @(t,mu) dphi(t,mu)./phi(t,mu);
                    
                    randomVariables{k}.distribution.characteristicfunction = phi;
                    randomVariables{k}.distribution.characteristicfunction_derivative = dphi;
                    randomVariables{k}.distribution.characteristicfunction_relativederivative = reldphi;
                    
                case 'logistic'
                    phi = @(t,mu,s)(exp(1i*mu(:).*t))./guarded_sinhc(pi*s(:).*t);
                    dphi = @(t,mu,s) guarded_logistic_derivative(t,mu,s);
                    reldphi = @(t,mu,s) dphi(t,mu,s)./phi(t,mu,s);
                    
                    randomVariables{k}.distribution.characteristicfunction = phi;
                    randomVariables{k}.distribution.characteristicfunction_derivative = dphi;
                    randomVariables{k}.distribution.characteristicfunction_relativederivative = reldphi;
                    
                case 'uniform'
                    phi = @(t,a,b)(guarded_expdiv(b(:).*t,a(:).*t,t.*(b(:)-a(:))));
                    dphi = @(t,a,b) (1 ./ (1i*(b(:)-a(:)) .* t.^2)) .* (t.*(1i*b(:).*exp(1i * b(:) .* t) - 1i*a(:) .* exp(1i * a(:) .* t)) - (exp(1i * b(:) .* t) - exp(1i * a(:) .* t)));
                    reldphi = @(t,a,b) dphi(t,a,b)./phi(t,a,b);
                    
                    randomVariables{k}.distribution.characteristicfunction = phi;
                    randomVariables{k}.distribution.characteristicfunction_derivative = dphi;
                    randomVariables{k}.distribution.characteristicfunction_relativederivative = reldphi;
                    
                case 'cauchy'
                    phi = @(t,mu,theta)(exp(1i*mu(:).*t - abs(t).*theta(:)));
                    dphi = @(t,mu,theta) ((1i*mu(:)-sign(t).*theta(:)).*exp(1i*mu(:).*t - abs(t).*theta(:)))
                    reldphi = @(t,mu,theta) dphi(t,mu,theta)./phi(t,mu,theta);
                    
                    randomVariables{k}.distribution.characteristicfunction = phi;
                    randomVariables{k}.distribution.characteristicfunction_derivative = dphi;
                    randomVariables{k}.distribution.characteristicfunction_relativederivative = reldphi;
                otherwise
            end
        end
    end
end

function y = guarded_sinhc(x)
y = sinh(x)./x;y(x==0)=1;

function y = guarded_expdiv(z1,z2,z3)
y = (exp(1i*z1)-exp(1i*z2))./(1i*z3);
y(z3==0) = 1;

function y = guarded_logistic_derivative(t,mu,s)

y = exp(1i*mu(:).*t).*(1i*mu(:).*(pi.*s(:).*t./(sinh(pi*s(:).*t))) + (pi*s(:).*(1./sinh(pi*s(:).*t))-pi^2*s(:).^2.*t./((tanh(pi*s(:).*t).*sinh(pi*s(:).*t)))));
if numel(t)==1 && t == 0
    % Simple case phi(t) for scalar t, so all elements to limit
    y = 1i*mu(:);
else
    % t is a matrix (i.e. a vectorized call computing several points at
    % once)
    % Two cases, a row vector t meaning vector phi(t)
    if size(t,1) == 1
        i = find(t == 0);
        y(:,i) = repmat(1i*mu,1,length(i));
    else
        [i,j] = find(t == 0);
        y(sub2ind(size(t),i,j)) = 1i*mu(i);
    end
end