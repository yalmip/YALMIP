function Ef = expected(f)
    
randomVariables = yalmip('getDistribution');
if length(randomVariables) == 0
    Ef = f;
    return
end

% Prune unused random declarations to avoid extra work
allVars = depends(f);
allwVars = [];
keep = false(1,length(randomVariables));
for i = 1:length(randomVariables)
    wVars = getvariables(randomVariables{i}.variables);
    if any(ismember(allVars, wVars))
        allwVars = [allwVars;wVars];    
        keep(i) = true;
    end
end
randomVariables = randomVariables(keep);

% Collect and normalize (put same distributions together, use covariance as
% second parameter in 'normal')
randomVariables = mergeDistributions(randomVariables);

% So we have decision variables x, and random variables w
w = recover(allwVars);
x = recover(setdiff(depends(f),allwVars));
% Code assumes vector so temporarily reshape
[n,m] = size(f);
f = reshape(f,[],1);
if linearin(f,w)
    % Expression is linear in w, we can simply replace w with expected
    % value and we do this using high-level code for now. 
    % FIXME: speed up simple cases
    Ef = f;
    for i = 1:length(randomVariables)
        Ef = replace(Ef,randomVariables{i}.variables,expect_from_distr(randomVariables{i}.distribution));
    end    
elseif degree(f) <= 2
    Ef = quadraticExpected(f,x,w,randomVariables);   
elseif degree(f,recover(allwVars)) <= 2
    % At least it is quadratic in w, so write as c(x)*p(w)    
    [c,p] = coefficients(f,w);
    % Note sure why coefficients adds higher degree stuff?
    [c,p] = pruneHigherThanQuadratic(c,p);
    Ep = quadraticExpected(p,[],w,randomVariables);
    Ef = c*Ep;
else
    error('Cannot compute this expectation')
end
Ef = reshape(Ef,n,m);

function [c_w,c_x,Q_xx,Q_xw,Q_ww,f_] = quadraticDecomposition(X,x,w)
xw = [x;w];
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
[Qs,cs,fs,dummy,h_] = vecquaddecomp(X,xw);
c_wTbase = [];
AAA = [];
ccc = [];
for i = 1:length(X)
    Q = Qs{i};
    c = cs{i};
    f = fs{i};    
    Q_ww{i} = Q(wind,wind);
    Q_xw{i} = Q(xind,wind);
    Q_xx{i} = Q(xind,xind);
    c_x{i} = c(xind);
    c_w{i} = c(wind); 
    f_{i} = f;
end


function [mu,Covariance] = expect_from_distr(distribution)
switch distribution.parameters{1}
    case {'normal','exponential'}
        mu = distribution.parameters{2};
    case 'data'    
        mu = mean(distribution.parameters{2},2);
    otherwise
end
if nargout == 2
    switch distribution.parameters{1}
        case 'normal'
            % Note, internal normalization so parameter is covariance!
            Covariance = distribution.parameters{3};  
        case 'exponential'
            Covariance = diag(distribution.parameters{2}.^2);  
        case 'data'
            Covariance = cov(distribution.parameters{2});              
        otherwise
    end
end

function Ef = quadraticExpected(f,x,w,randomVariables)
% Fully quadratic case can be done semi-efficiently
[c_w,c_x,Q_xx,Q_xw,Q_ww,f_] = quadraticDecomposition(f,x,w);
Ew = [];
Cov = [];
localwVars = [];
for i = 1:length(randomVariables)
    [Ew_,Cov_] = expect_from_distr(randomVariables{i}.distribution);
    localwVars = [localwVars;getvariables(randomVariables{i}.variables)];
    Ew = [Ew;Ew_];
    Cov = blkdiag(Cov,Cov_);
    % FIXME: map to ordering in w
end
Ef = [];
for i = 1:length(c_w)
    if isempty(x)
        Ef = [Ef;f_{i} + c_w{i}'*Ew + Ew'*Q_ww{i}*Ew + trace(Q_ww{i}*Cov)];
    else
        Ef = [Ef;x'*Q_xx{i}*x + c_x{i}'*x + c_w{i}'*Ew + 2*x'*Q_xw{i}*Ew+f_{i} + Ew'*Q_ww{i}*Ew + trace(Q_ww{i}*Cov)];
    end
end

