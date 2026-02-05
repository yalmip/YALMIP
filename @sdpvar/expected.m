function Ef = expected(f)
    
% Gather all distribution definitions, and prune them for the ones used in
% this expressions
randomVariables = yalmip('getDistribution');

keep = zeros(1,length(randomVariables));
f_variables = depends(f);
for i = 1:length(randomVariables)
    if any(ismember(f_variables,getvariables(randomVariables{i}.variables)))
        keep(i) = 1;
    end
end
randomVariables = {randomVariables{find(keep)}};

% Scalar definitions of same distribution bunched together
% This merges various normal/gaussian models, beware
[randomVariables,map] = mergeDistributions(randomVariables);

for i = 1:length(randomVariables)
    randomVariables{i}.distribution = combineElementMixtures(randomVariables{i}.distribution,randomVariables{i}.id);
end

% Get a variable representing them
allwVars = [];
for i = 1:length(randomVariables)
    allwVars = [allwVars getvariables(randomVariables{i}.variables)];
end
w = recover(allwVars(:));

if all(degree(f,w) <= 1)
    % Simple case, expression is linear in random variable
    Ef = f;
    for i = 1:length(randomVariables)
        Ef = replace(Ef,randomVariables{i}.variables,expect_from_distr(randomVariables{i}.distribution,randomVariables{i}.id));
    end
    return
elseif length(randomVariables) == 1 && degree(f) <= 2
    x = recover(setdiff(depends(f),allwVars));
    [c_w,c_x,Q_xx,Q_xw,Q_ww,f_] = quadraticDecomposition(f,x,w);         
    [Ew,S] = expect_from_distr(randomVariables{1}.distribution);  
    Eww = S + Ew*Ew';
    Ef = [];
    if isempty(x)
        for i = 1:length(c_w)
            Ef = [Ef;c_w{i}'*Ew + f_{i} + trace(Q_ww{i}*Eww)];
        end
    else
    for i = 1:length(c_w)
        Ef = [Ef;x'*Q_xx{i}*x + c_x{i}'*x + c_w{i}'*Ew + 2*x'*Q_xw{i}*Ew+f_{i} + trace(Q_ww{i}*Eww)];
    end
    end
    Ef = reshape(Ef,size(f));
    
elseif degree(f) <= 2    
    
    for i = 1:length(randomVariables)
        w = randomVariables{i}.variables
        x = recover(setdiff(depends(f),getvariables(w)));
        [c_w,c_x,Q_xx,Q_xw,Q_ww,f_] = quadraticDecomposition(f,x,w);
        [Ew,S] = expect_from_distr(randomVariables{1}.distribution);
        Eww = S + Ew*Ew';
        Ef = [];
        if isempty(x)
            for i = 1:length(c_w)
                Ef = [Ef;c_w{i}'*Ew + f_{i} + trace(Q_ww{i}*Eww)];
            end
        else
            for i = 1:length(c_w)
                Ef = [Ef;x'*Q_xx{i}*x + c_x{i}'*x + c_w{i}'*Ew + 2*x'*Q_xw{i}*Ew + f_{i} + trace(Q_ww{i}*Eww)];
            end
        end
        f = reshape(Ef,size(f));
    end       
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

function [mu,S] = expect_from_distr(distribution,id)

if length(unique(id)) > 1 && ~isempty(distribution.mixture)
    % We have a concatenation of mixtures, so we must generate the new
    % mixture by combinatorial combinations
    
end
    
switch func2str(distribution.generator)
    case 'random'
        parameters = {distribution.parameters{2:end}};
        if ~isempty(distribution.mixture)            
            mu = 0;
            S  = 0;
            for i = 1:size(distribution.mixture,2)
                parameters_i = parameter_components(parameters,i);
                [mu_i,S_i] = extract_mu_sigma(distribution.parameters{1},parameters_i);
                mu = mu + mu_i.*distribution.mixture(:,i);
                S = S + (S_i + mu_i*mu_i').*distribution.mixture(:,i);
            end            
            S = S-mu*mu';
        else            
            [mu,S] = extract_mu_sigma(distribution.parameters{1},parameters);
        end
                
    otherwise
        error('FIXME: expected all cases')
end

function parameters_i = parameter_components(parameters,i)
for k = 1:length(parameters)
    parameters_i{k} = parameters{k}{i};
end

function [mu,S] = extract_mu_sigma(name,parameters)
switch name
    case 'normal'
        % Careful. We redefine the normal, second parameters is (co)variance
        mu = parameters{1};
        S = parameters{2};
        if min(size(S))==1 && max(size(mu))>1
            S = diag(S);
        end
    case 'exponential'
        mu = parameters{1};
        S = diag(parameters{1}.^2);
    case 'logistic'
        mu = parameters{1};
        S = diag(parameters{2}.^2*pi^2/3);
    otherwise
        error('FIXME: expected all cases')
end
