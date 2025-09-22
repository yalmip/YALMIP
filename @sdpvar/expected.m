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
        Ef = replace(Ef,randomVariables{i}.variables,expect_from_distr(randomVariables{i}.distribution));
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
return

function [c_w,c_x,Q_xx,Q_xw,Q_ww,f_] = quadraticDecomposition(X,x,w)
xw = [x;w];
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
    Q_ww{i} = Q(wind,wind);
    Q_xw{i} = Q(xind,wind);
    Q_xx{i} = Q(xind,xind);
    c_x{i} = c(xind);
    c_w{i} = c(wind); 
    f_{i} = f;
end

function [mu,S] = expect_from_distr(distribution)

switch func2str(distribution.generator)
    case 'random'
        switch distribution.parameters{1}
            case 'normal'
                mu = distribution.parameters{2};
                S = distribution.parameters{3};
                if min(size(S))==1 && max(size(mu))>1
                    S = diag(S);
                end
            case 'exponential'                           
                 mu = distribution.parameters{2};
                 S = diag(distribution.parameters{2}.^2);
            case 'logistic'
                 mu = distribution.parameters{2};
                 S = diag(distribution.parameters{3}.^2*pi^2/3);
            otherwise
        end
    otherwise
        error
end

