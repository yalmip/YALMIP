function Ef = expected(f)
    
randomVariables = yalmip('getDistribution');
[randomVariables,map] = mergeDistributions(randomVariables);
allwVars = [];

for i = 1:length(randomVariables)
    allwVars = [allwVars;getvariables(randomVariables{i}.variables)];
end

if all(degree(f,recover(allwVars)) <= 1)
    Ef = f;
    for i = 1:length(randomVariables)
        Ef = replace(Ef,randomVariables{i}.variables,expect_from_distr(randomVariables{i}.distribution));
    end
    return
elseif length(randomVariables) == 1 && isequal(randomVariables{1}.distribution.parameters{1},'normal')

    w = recover(allwVars);
    x = recover(setdiff(depends(f),allwVars));
    [c_w,c_x,Q_xx,Q_xw,Q_ww,f_] = quadraticDecomposition(f,x,w);         
    [Ew,Eww] = expect_from_distr(randomVariables{1}.distribution);  
    Ef = [];
    for i = 1:length(c_w)
        Ef = [Ef;x'*Q_xx{i}*x + c_x{i}'*x + c_w{i}'*Ew + 2*x'*Q_xw{i}*Ew+f_{i} + trace(Q_ww{i}*Eww)];
    end
    Ef = reshape(Ef,size(f));
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

mu = distribution.parameters{2};
S = distribution.parameters{3};
if min(size(S))==1 && max(size(mu))>1
    S = diag(S);
end


