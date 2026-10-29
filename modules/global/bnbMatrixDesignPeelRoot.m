function [xBest,bestCost] = bnbMatrixDesignPeelRoot(p,xBest,bestCost)
%BNBMATRIXDESIGNPEELROOT Find a cheaper feasible ordered-level SDP design.
% This root heuristic only changes an incumbent after complete model checks.
% The caller checks the cheap group-count gate before entering this function.
groups=p.oneHotLevelGroups;
if ~iscell(groups) || numel(groups)<16 || numel(groups)>64 || ...
        numel(p.K.s)~=1 || p.K.f || any(p.K.c) || ...
        any(p.K.q) || any(p.K.r) || ...
        any(p.K.e) || any(p.K.p) || ~isequal(p.K.m,0) || ...
        ~isempty(p.K.sos.type) || ...
        nnz(p.Q) || any(p.variabletype) || ~isempty(p.evalMap) || ...
        ~isempty(p.nonlinear) || ~isempty(p.semicont_variables)
    return
end
n=p.K.s(1);
if n<1 || n>48 || n~=round(n), return; end
if ~isempty(p.K.rank) && ...
        (numel(p.K.rank)~=1 || any(p.K.rank(:)<n))
    return
end
timer=tic;
if sum(cellfun(@numel,groups))>256 || ...
        ~isfinite(p.f) || any(~isfinite(p.c)) || ...
        nnz(p.F_struc)>2e5 || ...
        any(~isfinite(nonzeros(p.F_struc)))
    return
end
nvar=numel(p.c);
members=cell2mat(cellfun(@(g)g(:),groups(:),'UniformOutput',false));
if isempty(members) || numel(unique(members))~=numel(members) || ...
        any(~isfinite(members)) || any(members~=round(members)) || ...
        any(members<1) || any(members>nvar) || ...
        ~all(ismember(members,p.binary_variables)) || ...
        ~all(ismember(members,p.integral_variables)) || ...
        any(~isfinite(p.lb(members))) || any(~isfinite(p.ub(members))) || ...
        any(p.lb(members)<-1e-9) || any(p.ub(members)>1+1e-9) || ...
        any(p.c(members)<-1e-9)
    return
end
top=startofSDPCone(p.K);
rows=top:top+n*n-1;
% Only PSD coefficient columns have a predictable spectral downdate.
for j=members(:)'
    if toc(timer)>0.5, return; end
    A=reshape(full(p.F_struc(rows,1+j)),n,n);
    A=(A+A')/2;
    if any(~isfinite(A(:))), return; end
    if min(eig(A)) < -1e-9*max(1,norm(A,'fro'))
        return
    end
end
x=zeros(nvar,1);
fixed=p.lb==p.ub;
x(fixed)=p.lb(fixed);
for k=1:numel(groups)
    group=groups{k}(:);
    x(group)=min(1,p.ub(group));
end
x=max(x,p.lb);
feasTol=min(1e-8,p.options.bnb.feastol);
if toc(timer)>0.5 || ~isfinite(feasTol) || feasTol<=0 || ...
        any(~isfinite(x)) || ...
        any(abs(x(p.integral_variables)-round(x(p.integral_variables)))> ...
            p.options.bnb.inttol) || ...
        ~checkfeasiblefast(p,x,feasTol)
    return
end
M=reshape(full(p.F_struc(rows,:)*[1;x]),n,n);
M=(M+M')/2;
if any(~isfinite(M(:))), return; end
checks=0;
while checks<3000 && toc(timer)<=0.5
    choice=0;
    bestScore=-inf;
    currentMin=min(eig(M));
    if ~isfinite(currentMin), return; end
    for k=1:numel(groups)
        if checks>=3000 || toc(timer)>0.5, break; end
        group=groups{k}(:);
        active=find(x(group)>0.5,1,'last');
        if isempty(active), continue; end
        member=group(active);
        if p.lb(member)>0.5 || p.c(member)<=1e-10, continue; end
        trial=x;
        trial(member)=0;
        checks=checks+1;
        if ~checkfeasiblefast(p,trial,feasTol), continue; end
        A=reshape(full(p.F_struc(rows,1+member)),n,n);
        T=M-(A+A')/2;
        if any(~isfinite(T(:))), continue; end
        spectralLoss=max(0,currentMin-min(eig(T)));
        score=p.c(member)/(1e-5+spectralLoss);
        if score>bestScore
            bestScore=score;
            choice=member;
        end
    end
    if ~choice, break; end
    x(choice)=0;
    M=reshape(full(p.F_struc(rows,:)*[1;x]),n,n);
    M=(M+M')/2;
    if any(~isfinite(M(:))), return; end
end
cost=computecost(p.f,p.c,p.Q,x,p);
if isfinite(cost) && cost<bestCost-1e-8 && ...
        checkfeasiblefast(p,x,feasTol)
    xBest=x;
    bestCost=cost;
end
end
