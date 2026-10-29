function [bestX,bestCost,info] = bnbGramPartitionRoot(p,budget,seed)
%BNBGRAMPARTITIONROOT Search in equivalence classes, not pair bits.
% Every returned point is checked in the full
% numerical model; the Gram structure only guides candidate generation.
if nargin<2, budget=0.2; end
if nargin<3, seed=739; end
bestX=[];
bestCost=inf;
info=struct('attempts',0,'feasible',0,'moves',0,'classes',0, ...
    'seed_checks',0,'seconds',0, ...
    'reason','unsupported');
if ~isscalar(budget) || ~isfinite(budget) || budget<=0, return; end
budget=min(budget,0.2);
t=tic;
if nnz(p.Q) || ~isempty(p.evalMap) || any(p.variabletype) || ...
        ~isempty(p.nonlinear) || ~isempty(p.semicont_variables) || ...
        numel(p.K.s)~=1 || p.K.f || any(p.K.c) || ...
        any(p.K.q) || any(p.K.r) || any(p.K.e) || any(p.K.p) || ...
        ~isequal(p.K.m,0) || ~isempty(p.K.sos.type) || ...
        (~isempty(p.K.rank) && ...
         (numel(p.K.rank)~=1 || any(p.K.rank(:)<p.K.s(1)))) || ...
        nnz(p.F_struc)>2e5 || size(p.F_struc,1)>5e4 || ...
        ~isfinite(p.f) || any(~isfinite(p.c)) || ...
        any(~isfinite(nonzeros(p.F_struc))) || ...
        ~isfinite(p.options.bnb.feastol) || p.options.bnb.feastol<=0
    info.reason='model/cone gate';
    return
end
% Uniform pair costs make a short swap search effective and keep it away
% from larger weighted partition instances where an early weak incumbent
% changed the BNB tree for the worse.
activeCosts=abs(p.c(abs(p.c)>1e-10));
if p.K.s(1)<24 || p.K.s(1)>30 || isempty(activeCosts) || ...
        max(activeCosts)-min(activeCosts)>1e-8*max(1,max(activeCosts))
    info.reason='bounded uniform-cost gate';
    return
end
gram=initGramTriangleSeparator(p,24);
if numel(gram)~=1
    info.reason='gram detector';
    return
end
g=gram{1};
pair=g.pair;
n=size(pair,1);
if n<24 || n>30 || nnz(triu(pair,1))~=n*(n-1)/2
    info.reason='pair coverage';
    return
end
vars=nonzeros(triu(pair,1));
if numel(unique(vars))~=numel(vars) || numel(vars)~=numel(p.c) || ...
        any(~ismember(vars,p.binary_variables)) || ...
        any(~ismember(vars,p.integral_variables)) || ...
        any(p.lb(vars)~=0) || any(p.ub(vars)~=1) || ...
        ~isfinite(g.otherEntry) || g.otherEntry>=0 || g.diagonal<=0
    info.reason='variable domain';
    return
end
% d+(r-1)b is the smallest eigenvalue for r distinct classes.
classLimit=1-g.diagonal/g.otherEntry;
if abs(classLimit-round(classLimit))<1e-6
    K=round(classLimit);
else
    K=floor(classLimit);
end
if K<2 || K>n || g.diagonal+K*g.otherEntry>=-1e-7 || ...
        g.diagonal+(K-1)*g.otherEntry<-1e-7
    info.reason=sprintf('class limit d=%.17g b=%.17g K=%d below=%.4g above=%.4g', ...
        g.diagonal,g.otherEntry,K,g.diagonal+(K-1)*g.otherEntry, ...
        g.diagonal+K*g.otherEntry);
    return
end
W=zeros(n);
for j=2:n
    for i=1:j-1
        v=pair(i,j);
        if g.equalValue==1
            W(i,j)=p.c(v);
        else
            W(i,j)=-p.c(v);
        end
        W(j,i)=W(i,j);
    end
end
oldRng=rng;
cleanup=onCleanup(@()rng(oldRng));
rng(seed,'twister');
feasTol=min(1e-8,p.options.bnb.feastol);
foundSeed=false;
for candidateK=K:-1:2
    sizes=floor(n/candidateK)*ones(candidateK,1);
    sizes(1:mod(n,candidateK))=sizes(1:mod(n,candidateK))+1;
    for seedTry=1:3
        if toc(t)>=budget, break; end
        labels=balancedGramLabels(sizes);
        labels=labels(randperm(n));
        x=decodeGramPartition(labels,pair,g.equalValue,numel(p.c));
        info.seed_checks=info.seed_checks+1;
        if checkfeasiblefast(p,x,feasTol)
            foundSeed=true;
            bestX=x;
            bestCost=computecost(p.f,p.c,p.Q,x,p);
            info.feasible=1;
            info.classes=candidateK;
            break
        end
    end
    if foundSeed || toc(t)>=budget, break; end
end
if ~foundSeed
    info.reason='no feasible class-count seed';
    info.seconds=toc(t);
    return
end
while toc(t)<budget && info.attempts<200
    info.attempts=info.attempts+1;
    labels=balancedGramLabels(sizes);
    labels=labels(randperm(n));
    for pass=1:100
        if toc(t)>=budget, break; end
        improvement=0;
        ai=0;
        bj=0;
        for i=1:n-1
            if toc(t)>=budget, break; end
            a=labels(i);
            ia=(labels==a);
            ia(i)=false;
            for j=i+1:n
                if labels(j)==a, continue; end
                b=labels(j);
                ib=(labels==b);
                ib(j)=false;
                delta=sum(W(j,ia)-W(i,ia))+sum(W(i,ib)-W(j,ib));
                if delta<improvement-1e-10
                    improvement=delta;
                    ai=i;
                    bj=j;
                end
            end
        end
        if ai==0, break; end
        tmp=labels(ai); labels(ai)=labels(bj); labels(bj)=tmp;
        info.moves=info.moves+1;
    end
    x=decodeGramPartition(labels,pair,g.equalValue,numel(p.c));
    cost=computecost(p.f,p.c,p.Q,x,p);
    if cost<bestCost && checkfeasiblefast(p,x,feasTol)
        bestX=x;
        bestCost=cost;
        info.feasible=info.feasible+1;
    end
end
info.seconds=toc(t);
if isempty(bestX), info.reason='no feasible candidate';
else, info.reason='checked incumbent'; end
end

function x=decodeGramPartition(labels,pair,equalValue,nvar)
x=zeros(nvar,1);
n=numel(labels);
for j=2:n
    for i=1:j-1
        x(pair(i,j))=double(labels(i)==labels(j))==equalValue;
    end
end
end

function labels=balancedGramLabels(sizes)
labels=zeros(sum(sizes),1);
last=0;
for k=1:numel(sizes)
    labels(last+1:last+sizes(k))=k;
    last=last+sizes(k);
end
end
