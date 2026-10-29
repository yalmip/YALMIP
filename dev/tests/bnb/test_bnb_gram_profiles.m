function tests = test_bnb_gram_profiles
tests = functiontests(localfunctions);
end

function test_exhaustive_partition_oracle(testCase)
% Enumerate labelled partitions independently of the profile DP.
n=7;
labels=mod(floor(bsxfun(@rdivide,(0:3^n-1)',3.^(0:n-1))),3);
sizes=[sum(labels==0,2) sum(labels==1,2) sum(labels==2,2)];
valid=all(sizes==0 | (sizes>=2 & sizes<=3),2);
testCase.verifyEqual(nnz(valid),630);
for encoding=0:1
    [p,ii,jj]=model(n,3,1,2,encoding);
    gram=initGramTriangleSeparator(p,3);
    [cuts,profile]=deriveGramProfileCuts(p,gram);
    testCase.assertEqual(numel(profile),1);
    testCase.verifyEqual(profile.equalPairs,[5 5]);
    testCase.verifyEqual(size(cuts,1),2);
    bits=double(labels(valid,ii)==labels(valid,jj));
    if ~encoding, bits=1-bits; end
    slack=[ones(nnz(valid),1) bits]*cuts';
    testCase.verifyGreaterThanOrEqual(min(slack(:)),-1e-12);
end
end

function test_root_profiles_and_both_encodings(testCase)
for encoding=0:1
    [p,ii,jj]=model(24,5,4,5,encoding);
    x=4/23*ones(numel(p.c),1);
    if ~encoding, x=1-x; end
    [cuts,info]=separateGramProfiles(p,x);
    testCase.verifyEqual(size(cuts,1),25);
    testCase.verifyEqual(info.records.degree,[5 5]);
    testCase.verifyEqual(info.records.equalPairs,[60 60]);
    labels=ceil((1:24)/6);
    bits=double(labels(ii)==labels(jj));
    if ~encoding, bits=1-bits; end
    testCase.verifyGreaterThanOrEqual(min(cuts*[1;bits(:)]),-1e-12);
    testCase.verifyEmpty(separateGramProfiles(p,bits(:)));

    p=model(32,6,4,5,encoding);
    x=5/31*ones(numel(p.c),1);
    if ~encoding, x=1-x; end
    [cuts,info]=separateGramProfiles(p,x);
    testCase.verifyEqual(size(cuts,1),1);
    testCase.verifyEqual(info.records.equalPairs,[70 70]);
    testCase.verifyLessThan(full(cuts*[1;x]),-9.9);
end
end

function test_rejects_unproved_domains_and_rows(testCase)
p=model(24,5,4,5,1);
x=4/23*ones(numel(p.c),1);
q=p; q.binary_variables(1)=q.binary_variables(2);
testCase.verifyEmpty(separateGramProfiles(q,x));
q=p; q.integral_variables(1)=q.integral_variables(2);
testCase.verifyEmpty(separateGramProfiles(q,x));
q=p; q.lb(1)=-1;
testCase.verifyEmpty(separateGramProfiles(q,x));
q=p; q.options.bnb.cuts=0;
testCase.verifyEmpty(separateGramProfiles(q,x));
q=p; q.F_struc(1,1)=-3; % This vertex has a different degree lower bound.
testCase.verifyEmpty(separateGramProfiles(q,x));
q=p; column=find(q.F_struc(1,2:end),1)+1;
q.F_struc(1,column)=1+1e-12; % No tolerance-based integer-row inference.
testCase.verifyEmpty(separateGramProfiles(q,x));
q=p; q.F_struc(1,1)=NaN;
testCase.verifyEmpty(separateGramProfiles(q,x));
gram=initGramTriangleSeparator(p,24);
gram{1}.exactPattern=false;
testCase.verifyEmpty(deriveGramProfileCuts(p,gram));
gram=initGramTriangleSeparator(p,24);
gram{1}.pair(1,2)=gram{1}.pair(1,3);
gram{1}.pair(2,1)=gram{1}.pair(1,3);
testCase.verifyEmpty(deriveGramProfileCuts(p,gram));
end

function [p,ii,jj]=model(n,K,L,U,equalValue)
% An exact uniform Gram block, with one independent bit per pair.
p=emptyNumericalModel;
p.K.l=2*n; p.K.s=n; p.K.r=[]; p.K.c=[]; p.K.m=0;
p.K.sos=struct('type',[]); p.K.rank=[];
nv=n*(n-1)/2;
p.c=zeros(nv,1); p.f=0; p.Q=sparse(nv,nv);
p.isbinary=true(1,nv); p.binary_variables=1:nv;
p.integral_variables=1:nv; p.lb=zeros(nv,1); p.ub=ones(nv,1);
p.variabletype=zeros(nv,1); p.evalMap=[]; p.nonlinear=[];
p.semicont_variables=[]; p.options.bnb.cuts=1;
p.F_struc=sparse(2*n+n*n,nv+1);
pair=zeros(n); ii=zeros(nv,1); jj=ii; at=0;
other=-1/(K-1);
if equalValue, base=other; coefficient=1-other;
else, base=1; coefficient=other-1; end
for j=2:n
    for i=1:j-1
        at=at+1; pair(i,j)=at; pair(j,i)=at;
        ii(at)=i; jj(at)=j;
        rows=2*n+[i+(j-1)*n j+(i-1)*n];
        p.F_struc(rows,1)=base;
        p.F_struc(rows,1+at)=coefficient;
    end
end
p.F_struc(2*n+(1:n+1:n*n),1)=1;
for vertex=1:n
    edges=nonzeros(pair(vertex,:));
    if equalValue
        p.F_struc(vertex,1)=-L; p.F_struc(vertex,1+edges)=1;
        p.F_struc(n+vertex,1)=U; p.F_struc(n+vertex,1+edges)=-1;
    else
        p.F_struc(vertex,1)=n-1-L; p.F_struc(vertex,1+edges)=-1;
        p.F_struc(n+vertex,1)=U-(n-1); p.F_struc(n+vertex,1+edges)=1;
    end
end
end
