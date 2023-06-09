function prob = appendMosekSDPdata(F_struc,K,prob)

prob.bardim = K.s;
prob.barc.subj = [];
prob.barc.subk = [];
prob.barc.subl = [];
prob.barc.val = [];
prob.bara.subi = [];
prob.bara.subj = [];
prob.bara.subk = [];
prob.bara.subl = [];
prob.bara.val = [];

C = F_struc(:,1);
A = -F_struc(:,2:end);

% -- Faster fix by Shahar
tops = [1 cumsum(K.s.^2)+1];
top = 1+K.f+K.l+sum(K.q)+3*K.e+sum(K.p);
[ii,jj,kk] = find(A(top:top + sum(K.s.^2)-1,:));
allcon = floor(interp1(tops,1:length(tops),ii,'linear'));
all_iilocal = ii-tops(allcon)'+1;
a = all_iilocal;
b = K.s(allcon);
allcol = ceil(a(:)./b(:))';
allrow = a(:)' - (allcol-1).*b(:)';
allvar = jj;
allval = kk;
% sort (for backward compatibility?)
[~,ind_sort] = sort(allcon);
allcon = allcon(ind_sort)';
allcol = allcol(ind_sort)';
allrow = allrow(ind_sort)';
allvar = allvar(ind_sort)';
allval = allval(ind_sort)';
% --

keep = find(allrow >= allcol);
allcol = allcol(keep);
allrow = allrow(keep);
allcon = allcon(keep);
allvar = allvar(keep);
allval = allval(keep);
prob.bara.subi = [prob.bara.subi allvar];
prob.bara.subj = [prob.bara.subj allcon];
prob.bara.subk = [prob.bara.subk allrow];
prob.bara.subl = [prob.bara.subl allcol];
prob.bara.val = [prob.bara.val allval];

for j = 1:length(K.s)
    n = K.s(j);
    Ci = C(top:top+n^2-1);
    Ci = tril(reshape(Ci,n,n));
    [k,l,val] = find(Ci);
    prob.barc.subj = [prob.barc.subj j*ones(1,length(k))];
    prob.barc.subk = [prob.barc.subk k(:)'];
    prob.barc.subl = [prob.barc.subl l(:)'];
    prob.barc.val = [prob.barc.val val(:)'];     
    top = top + n^2;
end