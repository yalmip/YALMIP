function prob = yalmip2SDPmosek(model)

% Extract to sedumi format
%model.C = model.F_struc(:,1);
%model.A = -model.F_struc(:,2:end);
model.b = -model.c;

K = model.K;

%prob.a = model.A(1:(K.f+K.l+sum(K.q)),:)';
%prob.c = model.C(1:(K.f+K.l+sum(K.q)));
prob.a = -model.F_struc(1:(K.f+K.l+sum(K.q)),2:end)';
prob.c = model.F_struc(1:(K.f+K.l+sum(K.q)),1);


prob.blx = [-inf(K.f,1);zeros(K.l,1);-inf(sum(K.q),1)];
prob.bux = [inf(K.f+K.l+sum(K.q),1)];

prob.bardim = model.K.s;
prob.blc = model.b;
prob.buc = model.b;

top = 1+K.f+K.l+sum(K.q);
prob.barc.subj = [];
prob.barc.subk = [];
prob.barc.subl = [];
prob.barc.val = [];
prob.bara.subi = [];
prob.bara.subj = [];
prob.bara.subk = [];
prob.bara.subl = [];
prob.bara.val = [];

tops = [1];
for j = 1:length(model.K.s)
    n = model.K.s(j);
    tops = [tops tops(end)+n^2];
end
[ii,jj,kk] = find(model.F_struc(top:top + sum(K.s.^2)-1,2:end));
%[ii,jj,kk] = find(model.A(top:top + sum(K.s.^2)-1,:));
cols = zeros(length(ii),1);
rows = zeros(length(ii),1);
allcol = [];
allrow = [];
allcon = [];
allvar = [];
allval = [];
for j = 1:length(model.K.s)    
    ind = find(ii>=tops(j) & ii<=tops(j+1)-1);
    iilocal = ii(ind)-tops(j)+1;
    col = ceil(iilocal/model.K.s(j));
    row = iilocal - (col-1)*model.K.s(j);
    allcol = [allcol col(:)'];
    allrow = [allrow row(:)'];
    allvar = [allvar jj(ind(:))'];
    allval = [allval -kk(ind(:))'];
    %allval = [allval kk(ind(:))'];
    allcon = [allcon repmat(j,1,length(col))];
end
keep = find(allrow >= allcol);
allcol = allcol(keep);
allrow = allrow(keep);
allcon = allcon(keep);
allvar = allvar(keep);
allval = allval(keep);
%allvar = jj(keep);
%allval = kk(keep);
prob.bara.subi = [prob.bara.subi allvar];
prob.bara.subj = [prob.bara.subj allcon];
prob.bara.subk = [prob.bara.subk allrow];
prob.bara.subl = [prob.bara.subl allcol];
prob.bara.val = [prob.bara.val allval];

for j = 1:length(model.K.s)
    n = model.K.s(j);
    %Ci = model.C(top:top+n^2-1);
    Ci = model.F_struc(top:top+n^2-1,1);
    Ci = tril(reshape(Ci,n,n));
    [k,l,val] = find(Ci);
    prob.barc.subj = [prob.barc.subj j*ones(1,length(k))];
    prob.barc.subk = [prob.barc.subk k(:)'];
    prob.barc.subl = [prob.barc.subl l(:)'];
    prob.barc.val = [prob.barc.val val(:)'];     
    top = top + n^2;
end

if any(K.q)    
    prob.cones.type   = [repmat(0,1,length(K.q))];
    top = startofSOCPCone(K);
    prob.cones.sub = [];
    prob.cones.subptr = [];
    for i = 1:length(K.q)
        prob.cones.subptr = [ prob.cones.subptr 1+length(prob.cones.sub)];
        prob.cones.sub = [prob.cones.sub top:top+K.q(i)-1];
        top = top + K.q(i);
    end
end
