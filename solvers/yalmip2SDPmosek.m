function prob = yalmip2SDPmosek(model)

% Extract to sedumi format
model.C = model.F_struc(:,1);
model.A = -model.F_struc(:,2:end);
model.b = -model.c;

K = model.K;

prob.a = model.A(1:(K.f+K.l+sum(K.q)),:)';
prob.c = model.C(1:(K.f+K.l+sum(K.q)));

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

for j = 1:length(model.K.s)
    n = model.K.s(j);
    Ci = model.C(top:top+n^2-1);
    Ci = tril(reshape(Ci,n,n));
    [k,l,val] = find(Ci);
    prob.barc.subj = [prob.barc.subj j*ones(1,length(k))];
    prob.barc.subk = [prob.barc.subk k(:)'];
    prob.barc.subl = [prob.barc.subl l(:)'];
    prob.barc.val = [prob.barc.val val(:)'];
    
    Ais = model.A(top:top+n^2-1,:);
    [elementIndex,variableIndex,val] = find(Ais);
    [row,column] = ind2sub([n n],elementIndex);
    keep = find(row >= column);
    row = row(keep);
    column = column(keep);  
    variableIndex = variableIndex(keep);
    val = val(keep);
    
    % Which LMI
    prob.bara.subj = [prob.bara.subj j*ones(1,length(val))];
    % Which variable
    prob.bara.subi = [prob.bara.subi variableIndex(:)'];
    % Row
    prob.bara.subk = [prob.bara.subk row(:)'];
    % Column
    prob.bara.subl = [prob.bara.subl column(:)'];
    % Value
    prob.bara.val = [prob.bara.val val(:)'];
    
    top = top + n^2;
end

if K.q(1)>0    
    prob.cones.type   = [repmat(0,1,length(K.q))];
    top = 1 + K.f + K.l;
    prob.cones.sub = [];
    prob.cones.subptr = [];
    for i = 1:length(K.q)
        prob.cones.subptr = [ prob.cones.subptr 1+length(prob.cones.sub)];
        prob.cones.sub = [prob.cones.sub top:top+K.q(i)-1];
        top = top + K.q(i);
    end
end
