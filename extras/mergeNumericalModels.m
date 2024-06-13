function p = mergeNumericalModels(p0,p1)

if ~isempty(p1.F_struc)
    if size(p1.F_struc,2) < size(p0.F_struc,2)
        p1.F_struc(end,size(p0.F_struc,2))=0;
    end
else
    p = p0;
    return;
end
if isempty(p0.F_struc) & isempty(p1.F_struc)
    p = p1;
    return
end

if (size(p0.F_struc,2) < size(p1.F_struc,2)) & (size(p0.F_struc,1)>0)
    p0.F_struc(end,size(p1.F_struc,2))=0;
end

p1.F_struc = sparse(p1.F_struc);
p0.F_struc = sparse(p0.F_struc);

% It is assumed p0 is the master model, i.e. p1 only contains constraints
p = p0;

if any(p1.K.f)
    p.F_struc = [p1.F_struc(1:p1.K.f,:);p.F_struc];
    p.K.f = p.K.f + p1.K.f;
end

if any(p1.K.l)
    p.F_struc = [p.F_struc(1:p.K.f,:);
        p1.F_struc(1+p1.K.f:p1.K.f+p1.K.l,:)
        p.F_struc(1+p.K.f:end,:)];
    p.K.l = p.K.l + p1.K.l;
end

if isfield(p1.K,'q') & nnz(p1.K.q)>0
      p.F_struc = [p.F_struc(1:p.K.f+p.K.l,:);
                   p1.F_struc(1+p1.K.f+p1.K.l:p1.K.f+p1.K.l+sum(p1.K.q),:)
                   p.F_struc(1+p.K.f+p.K.l:end,:)];
    p.K.q = [p1.K.q p.K.q];
    p.K.q = p.K.q(find(p.K.q));
end

if isfield(p1.K,'e') & nnz(p1.K.e)>0
      p.F_struc = [p.F_struc(1:p.K.f+p.K.l+sum(p.K.q),:);
                   p1.F_struc(1+p1.K.f+p1.K.l+sum(p1.K.q):p1.K.f+p1.K.l+sum(p1.K.q)+3*p1.K.e,:)
                   p.F_struc(1+p.K.f+p.K.l+sum(p.K.q):end,:)];
   p.K.e = p.K.e + p1.K.e;
end

if isfield(p1.K,'s')
    if p1.K.s > 0
        error('FIXME')
        %     p.F_struc = [p.F_struc;p1.F_struc(1+p1.K.f+p1.K.l+sum(p1.K.q):end,:)];
        %     p.K.s = [p.K.s p1.K.s];p.K.s(p1.K.s==0)=[];
    end
end

