function p = mergeNumericalModels(p0,p1)

% It is assumed p0 is the master model, i.e. p1 only contains constraints
p = p0;

if p1.K.f > 0
    p.F_struc = [p1.F_struc(1:K.f,:);p.F_struc];
    p.K.f = p.K.f + p1.K.f;
end

if p1.K.l > 0
    p.F_struc = [p.F_struc(1:p.K.f,:);
        p1.F_struc(1+p1.K.f:p1.K.f+p1.K.l,:)
        p.F_struc(1+p.K.f:end,:)];
    p.K.l = p.K.l + p1.K.l;
end

% if p1.K.q > 0
% end

% if p1.K.s > 0
%     p.F_struc = [p.F_struc;p1.F_struc(1+p1.K.f+p1.K.l+sum(p1.K.q):end,:)];
%     p.K.s = [p.K.s p1.K.s];p.K.s(p1.K.s==0)=[];
% end

