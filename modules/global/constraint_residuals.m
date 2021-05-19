function [res,X,detailed] = constraint_residuals(p,x)
res= [];
X = [];
detailed.f = [];
detailed.l = [];
detailed.q = [];
detailed.e = [];
detailed.p = [];
detailed.s = [];
detailed.b = [];
if ~isempty(p.F_struc)
    r = p.F_struc*[1;x];
end
if any(p.K.f)
    detailed.f = -abs(r(1:p.K.f));
    res = detailed.f;
end
if any(p.K.l)
    top = startofLPCone(p.K);
    detailed.l = r(top:top+p.K.l-1);
    res = [res;detailed.l];
end
if any(p.K.q)
    top = startofSOCPCone(p.K);
    for i = 1:length(p.K.q)
        zi = r(top:top + p.K.q(i)-1);top = top + p.K.q(i);
        detailed.q = [detailed.q;zi(1)^2-zi(2:end)'*zi(2:end)];
    end
    res = [res;detailed.q];
end
if any(p.K.e)
    top = startofEXPCone(p.K);
    for i = p.K.e
        zi = r(top:top + 2);top = top + 3;
        if zi(2)==0
            detailed.e = [detailed.e;zi(3)^2];
        else
            detailed.e = [detailed.e;zi(3)^2-zi(2)*exp(zi(1)/zi(2))];
        end
    end
    res = [res;detailed.e];
end
if any(p.K.p)
    top = startofPOWCone(p.K);
    for i = length(p.K.p)
        zi = r(top:top + p.K.p(i)-1);top = top + p.K.p(i);
        a = zi(end);        
        detailed.p = [detailed.p;zi(1)^a*zi(2)^(1-a) - norm(zi(3:end-1))];
    end
    res = [res;detailed.p];
end
if any(p.K.s)
    top = startofSDPCone(p.K);
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X{i} = r(top:top+n^2-1);top = top+n^2;
        X{i} = reshape(X{i},n,n);
        detailed.s = [detailed.s;min(eig(X{i}))];
    end
    res = [res;detailed.s];
end
detailed.b = min([p.ub-x;x-p.lb]);
res = [res;detailed.b];