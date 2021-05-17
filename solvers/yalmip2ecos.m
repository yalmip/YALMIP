function ecosmodel = yalmip2ecos(model)

% Retrieve needed data
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
ub      = model.ub;
lb      = model.lb;

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

if any(K.f)
    b = -full(F_struc(1:K.f,1));
    A = F_struc(1:K.f,2:end);
    h = full(F_struc(1+K.f:end,1));        
    G = -F_struc(1+K.f:end,2:end);
else
    A = [];
    b = [];
    h = full(F_struc(:,1));
    G = -F_struc;
    G(:,1)=[];  
end

dims.l = K.l;
if ~any(K.q)
    dims.q = [];
else
    dims.q = K.q;
end

ecosmodel.A = A;
ecosmodel.b = b;
ecosmodel.G = G;
ecosmodel.h = h;
ecosmodel.c = c;
ecosmodel.dims = dims;
