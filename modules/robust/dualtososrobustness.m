function SOSModel = dualtososrobustness(UncertainConstraint,UncertaintySet,UncertainVariables,DecisionVariables,p_tau_degree,localizer_tau_degree,Z_degree)
  
[E,F] = getEFfromSET(UncertaintySet);
[F0,Fz,Fx] = getFzxfromSET(UncertainConstraint,UncertainVariables,DecisionVariables);

if is(UncertainConstraint,'sdp')
    n = length(F0);
    v = sdpvar(n,1);
    d = v'*F0*v;
    b = [];for i = 1:length(Fx);b = [b;v'*Fx{i}*v];end
    c = [];for i = 1:length(Fz);c = [c;v'*Fz{i}*v];end
elseif is(UncertainConstraint,'socp')
    n = length(F0)-1;
    v = sdpvar(n,1);
    d = [1 v']*F0;
    b = [];for i = 1:length(Fx);b = [b;[1 v']*Fx{i}];end
    c = [];for i = 1:length(Fz);c = [c;[1 v']*Fz{i}];end
end
[Z,coeffs] = createDualParameterization(UncertaintySet,v,Z_degree);

Zblock = blkdiag(Z{:});
coeffs = [DecisionVariables;coeffs(:)]

D = [];
for i = 1:length(F)
    D = [D, coefficients(trace(Zblock'*F{i})-c(i),v)==0];
end

gv = (1-v'*v);
for i = 1:length(Z)       
    if is(UncertaintySet(i),'sdp')
        [tau,coefftau] = polynomial(v,localizer_tau_degree);
         coeffs = [coeffs;coefftau];
        D=[D,sos(Z{i}-eye(length(Z{i}))*tau*gv)];
    elseif is(UncertaintySet(i),'socp')
        u = sdpvar(length(Z{i})-1,1);
        [tau,coefftau] = polynomial(v,localizer_tau_degree);
        coeffs = [coeffs;coefftau];
        D = [D, sos([1 u']*Z{i}-tau*(1-u'*u))];
    elseif is(UncertaintySet(i),'elementwise')
         [tau,coefftau] = polynomial(v,localizer_tau_degree);
         coeffs = [coeffs;coefftau];
        D = [D, sos(Z{i}-tau*gv)];
   elseif is(UncertaintySet(i),'equality')
         [tau,coefftau] = polynomial(v,localizer_tau_degree);
         coeffs = [coeffs;coefftau];          
    end
end

[tau,coefftau] = polynomial(v,p_tau_degree);
p = -(trace(Zblock'*E) + b'*DecisionVariables + d) - tau*gv;
coeffs = [coeffs;coefftau];
SOSModel = compilesos([D, sos(p)],[],sdpsettings('sos.model',2,'sos.scale',0),coeffs);
