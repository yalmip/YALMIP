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

coeffs = [DecisionVariables;coeffs(:)]
Zblock = blkdiag(Z{:});

D = [];
for i = 1:length(F)
    D = [D, coefficients(trace(Zblock'*F{i})-c(i),v)==0];
end

[trivialFixed,thevalue] = fixedvariables(D);
while ~isempty(trivialFixed) && length(D)>0
    D = replace(D,trivialFixed,thevalue);
  %  if ~isempty(D)
  %      D = 
  %  D = sdpvar(replace(D,trivialFixed,thevalue))==0;
    for i = 1:length(Z)
        Z{i} = replace(Z{i},trivialFixed,thevalue);
    end
    Zblock = replace(Zblock, trivialFixed,thevalue);
    if length(D)>0
        [trivialFixed,thevalue] = fixedvariables(D);
    end
end


% At this point Z is a function of v where v was used to scalarize the
% uncertain constraint. Now we must ensure Z{i}(v) in cone
gv = (1-v'*v);
for i = 1:length(Z)       
    if is(UncertaintySet(i),'sdp')
        % We use the matrix sos approach
        [tau,coefftau] = polynomial(v,localizer_tau_degree);
        coeffs = [coeffs;coefftau];
        D=[D,sos(Z{i}-eye(length(Z{i}))*tau*gv)];
    elseif is(UncertaintySet(i),'socp')
        % To get a SOS condition on dual Z{i}(v) in socp, we have to
        % introduce a new variable to scalarize the socp
        u = sdpvar(length(Z{i})-1,1);
        [tau,coefftau] = polynomial(v,localizer_tau_degree);
        coeffs = [coeffs;coefftau];
        D = [D, sos([1 u']*Z{i}-tau*(1-u'*u))];        
    elseif is(UncertaintySet(i),'elementwise')
         [tau,coefftau] = polynomial(v,localizer_tau_degree);
         coeffs = [coeffs;coefftau];
        D = [D, sos(Z{i}-tau*gv)];
   elseif is(UncertaintySet(i),'equality')
       % No constraints on dual         
    end
end

[tau,coefftau] = polynomial(v,p_tau_degree);
p = -(trace(Zblock'*E) + b'*DecisionVariables + d) - tau*gv;
coeffs = [coeffs;coefftau];
SOSModel = compilesos([D, sos(p)],[],sdpsettings('sos.model',2,'sos.scale',0),coeffs);

function [z,val] = fixedvariables(D)

Base = getbase(sdpvar(D));
A = -Base(:,2:end);
b = Base(:,1);
v = getvariables(D);
z = [];
val = [];
for i = 1:size(A,1)
    j = find(A(i,:));
    if length(j)==1
        z = [z v(j)];
        val = [val b(i)/A(i,j)];
    end
end
z = recover(z);
