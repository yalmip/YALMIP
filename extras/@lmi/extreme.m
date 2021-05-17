function V = extreme(P)
% extreme Computes vertices (very rudimentary implementation)
%
% V = extreme(F)
%
% F : Constraint object defining a polytope
%
% See also CHEBYBALL, BOUNDINGBOX

[p,recoverdata,solver,diagnostic,F] = compileinterfacedata(P,[],[],[],sdpsettings,0);
x = recover(depends(P));

if anyCones(p.K) || any(p.variabletype) || any(p.binary_variables) || any(p.integer_variables)
  error('EXTREME can only be applied to LP-representable constraints.')
end

if p.K.f > 0
    error('EXTREME currently not working with equalities')
end

A = -p.F_struc(:,2:end);% Hx < K
b = p.F_struc(:,1);
if ~all(b) > 0
    [x_c,R] = chebyball(A*x <= b);
else
    x_c = zeros(size(A,2),1);
end

V = vertexenumerate(full(A),full(b),x_c);