function V = vertex(P,xin)
% VERTEX Computes vertices (very rudimentary implementation)
%
% V = VERTEX(P)
%
% P : Constraint object defining a polytope
%
% With a second argument x the projection of the vertices
% to the x-space is returned.
%
% See also CHEBYBALL, BOUNDINGBOX, POLYTOPE

[p,recoverdata,solver,diagnostic,F] = compileinterfacedata(P,[],[],[],sdpsettings,0);

if anyCones(p.K) || any(p.variabletype) || any(p.binary_variables) || any(p.integer_variables)
  error('VERTEX can only be applied to LP-representable constraints.')
end

A = -p.F_struc(p.K.f + 1:end,2:end);
b = p.F_struc(p.K.f + 1:end,1);
    Aeq = -p.F_struc(1:p.K.f,2:end);
    beq = p.F_struc(1:p.K.f,1);
if ~isempty(Aeq)
    % A*x <= b, Aeq*x = Beq, x = x0 + N*z
    % A*N*z < b-A*x0
    x0 = Aeq\beq;
    N = null(full(Aeq));
    b = b-A*x0;
    A = A*N;    
end

if all(b > 0)    
	x_c = zeros(size(A,2),1);
else
    x = sdpvar(size(A,2),1);
    if isempty(Aeq)
        [x_c,R] = chebyball(A*x <= b);
    else        
        [x_c,R] = chebyball([A*x <= b]);
    end
end
if R < 0
    error('The polytope appears empty');
end

V = vertexenumerate(full(A),full(b),x_c);
if ~isempty(Aeq)
    V = repmat(x0,1,size(V,2)) + N*V;
end

if nargin == 2
    for i = 1:length(xin)
        xi(i) = getvariables(xin(i));
        xi(i) = find(xi(i) == recoverdata.used_variables);        
    end
    V = V(xi,:);
    V = unique(V','rows')';
end
    