function varargout = boundingbox(F,ops,x)
%BOUNDINGBOX Computes bounding box of a constraint
%
% If only one output is requested, only the symbolic model is returned
%  B = boundingbox(F)
%
% If three outputs are requested, numerical values are returned too
%  [B,L,U] = boundingbox(F)
%
% A second argument can be used to specificy solver settings
%  B = boundingbox(F,sdpsettings('solver','cplex'))
%
% A third argument can (should!) be used to obtain a bounding box in a
% particular set of variables, i.e. the bounding box of a projection.
% Unless you specify which variables you are computing the bounding box
% w.r.t, it will be very hard for you to understand which variables the
% values in L and U relate to
%  [B,L,U] = boundingbox(F,[],x)
% B will now be the box [L <= x <= U] (infinite bounds not included)

if nargin < 3
    x = recover(depends(F));
else
    x = x(:);
    if ~isa(x,'sdpvar')
        error('The third argument should be an SDPVAR obeject');
    end
end

if nargin < 2
    ops = sdpsettings('verbose',0);
end
if isempty(ops)
    ops = sdpsettings('verbose',0);
end
if ~isa(ops,'struct')
    error('The second argument should be an SDPSETTINGS struct (or empty)');
end

sol = solvesdp(F,[x;-x],ops);
n = length(x);
for i = 1:n
    xi = x(i);
    if isa(xi,'sdpvar')
        if sol.problem(i)==0
           % selectsolution(i);
            L(i,1) = double(xi,i);
        else
            L(i,1) = -inf;
        end
        if sol.problem(n+i)==0
           % selectsolution(n+i);
            U(i,1) = double(xi,n+i);
        else
            U(i,1) = inf;
        end
    else
        L(i,1) = xi;
        U(i,1) = xi;
    end
end

% Only add finite bounds
Lf = find(~isinf(L));
Uf = find(~isinf(U));
B = [];
if ~isempty(Lf)
    xLf = x(Lf);
    if isa(xLf,'sdpvar')
        B = [B, (xLf >= L(Lf)):'Finite lower bounds'];
    end
end
if ~isempty(Uf)
    xUf = x(Uf);
    if isa(xUf,'sdpvar')
        B = [B, (xUf <= U(Uf)):'Finite upper bounds'];
    end
end
varargout{1} = B;
varargout{2} = L;
varargout{3} = U;