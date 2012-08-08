function varargout = boundingbox(F,ops,x)
%BOUNDINGBOX Computes bounding box of a constraint object
%
% If only one output is requested, only the symbolic model is returned
% B = boundingbox(F)
%
% If three outputs are requrested, the symbolic model is appended
% [B,L,U] = boundingbox(F)
%
% A second argument can be used to specificy solver settings
% B = boundingbox(F,sdpsettings('solver','cplex'))
%
% A third argument can be used to obtain a bounding box in a particular 
% set of variables, i.e. the bounding box of a projection
% B = boundingbox(F,[],x)

% Author Johan Löfberg
% $Id: boundingbox.m,v 1.1 2004-12-08 00:07:15 johanl Exp $

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

% sol = solvesdp(F,[x;-x],ops);
% n = length(x);
% for i = 1:n
%     if sol.problem(i)==0
%         selectsolution(i);
%         L(i,1) = double(x(i));
%     else
%         L(i,1) = -inf;
%     end
%     if sol.problem(n+i)==0
%         selectsolution(length(x)+i);
%         U(i,1) = double(x(i));
%     else
%         U(i,1) = inf;
%     end
% end
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
    B = [B, (x(Lf) >= L(Lf)):'Finite lower bounds'];
end
if ~isempty(Uf)
    B = [B, (x(Uf) >= L(Uf)):'Finite upper bounds'];
end
%B = [L <= x <= U];
switch nargout
    case 0
        B
    case 1
        varargout{1} = B;
    case 2
        varargout{1} = B;
        varargout{2} = L;
    case 3
        varargout{1} = B;
        varargout{2} = L;
        varargout{3} = U;
    otherwise
end