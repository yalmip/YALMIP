function varargout = chebyball(F,ops)
%CHEBYBALL Computes Chebyshev ball of a constraint object
%
% If two outputs are requested, the numerical data is returned
% [xc,R] = chebyball(F)
%
% If three outputs are requrested, the symbolic model (x-xc)'(x-xc)<R^2 is
% appended
% [xc,R,C] = chebyball(F)
%
% If only one output is requested, only the symbolic constraint is returned
% C = chebyball(F)

% Author Johan Löfberg
% $Id: chebyball.m,v 1.1 2004-12-08 00:07:15 johanl Exp $

[model,recoverdata,diagnostic,p] = export(F,[],[],[],[],0);
if p.K.q(1) > 0 | p.K.s(1) > 0 | any(p.variabletype) | ~isempty(p.binary_variables) | ~isempty(p.integer_variables)
    error('Polytope can only be applied to linear elementwise constraints.')
end

A = -p.F_struc(:,2:end);
b = p.F_struc(:,1);

x = recover(p.used_variables);
r = sdpvar(1);

if nargin < 2
    ops = sdpsettings('verbose',0);
end

sol = solvesdp(A*x+r*sqrt(sum(A.^2,2))<=b,-r,ops);

xc = double(x);
R = double(r);
C = (x-xc)'*(x-xc) <= R^2;
if sol.problem == 1;
    R = 0;
elseif sol.problem == 2
    R = inf;
end

switch nargout
    case 0
    case 1
        varargout{1} = C;
    case 2
        varargout{1} = xc;
        varargout{2} = R;
    case 3
        varargout{1} = xc;
        varargout{2} = R;
        varargout{3} = C;
    otherwise
end