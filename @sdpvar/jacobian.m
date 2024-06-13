function dfdx = jacobian(f,x,d)
% JACOBIAN Jacobian of scalar, vector, or matrix
%
% J = JACOBIAN(p)    Jacobian w.r.t all variables in p
% J = JACOBIAN(p,x)  Jacobian w.r.t the SDPVAR variables x
%
% See also INT, HESSIAN, LINEARIZE, SDISPLAY

[n,m] = size(f);
if min([n m]) > 1
    f = reshape(f,[],1);
end

switch nargin
    case 1
        dfdx = shadowjacobian(f);
    case 2
        if length(getvariables(x) > length(x(:)))
            % typical case is jacobian(f(X),X) with X a symmetric matrix
            % Shadowjacobian assumes elements of x are independent
            dfdx = shadowjacobian(f,recover(x));
            x_indep = recover(x);
            dfdx = map_to_original(dfdx,x,x_indep);
        else
            dfdx = shadowjacobian(f,x(:));
        end
    case 3
        if d == 1
            dfdx = shadowjacobian(f,x(:));
        else
            dfdx = jacobian(jacobian(f,x(:)),x(:),d-1);
        end
    otherwise
        error('Too many input arguments.');
end

if min([n m]) > 1
    dfdx = reshape(full(dfdx),n,m,[]);
end