function dfdx = jacobian(f,x)
% JACOBIAN Jacobian of scalar or vector
%
% J = JACOBIAN(p)    Jacobian w.r.t all variables in p
% J = JACOBIAN(p,x)  Jacobian w.r.t the SDPVAR variables x
%
% See also INT, HESSIAN, LINEARIZE, SDISPLAY

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
    otherwise
        error('Too many input arguments.');
end
