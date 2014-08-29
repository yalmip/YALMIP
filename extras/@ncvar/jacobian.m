function dfdx = jacobian(f,x)
% JACOBIAN Jacobian of scalar or vector
%
% J = JACOBIAN(p)    Jacobian w.r.t all variables in p
% J = JACOBIAN(p,x)  Jacobian w.r.t the SDPVAR variables x
%
% See also SDPVAR, HESSIAN, LINEARIZE

switch nargin    
    case 1        
        dfdx = shadowjacobian(f);
    case 2
        dfdx = shadowjacobian(f,x);
    otherwise
        error('Too many input arguments.');
end
