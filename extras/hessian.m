function d2fdx2 = hessian(f,x)
% HESSIAN Hessian of scalar SDPVAR object
%
% J = HESSIAN(p)    Hessian w.r.t all variables in p
% J = HESSIAN(p,x)  Hessian w.r.t the SDPVAR variables x
%
% See also SDPVAR, JACOBIAN, LINEARIZE

% Author Johan Löfberg
% $Id: hessian.m,v 1.2 2004-07-02 08:17:31 johanl Exp $


if nargin==1
    if isa(f,'sdpvar')
        x = recover(depends(f));
    else
        x = 0;
    end
else
    if length(getvariables(x))<length(x)
        error('x should be a vector of scalar independant variables');
    end
end

if prod(size(f))>1
   error('Hessian only defined for scalars.')
end 

if isa(f,'double')
    d2fdx2 = zeros(length(x));
    return
end

d2fdx2 = jacobian(jacobian(f,x)',x);
