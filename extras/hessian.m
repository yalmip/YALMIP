function d2fdx2 = hessian(f,x)
% HESSIAN Hessian of scalar polynomial SDPVAR object
%
% J = HESSIAN(p)    Hessian w.r.t all variables in p
% J = HESSIAN(p,x)  Hessian w.r.t the SDPVAR variables x
%
% See also INT, JACOBIAN, LINEARIZE, SDISPLAY

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

if any(ismember([getvariables(f) depends(f)],yalmip('extvariables')))
    error('Hessian is only applicable to polynomial expressions');
end

if isa(f,'double')
    d2fdx2 = zeros(length(x));
    return
end

d2fdx2 = jacobian(jacobian(f,x)',x);
