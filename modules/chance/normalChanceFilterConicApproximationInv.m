function [conicApproximation] = normalChanceFilterConicApproximationInv(gamma)
%NORMALCHANCEFILTERCONICAPPROXIMATIONINV Exponential conic representable
%approximation of the inverse of the probit function.

if isa(gamma,'sdpvar')
    aa = 0.050229622348771;
    bb = 7.573772400040184e+04;
    cc = 2.732774841525416;
    kk = 0.150527341988232;
    conicApproximation = aa*lambertw(bb*gamma)+kk+cc*gamma;
else
    % The approximation will be the exact value if gamma is fixed
    conicApproximation = inv(icdf('normal',1-gamma,0,1));
end
end

