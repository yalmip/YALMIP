function [conicApproximation] = normalChanceFilterConicApproximationRoot(gamma)
%NORMALCHANCEFILTERCONICAPPROXIMATIONROOT Exponential conic representable
%approximation of the square root of probit function.

if isa(gamma,'sdpvar')
    aa = -0.163460523135549;
    bb = 1.996987289085205e+03;
    cc = -1.232492830923356;
    kk = 1.898392103622973;
    conicApproximation = aa*lambertw(bb*gamma)+kk+cc*gamma;
else
    % The approximation will be the exact value if gamma is fixed
    conicApproximation = sqrt(icdf('normal',1-gamma,0,1));
end
end

