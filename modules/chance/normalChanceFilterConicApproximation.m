function [conicApproximation] = normalChanceFilterConicApproximation(gamma)
%NORMALCHANCEFILTERCONICAPPROXIMATION Exponential conic representable
%approximation of the probit function.

% One upper bound...
if isa(gamma,'sdpvar')
    aa = 0.499492956059166;
    bb = 8.082867432374761e+03;
    cc = -1.475743096725997;
    kk =  3.965651977413067;
    conicApproximation = -aa*lambertw(bb*gamma)+kk+cc*gamma;
else
    % The approximation will be the exact value if gamma is fixed
    conicApproximation = icdf('normal',1-gamma,0,1);
end
end

