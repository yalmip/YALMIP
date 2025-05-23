function [conicApproximation] = normalChanceFilterConicApproximationLog(gamma)
%NORMALCHANCEFILTERCONICAPPROXIMATIONLOG Exponential conic representable
%approximation of the log of the probit function.

if isa(gamma,'sdpvar')
    aa = -0.196671288384826;
    bb = 1.674000760107396e+03;
    cc = -2.899789120102472;
    kk = 1.283146603876050;
    conicApproximation = aa*lambertw(bb*gamma)+kk+cc*gamma;
else
    % The approximation will be the exact value if gamma is fixed
    conicApproximation = log(icdf('normal',1-gamma,0,1));
end
end

