function [newConstraint] = normalChanceFilterConicFormulation(b,c,e,theMean,gamma)
%NORMALCHANCEFILTERCONICAPPROXIMATION Exponential cone formulation using
%an approximation of the probit function.

sdpvar y z % epigraph variable for asymptotic bahavior

% One upper bound...
if isa(gamma,'sdpvar')
    
    % Approximation to address the asymptotic behavior of probit
    cc = -1.325786346642915e+03;
    kk = 2.042472559346282;
    dd = -0.196422986236589;
    Phi_InverseApproximationTail = kk+cc*gamma+dd*y;
    
    % Approximation for higher risks
    aa = -0.499492956059166;
    bb = 8.082867432374761e+03;
    cc = -1.475743096725997;
    kk =  3.965651977413067;
    Phi_InverseApproximation = aa*lambertw(bb*gamma)+kk+cc*gamma;

    newConstraint = [b + c'*theMean >= z*norm(e),
        z >= Phi_InverseApproximation,
        z >= Phi_InverseApproximationTail,
        expcone([y;1;gamma])];

else
    % The approximation will be the exact value if gamma is fixed
    Phi_InverseApproximation = icdf('normal',1-gamma,0,1);
    newConstraint =  b + c'*theMean >= Phi_InverseApproximation*norm(e);
end
end

