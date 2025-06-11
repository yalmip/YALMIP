function [newConstraint] = normalChanceFilterConicFormulation(b,c,e,theMean,gamma)
%NORMALCHANCEFILTERCONICAPPROXIMATION Exponential cone formulation using
%an approximation of the probit function.

sdpvar y % epigraph variable for asymptotic bahavior

% One upper bound...
if isa(gamma,'sdpvar')
    aa = -0.261576448571477;
    bb = 3.488899999998545e+05;
    cc = -1.868722313969487;
    kk =  3.405748225416125;
    dd = -0.128360154389090;

    Phi_InverseApproximation = aa*lambertw(bb*gamma)+kk+cc*gamma+dd*y;
   
    newConstraint = [b + c'*theMean >= Phi_InverseApproximation*norm(e),
        expcone([y;1;gamma])];
else
    % The approximation will be the exact value if gamma is fixed
    Phi_InverseApproximation = icdf('normal',1-gamma,0,1);
    newConstraint =  b + c'*theMean >= Phi_InverseApproximation*norm(e);
end
end

