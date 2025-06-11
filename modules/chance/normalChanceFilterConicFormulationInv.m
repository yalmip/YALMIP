function [newConstraint] = normalChanceFilterConicFormulationInv(c,b,gamma)
%NORMALCHANCEFILTERCONICFORMULATIONINV Exponential cone formulation of the
%disjoint contraints using an approximation of 1/probit.

if isa(gamma,'sdpvar')
    aa = 0.050229622348771;
    bb = 7.573772400040184e+04;
    cc = 2.732774841525416;
    kk = 0.150527341988232;
    invPhi_InverseApproximation = aa*lambertw(bb*gamma)+kk+cc*gamma;
else
    % The approximation will be the exact value if gamma is fixed
    invPhi_InverseApproximation = inv(icdf('normal',1-gamma,0,1));
end

% probability(b(x) + c(x)'*w >= 0)...

% separate c0 and ci
data = getbase(c);
c0 = data(:,1);
ci = data(:,2:end);

a = 0;
if isa(c,'sdpvar')
    % compute ri = ||c0+ci|| corresponding to the feedbacks
    for i=1:size(ci,2)
        r{i} = norm(c0+ci(:,i));
    end
    sd = recover(c); % recover the binary variables
    for i=1:size(ci,2)
        a = a+sd(i)*sqrt(r{i});
    end
else
    r = norm(c0);
    a = sqrt(r);
end

sdpvar t z; % epigraph variables

newConstraint = [norm([b-t,2*a]) <= b + t,
    t <= invPhi_InverseApproximation];
end