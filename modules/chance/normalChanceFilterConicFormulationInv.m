function [newConstraint] = normalChanceFilterConicFormulationInv(c,b,gamma)
%NORMALCHANCEFILTERCONICFORMULATIONINV Exponential cone formulation of the
%disjoint contraints using an approximation of 1/probit.

% probability(b(x) + c(x)'*w >= 0)...

% separate c0 and ci
data = getbase(c);
c0 = data(:,1);
ci = data(:,2:end);

% pre-processing
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

sdpvar t y; % epigraph variables

% approximating and formulating
if isa(gamma,'sdpvar')
    aa = 0.045463730850896;
    bb = 1.186500274690110e+03;
    cc = 2.676791341187682;
    kk = 0.374878600022031;
    dd = 0.012244142944388;

    invPhi_InverseApproximation = aa*lambertw(bb*gamma)+kk+cc*gamma+dd*y;

    newConstraint = [norm([b-t,2*a]) <= b + t,
        t <= invPhi_InverseApproximation,
        expcone([y;1;gamma])];
else
    % The approximation will be the exact value if gamma is fixed
    invPhi_InverseApproximation = inv(icdf('normal',1-gamma,0,1));

    newConstraint = [norm([b-t,2*a]) <= b + t,
        t <= invPhi_InverseApproximation];
end
end