function [newConstraint] = normalChanceFilterConicFormulationLog(c,b,gamma)
%NORMALCHANCEFILTERCONICFORMULATIONLOG Exponential cone formulation of the
%disjoint contraints using an approximation of log(probit).

if isa(gamma,'sdpvar')
    aa = -0.196671288384826;
    bb = 1.674000760107396e+03;
    cc = -2.899789120102472;
    kk = 1.283146603876050;
    logPhi_InverseApproximation = aa*lambertw(bb*gamma)+kk+cc*gamma;
else
    % The approximation will be the exact value if gamma is fixed
    logPhi_InverseApproximation = log(icdf('normal',1-gamma,0,1));
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
        a = a+sd(i)*log(r{i});
    end
else
    r = norm(c0);
    a = log(r);
end

sdpvar z t; % epigraph variables

newConstraint = [a + logPhi_InverseApproximation <= t,
    expcone([t;1;z]),
    z == b];
end