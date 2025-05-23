function [newConstraint] = normalChanceFilterConicFormulationInv(c,b,invPhi_Inverse)
%NORMALCHANCEFILTERCONICFORMULATIONINV Exponential cone formulation of the
%disjoint contraints using 1/probit.

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
    t <= invPhi_Inverse];
end