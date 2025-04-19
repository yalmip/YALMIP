function [newConstraint] = normalChanceFilterConicFormulationRoot(c,b,rootPhi_Inverse)
%NORMALCHANCEFILTERCONICFORMULATIONROOT Exponential cone formulation of the
%disjoint contraints using sqrt(probit).

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
        a = a+sd(i)*inv(r{i});
    end

else
    r = norm(c0);
    a = inv(r);
end

sdpvar t; % epigraph variable

newConstraint = [norm([b-a,2*t]',2) <= b+a,
    rootPhi_Inverse <= t];
end

