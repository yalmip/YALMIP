function [newConstraint] = normalChanceFilterConicFormulationRoot(c,b,gamma)
%NORMALCHANCEFILTERCONICFORMULATIONROOT Exponential cone formulation of the
%disjoint contraints using sqrt(probit).

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
        a = a+sd(i)*inv(r{i});
    end

else
    r = norm(c0);
    a = inv(r);
end

sdpvar t y; % epigraph variable

% approximating and formulating
if isa(gamma,'sdpvar')
    aa = -0.099181619811223;
    bb = 2.746525371530841e+03;
    cc = -1.305888804476780;
    kk = 1.579829567439729;
    dd = -0.043478545831448;
    rootPhi_Inverse = aa*lambertw(bb*gamma)+kk+cc*gamma+dd*y;

    newConstraint = [norm([b-a,2*t]',2) <= b+a,
        rootPhi_Inverse <= t,
        expcone([y;1;gamma])];
else
    % The approximation will be the exact value if gamma is fixed
    rootPhi_Inverse = sqrt(icdf('normal',1-gamma,0,1));

    newConstraint = [norm([b-a,2*t]',2) <= b+a,
        rootPhi_Inverse <= t];
end
end