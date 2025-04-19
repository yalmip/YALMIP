function newConstraint = normalChanceFilter(b,c,distribution,gamma,w,options,isDisjointProblem)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};

if isa(c,'double') && (~strcmpi(options.chance.expcone,''))
    % c is constant, so no reason really to use log tricks etc?
    % for now, just revert to standard expcone approximation
    options.chance.expcone = 'yes';
end
  
if strcmpi(options.chance.expcone,'yes')
    if isa(c,'sdpvar')
        error('Cannot have decision variables multplying uncertainty when using expcone approximation of inverse cdf')
    end
    Phi_Inverse = normalChanceFilterConicApproximation(gamma);
elseif strcmpi(options.chance.expcone,'root') && isDisjointProblem
    rootPhi_Inverse = normalChanceFilterConicApproximationRoot(gamma);
elseif strcmpi(options.chance.expcone,'log') && isDisjointProblem
    logPhi_Inverse = normalChanceFilterConicApproximationLog(gamma);
elseif strcmpi(options.chance.expcone,'inv') && isDisjointProblem
    invPhi_Inverse = normalChanceFilterConicApproximationInv(gamma);
else
    % Just go for a general nonlinear model and hope for the best
    Phi_Inverse = icdf('normal',1-gamma,0,1);
end
if min(size(covariance))==1
    covariance = diag(covariance);
end
if isa(covariance,'sdpvar')
    error('Covariance cannot be an SDPVAR in normal distribution. Maybe you meant to use factorized covariance in ''normalf''');
end
e = chol(covariance)*c;

if strcmpi(options.chance.expcone,'root') && isDisjointProblem
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

elseif strcmpi(options.chance.expcone,'log') && isDisjointProblem
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

    newConstraint = [a + logPhi_Inverse <= t,
        expcone([t;1;z]),
        z == b];
elseif strcmpi(options.chance.expcone,'inv') && isDisjointProblem
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

elseif isa(Phi_Inverse,'sdpvar')
    % If we have a general nonlinear model, we should not use norm as it
    % is intended to be socp-represented. Avoid that by using the callback
    % version of norm
    % Note, this cannot happen when expcone is used, stopped above
    newConstraint = b + c'*theMean >= Phi_Inverse*norm_callback(e);
else
    newConstraint =  b + c'*theMean >= Phi_Inverse*norm(e);
end