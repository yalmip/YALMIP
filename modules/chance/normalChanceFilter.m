function newConstraint = normalChanceFilter(b,c,distribution,gamma,w,options)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
isDisjointProblem = 1;% FAKEIT
if isa(gamma,'sdpvar') && strcmpi(options.chance.expcone,'yes')     
    if isa(c,'sdpvar')
        error('Cannot have decision variables multplying uncertainty when using expcone approximation of inverse cdf')
    end
    % One upper bound...
    aa = 0.499492956059166;
    bb = 8.082867432374761e+03;
    cc = -1.475743096725997;
    kk =  3.965651977413067;   
    Phi_Inverse = -aa*lambertw(bb*gamma)+kk+cc*gamma;
elseif isa(gamma,'sdpvar') && strcmpi(options.chance.expcone,'root') && isDisjointProblem     
    
    RootPhi_Inverse = 0;
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
if isa(Phi_Inverse,'sdpvar') && strcmpi(options.chance.expcone,'root') && isDisjointProblem     
    % probability(b(x) + c(x)'*w >= 0)...
    data = getbase(c);
    c0 = data(:,1);
    ci = data(:,2:end);
    
    
    
elseif isa(Phi_Inverse,'sdpvar')
    % If we have a general nonlinear model, we should not use norm as it 
    % is intended to be socp-represented. Avoid that by using the callback
    % version of norm
    % Note, this cannot happen when expcone is used, stopped above
    newConstraint = b + c'*theMean >= Phi_Inverse*norm_callback(e);    
else
    newConstraint =  b + c'*theMean >= Phi_Inverse*norm(e);
end