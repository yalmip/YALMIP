function newConstraint = normalChanceFilter(b,c,distribution,gamma,w,options)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
if isa(gamma,'sdpvar') && strcmpi(options.chance.expcone,'yes')     
    if isa(c,'sdpvar')
        error('Cannot have decision variables multplying uncertainty when using expcone approximation of inverse cdf')
    end
    % One upper bound...
    aa = 4.274819565452955e-01;
    bb = 5.101801418186901e+06;
    kk = 6.064075380897068e+00;
    cc = -1.655121254255513e+00;
    Phi_Inverse = -aa*lambertw(bb*gamma)+kk+cc*gamma;
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
if isa(Phi_Inverse,'sdpvar')  
    % If we have a general nonlinear model, we should not use norm as it 
    % is intended to be socp-represented. Avoid that by using the callback
    % version of norm
    % Note, this cannot happen when expcone is used, stopped above
    newConstraint = b + c'*theMean >= Phi_Inverse*norm_callback(e);    
else
    newConstraint =  b + c'*theMean >= Phi_Inverse*norm(e);
end