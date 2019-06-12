function newConstraint = normalChanceFilter(b,c,distribution,gamma,w,options)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
Phi_Inverse = icdf('normal',1-gamma,0,1);
if min(size(covariance))==1
    covariance = diag(covariance);
end
if isa(covariance,'sdpvar')
    error('Covariance cannot be an SDPVAR in normal distribution. Maybe you meant to use factorized covariance in ''normalf''');
end
e = chol(covariance)*c;
if isa(Phi_Inverse,'sdpvar')    
    newConstraint = b + c'*theMean >= Phi_Inverse*norm_callback(e);    
else
    newConstraint =  b + c'*theMean >= Phi_Inverse*norm(e);
end