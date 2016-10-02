function newConstraint = normalChanceFilter(b,c,distribution,confidencelevel,w,options)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
gamma = icdf('normal',confidencelevel,0,1);
if min(size(covariance))==1
    covariance = diag(covariance);
end
if isa(covariance,'sdpvar')
    error('Covariance cannot be an SDPVAR in normal distribution. Maybe you meant to use factorized covariance in ''normalf''');
end
e = chol(covariance)*c;
if isa(gamma,'sdpvar')
    newConstraint = b + c'*theMean >= gamma*norm_callback(e);
else
    newConstraint =  b + c'*theMean >= gamma*norm(e);
end