function newConstraint = normalfactorizedChanceFilter(b,c,distribution,confidencelevel);
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
gamma = icdf('normal',confidencelevel,0,1);
e = covariance*c;
if isa(gamma,'sdpvar')
    newConstraint = b + c'*theMean >= gamma*norm_callback(e);
else
    newConstraint = b + c'*theMean >= gamma*norm(e);
end