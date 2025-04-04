function newConstraint = normalfactorizedChanceFilter(b,c,distribution,gamma,w,options);
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
Phi_Inverse = icdf('normal',1-gamma,0,1);
e = covariance*c;
if isa(Phi_Inverse,'sdpvar')
    newConstraint = b + c'*theMean >= Phi_Inverse*norm_callback(e);
else
    newConstraint = b + c'*theMean >= Phi_Inverse*norm(e);
end