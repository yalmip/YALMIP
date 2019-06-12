function model = momentChanceFilter(b,c,distribution,gamma,w,options)
% Chance filter for distribution only specified by mean and variance
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};

S = chol(covariance);
rho = sqrtm(1-gamma);
e = [S*c;b + c'*theMean];
if isa(rho,'sdpvar')
    model =  b + c'*theMean >= rho*norm_callback(e);
else
    model =  b + c'*theMean >= rho*norm(e);
end