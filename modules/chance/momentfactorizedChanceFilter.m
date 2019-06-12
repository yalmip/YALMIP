function model = momentfactorizedChanceFilter(b,c,distribution,gamma,w,options)
% Chance filter for distribution only specified by mean and factorized variance
theMean    = distribution.parameters{2};
R          = distribution.parameters{3};

rho = sqrtm(1-gamma);
e = [R*c;b + c'*theMean];
if isa(rho,'sdpvar')
    model =  b + c'*theMean >= rho*norm_callback(e);
else
    model =  b + c'*theMean >= rho*norm(e);
end