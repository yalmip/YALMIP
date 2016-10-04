function model = momentfactorizedChanceFilter(b,c,distribution,confidencelevel,w,options)
% Chance filter for distribution only specified by mean and factorized variance
theMean    = distribution.parameters{2};
R          = distribution.parameters{3};

gamma = sqrtm(confidencelevel);
e = [R*c;b + c'*theMean];
if isa(gamma,'sdpvar')
    model =  b + c'*theMean >= gamma*norm_callback(e);
else
    model =  b + c'*theMean >= gamma*norm(e);
end