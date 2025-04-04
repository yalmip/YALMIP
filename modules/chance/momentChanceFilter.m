function model = momentChanceFilter(b,c,distribution,gamma,w,options)
% Chance filter for distribution only specified by mean and variance
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
issymmetric = 0;
isunimodal = 0;
for k = 4:length(distribution.parameters)
    if isequal(distribution.parameters{k},'symmetric')
        issymmetric = 1;
    elseif isequal(distribution.parameters{k},'symmetric-unimodal')
        isunimodal = 1;
    end
end

R = chol(covariance);
rho = sqrt((1-gamma)/gamma);
if issymmetric
    rho = sqrt(1/(2*gamma));
end
if isunimodal
    rho = sqrt((2/3)*log(1./gamma));
end
    
e = R*c;
if isa(rho,'sdpvar')
    model =  b + c'*theMean >= rho*norm_callback(e);
else
    model =  b + c'*theMean >= rho*norm(e);
end