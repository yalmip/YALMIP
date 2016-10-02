function model = momentChanceFilter(b,c,distribution,confidencelevel)
% Chance filter for distribution only specified by mean and variance
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
X = covariance + theMean*theMean';
S = chol(X);
gamma = sqrtm(confidencelevel);
e = [S*c+b*(inv(S')*theMean);b*sqrtm(1-theMean'*inv(X)*theMean)];
model =  b + c'*theMean >= gamma*norm_callback(e);