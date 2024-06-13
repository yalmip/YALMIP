function newConstraint = normalChanceFilter(b,c,distribution,gamma,w,options)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
if isa(gamma,'sdpvar') && strcmpi(options.chance.expcone,'yes') 
    a = 4.274819565452955e-01;
    b = 5.101801418186901e+06;
    k = 6.064075380897068e+00;
    c = -1.655121254255513e+00;
    Phi_Inverse = -a*lambertw(b*gamma)+k+c*gamma;
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
    newConstraint = b + c'*theMean >= Phi_Inverse*norm_callback(e);    
else
    newConstraint =  b + c'*theMean >= Phi_Inverse*norm(e);
end