function newConstraint = nonsymmetricUnivariateChanceFilter(b,c,distribution,gamma,options,funcs,x)

% General case must be handled via characteristic framework
if (length(c) > 1) || isa(c,'sdpvar') || any(strcmpi(options.chance.characteristic,{'yes','on'})) || isequal(options.chance.characteristic,1)
    newConstraint = [characteristic_cdf(x,funcs,distribution) >= 1-gamma];
    return
end

if c>=0
    % Such as b+c*w >= 0 i.e. -b/c <= w
    newConstraint = -b <= c*icdf(distribution.parameters{1},gamma,distribution.parameters{2:end});
else
    newConstraint = b >= -c*icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end});
end
   