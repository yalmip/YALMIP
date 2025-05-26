function distribution = setupCharacteristics(name, distribution)

% REMEMBER, these have to instantiated in modelchance also
switch name
    case 'exponential'
        % N.B parameterized in scale parameter mu = 1/lambda in
        % Statistics toolbox and thus our parameter
        phi = @(t,mu)1./(1-t*1i.*mu(:));
        dphi = @(t,mu)1i*(mu(:))./(1-1i*t.*mu(:)).^2;
        distribution.characteristicfunction = phi;
        distribution.characteristicfunction_derivative = dphi;
        distribution.characteristicfunction_relativederivative = @(t,mu) dphi(t,mu)./phi(t,mu);
    case 'laplace'
        distribution.characteristicfunction = @(t,mu,b)exp(1i*mu(:)*t)./(1+b(:).^2.*t.^2);
        % FIXME!
        distribution.characteristicfunction_derivative = [];
        distribution.characteristicfunction_relativederivative = [];
    case 'logistic'
        phi = @(t,mu,s)(exp(1i*mu(:).*t))./guarded_sinhc(pi*s(:).*t);
        dphi = @(t,mu,s) guarded_logistic_derivative(t,mu,s);
        distribution.characteristicfunction = phi;
        distribution.characteristicfunction_derivative = dphi;
        distribution.characteristicfunction_relativederivative = @(t,mu,s) dphi(t,mu,s)./phi(t,mu,s);
        
    case 'uniform'
        phi = @(t,a,b)(guarded_expdiv(b(:).*t,a(:).*t,t.*(b(:)-a(:))));
        dphi = @(t,a,b) (1 ./ (1i*(b(:)-a(:)) .* t.^2)) .* (t.*(1i*b(:).*exp(1i * b(:) .* t) - 1i*a(:) .* exp(1i * a(:) .* t)) - (exp(1i * b(:) .* t) - exp(1i * a(:) .* t)));        
        distribution.characteristicfunction = phi;
        distribution.characteristicfunction_derivative = dphi;        
        distribution.characteristicfunction_relativederivative = @(t,a,b) dphi(t,a,b)./phi(t,a,b);

	case 'normal'
        phi = @(t,mu,sigma) exp(1i*t.*mu(:) - 0.5*sigma(:).^2.*t.^2);
        dphi = @(t,mu,sigma) exp(1i*mu(:) - sigma(:).^2.*t).*exp(1i*t.*mu(:) - 0.5*sigma(:).^2.*t);
        distribution.characteristicfunction = phi;
        distribution.characteristicfunction_derivative = dphi;        
        distribution.characteristicfunction_relativederivative = @(t,a,b) dphi(t,a,b)./phi(t,a,b);
        
    otherwise
        distribution.characteristicfunction = [];
        distribution.characteristicfunction_derivative
        distribution.characteristicfunction_derivative = [];
end

function y = guarded_sinhc(x)
y = sinh(x)./x;y(x==0)=1;

function y = guarded_expdiv(z1,z2,z3)
y = (exp(1i*z1)-exp(1i*z2))./(1i*z3);
y(z3==0) = 1;

function y = guarded_logistic_derivative(t,mu,s)

y = exp(1i*mu(:).*t).*(1i*mu(:).*(pi.*s(:).*t./(sinh(pi*s(:).*t))) + (pi*s(:).*(1./sinh(pi*s(:).*t))-pi^2*s(:).^2.*t./((tanh(pi*s(:).*t).*sinh(pi*s(:).*t)))));
if numel(t)==1 && t == 0
    % Simple case phi(t) for scalar t, so all elements to limit
    y = 1i*mu(:);
else
    % t is a matrix (i.e. a vectorized call computing several points at
    % once)
    % Two cases, a row vector t meaning vector phi(t)
    if size(t,1) == 1
        i = find(t == 0);
        y(:,i) = repmat(1i*mu,1,length(i));
    else
        [i,j] = find(t == 0);
        y(sub2ind(size(t),i,j)) = 1i*mu(i);
    end
end