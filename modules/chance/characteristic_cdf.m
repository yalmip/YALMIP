function varargout = characteristic_cdf(varargin)
% Internal function to supply numerical evaluation of CDF and its
% derivative based on characteristic function

switch class(varargin{1})
    
    case 'double'
        % Asked to evaluate cdf at point x
        % The function is Probability( h(x) + g'(x)*w <= 0 )
        % i.e.  Probability( z <= -h(x)) where z = g'(x)*w
        x = varargin{1}(:);
        funcs = varargin{2};
        distribution = varargin{3};
        
        phi = distribution.characteristicfunction;
        
        % Evaluate terms at x
        g = funcs.g(x);
        h = funcs.h(x);
        
        % Define char. func for linear combination
        phi = @(t) prod(phi(g(:)*t),1);        
        
        % Compute the cdf
        varargout{1} = compute_cdf_using_phi(-h, phi);
        
    case 'sdpvar'       
        varargout{1} = yalmip('define',mfilename,varargin{:});    
        
    case 'char'
        
        operator = struct('convexity','none','monotonicity','increasing','definiteness','positive','model','callback');        
        operator.range = [0 1];
        
        % Create a derivative callback 
        funcs = varargin{4};
        distribution = varargin{5};
        phi = distribution.characteristicfunction;
        dphi = distribution.characteristicfunction_derivative;
        % Create a function which computes gradient at x
        operator.derivative = @(x)compute_dcdf_using_phi(x,funcs.h,funcs.dh,funcs.g,funcs.dg,phi,dphi);
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
end

function cdf = compute_cdf_using_phi(y,phi)
% Perform the inverse Fourier transform to obtain the CDF
% for Probability(z <= y) where z has characterstic function phi(t)

% Here, brutally naive implementation. This is really ripe for issues due
% to singularities, so it should be guarded somehow with known value of the
% limit of phi(t)/t as t->0, and integral approximation is without any
% thought
integrand = @(t) imag(phi(t) .* exp(-1i * t * y)./t);
cdf =  .5-integral(integrand,0, 100)/pi;


function dcdf = compute_dcdf_using_phi_finite_difference(x,h,dh,g,dg, phi,dphi)
h0 = h(x);
dh0 = dh(x);
g0 = g(x);
dg0 = dg(x);


% Silly numerical differentiation (which is what a solver does anyway)
phi0 = @(t) prod(phi(g0(:)*t),1);  
cdf0 = compute_cdf_using_phi(-h0,phi0);
eps = 1e-8;
dcdf = [];
for k = 1:length(x)
    x_ = x;x_(k) = x_(k)+eps;
    g_ = g(x_);
    h_ = h(x_);    
    phi_ = @(t) prod(phi(g_(:)*t),1);  
    cdf_ = compute_cdf_using_phi(-h_,phi_);
    dcdf = [dcdf;(cdf_-cdf0)/eps];
end


%% -------------------------------------------------------------------------------------- 
function dcdf = compute_dcdf_using_phi(x,h,dh,g,dg,phi,dphi)

h0 = h(x);
dh0 = dh(x);
g0 = g(x);
dg0 = dg(x);

% compute phi_z(t)
phi_z = @(t) prod(phi(g0(:)*t),1);

% compute exponential part
exp_ith = @(t) exp(1i*t*h0);

% compute f_z(-h0)
integrand1 = @(t) real(exp_ith(t) .* phi_z(t));

% compute upper limit of t
epsilon = 1e-6;
t_low = 0;
t_high = 10;  
max_trial = 20; 
for k = 1:max_trial
    if abs(phi_z(t_high)) < epsilon
        break;
    else
        t_high = t_high*2; 
    end
end
if abs(phi_z(t_high)) >= epsilon
    error('cannot find suitable t_high for phi_z(t) < epsilon，please increase max_trial or adjust initial value');
end
options = optimset('TolX', 1e-3, 'Display', 'off');
t_max = fzero(@(t) abs(phi_z(t)) - epsilon, [t_low,t_high], options);
pdf_val = (1/pi)*integral(integrand1,0,t_max,'RelTol',1e-4,'AbsTol',epsilon);

% compute dphi_z/dg
num_j = length(g0);
dphi_p1 = @(t) [];
for j = 1:num_j
    dphi_p1 = @(t) [dphi_p1(t),dphi(g0(j)*t)];
end
integrand2 = @(t) imag(exp_ith(t) .* phi_z(t) .* diag(dphi_p1(t)) ./ phi(g0(:)*t));
terms = (-1/pi)*integral(integrand2,0,t_max,'RelTol',1e-4,'AbsTol',epsilon,'ArrayValued',true);

% commpute the whole derivative
dcdf = (-pdf_val.*dh0') + terms'*dg0;

dcdf_check = compute_dcdf_using_phi_finite_difference(x,h,dh,g,dg, phi,dphi);
[dcdf(:) dcdf_check(:)]

