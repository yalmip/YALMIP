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
        % First parameter is name, so remove that
        parameters = {distribution.parameters{2:end}};
        mixtureweights = distribution.mixture;
        
        % Evaluate terms at x
        g = funcs.g(x);
        h = funcs.h(x);
                      
        % Silly case Probability(h + 0*w <= 0)
        if nnz(g)==0
            varargout{1} = h <= 0;
            return
        end
        
        % Define char. func for linear combination z = g(x)^Tw, and insert
        % the numerical parameters (means, shapes etc)  into the function
        % definition of the characteristic functions
        if ~isa(parameters{1},'cell')
            if length(g) == 1
                phi_z = @(t) phi(g(:)*t,parameters{:});
            else
                phi_z = @(t) prod(phi(g(:)*t,parameters{:}),1);
            end
        else
            % This is a mixture!
            % Create mixture characteristic function...
            phi_mixture = @(t) createMixtureSum(phi,g(:)*t,mixtureweights,parameters{:});
            phi_z = @(t) prod(phi_mixture(t),1);
        end
        
        % Compute the cdf. We send the name of the distribution to make 
        % choices about integration method
        varargout{1} = compute_cdf_using_phi(-h, phi_z,distribution.parameters{1});
        
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
        % First parameter is name, so remove that
        parameters = {distribution.parameters{2:end}};
        mixtureweights = distribution.mixture;
        % Create a function which computes gradient at x
        operator.derivative = @(x)compute_dcdf_using_phi(x,funcs.h,funcs.dh,funcs.g,funcs.dg,phi,dphi,parameters,mixtureweights,distribution.parameters{1});
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        % Clear the memoization used to speed up stuff during optimization
        cacheData('clear')
    otherwise
end

function cdf = compute_cdf_using_phi(y,phi_z,distribution)
% Perform the inverse Fourier transform to obtain the CDF
% for Probability(z <= y) where z has characterstic function phi(t)
% Later this will be computed together with derivatives etc but this
% requires some global logic to sync those when solver wants stuff

% We cannot use rotation on uniform distributions
if ~any(ismember(distribution,{'uniform'}))  

    % Extract expected value (can be done analytically but hack for now)
    mu = imag((phi_z(1e-6)-phi_z(-1e-6))/(2*1e-6)); 
    
    % We essentially have e^(i*mu*t)*phi_0(t)*exp(i*(-y)*t)
    % Use contour deformation. First figure out a rotation to damp    
    if y-mu == 0
        alpha = 0;
    else
        alpha = sign(mu-y)*30*pi/180;
    end
    q = exp(1i * alpha);
    
    % We use a singularity removal to deal with 1/t and it put back
    % afterwards, phi_z(0) = 1 
    % The 'q' from dt and 'q' from 1/t cancel out in the ((phi-1)/t)dt term.    
    integrand = @(r) imag((phi_z(r*q) .* exp(-1i*y*r*q) - 1) ./ r);
    % Expontial is damped to decay via term exp(-sin(alpha)*(mu-y)*r)
    % so we only integrate over interval where this is numerically non-zero
    U = log(10^-20)/(sin(alpha)*(y-mu));
    if U > 1e4
        U = inf;
    end
    I_num = integral(integrand, 0, U, 'RelTol', 1e-12);
    % We integrated along a small circle at the origin when removing the
    % singulatity, the arc contribution is i * alpha 
    I_total = I_num + imag(1i * alpha);
    cdf = 0.5 - I_total/pi;
else
    % Vanilla integration. Has to be done more carefully
    integrand = @(t) imag(phi_z(t) .* exp(-1i * t * y)./t);
    I_builtin = integral(integrand,0, 100);
    cdf = min(1,max(0,.5 - I_builtin/pi));
end

function dFdx = compute_dcdf_using_phi(x,h,dh,g,dg,phi,dphi,parameters,mixtureweights,distribution)
% Compute derivative of Probability(h(x)+g(x)'*w <= 0)

% Evaluate the operators
h_x  = h(x);
dh_x = dh(x);
g_x  = g(x);
dg_x = dg(x);

% Create an md5 identifier of these arguments, used to cache some data in
% the integration routines over iterations
functionHash = makeArgHash({h,dh,g,dg,phi,dphi,parameters,mixtureweights});

% Create characteristic functions etc
if ~isa(parameters{1},'cell')
    phi_z   = @(t) prod(phi(g_x(:)*t,parameters{:}),1);
    phi_    = @(t) phi(g_x(:)*t,parameters{:});
    dphi_   = @(t) dphi(g_x(:)*t,parameters{:});
else
    % This is a mixture. The characteristic of the mixture has to be
    % created etc
    [phi_,dphi_,phi_z] = characteristic_mix(g_x,mixtureweights,parameters,phi,dphi);
end

% Compute f_z(-h0) using Gil-Pelaez
% this will be moved to be computed together with derivative instead but
% for now we do the double work to keep code simple. Also, for now we use
% built-in integral
exp_ith = @(t) exp(1i*t*h_x);
integrand_pdf = @(t) real(exp_ith(t) .* phi_z(t));
if nnz(dh_x)==0 || (nnz(g_x)==0)
    pdf_val = 0;
else    
    if ~any(ismember(distribution,{'uniform'}))      
        I_deformation  = integral_lazy_rotation(phi_z,h_x,30*pi/180,@(z)real(z));              
        pdf_val = max(0,real(I_deformation)/pi);     
    else        
        % Vanilla integration for uniform
        I_builtin = integral(integrand_pdf,0,100,'AbsTol',1e-6);
        pdf_val = max(0,real(I_builtin)/pi);                       
    end
end

if nnz(dg_x)==0 || nnz(g_x)==0
    dFdx = (-pdf_val.*dh_x');
else
    if ~any(ismember(distribution,{'uniform'}))    
        theta = (30*pi/180)*sign(h_x);
        z = exp(1i*theta);
        dFdg = ((-1/pi)*dcdfintegralEvaluator(@(t)phi_(t*z),@(t)dphi_(t*z),@(t)(z*exp_ith(t*z)),functionHash,dg_x,@(z)imag(z)));
    else
        dFdg = (-1/pi)*dcdfintegralEvaluator(phi_,dphi_,exp_ith,functionHash,dg_x,@(z)(imag(z)));
    end
    dFdx = (-pdf_val.*dh_x') + dFdg'*dg_x;
end

function GK715 = setupGK715
GK715.Nodes = [ ...
    -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, ...
    -0.7415311855993944, -0.5860872354676911, -0.4058451513773972, ...
    -0.2077849550078985, 0, 0.2077849550078985, ...
    0.4058451513773972, 0.5860872354676911, 0.7415311855993944, ...
    0.8648644233597691, 0.9491079123427585, 0.9914553711208126];
GK715.Weights_15 = [ ...
    0.02293532201052922, 0.06309209262997855, 0.1047900103222502, ...
    0.1406532597155259, 0.1690047266392679, 0.1903505780647854, ...
    0.2044329400752989, 0.2094821410847278, 0.2044329400752989, ...
    0.1903505780647854, 0.1690047266392679, 0.1406532597155259, ...
    0.1047900103222502, 0.06309209262997855, 0.02293532201052922];
GK715.Weights_7 = [ ...
    0.1294849661688697, 0.2797053914892767, 0.3818300505051189, ...
    0.4179591836734694, 0.3818300505051189, 0.2797053914892767, ...
    0.1294849661688697];

function C = createMixtureSum(f, g_t,mixtureweights,varargin)
% Characteristic functions sum_j alpha_j phi(g*t,parameters_j)
C = [];
% varargin are all distribution parameters, each containing a cell with
% the indvidual component parameters
% Every element can have different mixture weights, but the number of
% components is the same for all elements
n_mixtures = length(varargin{1});
C = 0;
for j = 1:n_mixtures
    % Extract parameters for this mixture component
    jthParameters = cellfun(@(c) c{j}, varargin, 'UniformOutput', false);
    % Evaluate characteristic function. Note that g(x)*t can be a row
    % vector as an integral evaluator might evaluate multiple points
    C = C + mixtureweights(:,j).*f(g_t,jthParameters{:});
end

function [phi_,dphi_,phi_z] = characteristic_mix(g0,mixtureweights,parameters,phi,dphi)
W = size(mixtureweights,2);
P = numel(parameters);
% Change the parameters cell from being structured by the parameter number
% to being structured by component number
pars_cell = cell(1,W);
for w = 1:W
    pars = cell(1,P);
    for p = 1:P
        v = parameters{p};
        if iscell(v)
            pars{p} = v{w};
        else
            pars{p} = v(w);
        end
    end
    pars_cell{w} = pars;
end
% Compute the characteristic function of the mixture distribution and its
% derivative
phi_ = @(t) phi_mix(t,g0,W,mixtureweights,pars_cell,phi);
dphi_ = @(t) dphi_mix(t,g0,W,mixtureweights,pars_cell,dphi);
phi_z = @(t) phi_z_mix(t,g0,W,mixtureweights,pars_cell,phi);

function phi_ = phi_mix(t,g0,W,mixtureweights,pars_cell,phi)
phi_ = 0;
gt = g0(:)*t;
for w = 1:W
    phi_ = phi_ + mixtureweights(:,w).*phi(gt,pars_cell{w}{:});
end

function dphi_ = dphi_mix(t,g0,W,mixtureweights,pars_cell,dphi)
dphi_ = 0;
gt = g0(:)*t;
for w = 1:W
    dphi_ = dphi_ + mixtureweights(:,w).*dphi(gt,pars_cell{w}{:});
end

function phi_z = phi_z_mix(t,g0,W,mixtureweights,pars_cell,phi)
phi_ = phi_mix(t,g0,W,mixtureweights,pars_cell,phi);
phi_z = prod(phi_,1);


function db = cacheData(delete)
persistent cache
if nargin > 0
    cache = [];
    return
end
if isempty(cache)
    cache = containers.Map('KeyType','char','ValueType','any');
end
db = cache;

function y = fixNaN(x)
% Replace NaN or Inf in real part with 0
re = real(full(x));
re(isnan(re)) = 0;
% Replace NaN or Inf in imaginary part with 0
im = imag(full(x));
im(isnan(im)) = 0;
y = complex(re, im);

function I = dcdfintegralEvaluator(phi_,dphi_,exp_ith,functionHash,dg_x,part)
% Compute integral on the transformed domain 0->1
% Could be various variants, but only one available now
Points = setupGK715;
N = 3;

if nargin < 6
    % Possibility to send operator to apply operator on resulting integral
    % terms
    part = @(z)z;
end

if nargin < 5
    % Not used at the moment, but it should be possible to exploit that we
    % know what the computed integrals are multiplied with
    dg_x = [];
end

% Protect from weird stuff (silently so a bit dangerous)
phi_   = @(t)fixNaN(phi_(t));
dphi_ = @(t)fixNaN(dphi_(t));


% We resuse subdivisions from earlier computations on the same integral
% computation (for other x though) in the hope that they are somewhat
% similiar over iterations
cache = cacheData();
try
    % Cache de-activated for now
    % oldSubDivisionData = cache(functionHash);
    % D = unique(oldSubDivisionData);
    oldSubDivisionData = [];
    D = linspace(0,1,N);
catch
    % We start with an initial sub-division
    oldSubDivisionData = [];
    D = linspace(0,1,N);
end

% Record final sub-division
Dnew = [];
I = 0;
Stabledivision = zeros(1,length(D)-1);
for i = 1:length(D)-1
    [Ii,n,Di] = adaptiveGK715(phi_,dphi_,exp_ith,D(i),D(i+1),Points,0,dg_x,part);
    if length(Di) == 2
        Stabledivision(i) = 1;
    end
    Dnew = [Dnew Di];
    I = I+Ii;
end
cache(functionHash) = Dnew;

function [I,n,D] = adaptiveGK715(phi_,dphi_,exp_ith,a,b,GK715,n,dg_x,part)
% We have sub-divided down to a->b. Map standard Gauss-Kronrad points to
% this interval and then map those to original domain 0->inf
center = (a+b)/2;
halfwidth = (b-a)/2;
s = GK715.Nodes*halfwidth + center;
% and now map to original domain using nonlinear transformation
t = (s./(1-s)).^2;
coordinateChangeJacobian = (2*s./(1-s).^3).*halfwidth;

% Counter to diagnose depth of recursion
n = n+1;

% Evaluate all functions in all points
PHI    = phi_(t);
dPHI   = dphi_(t);
PHIZ   = prod(PHI,1);
RelPHI = dPHI./PHI;
% Careful edge-case when phi is 0
RelPHI(isnan(RelPHI))=0;

% Product characteristic, exponential and transformation can be combined
Weight = (PHIZ.*exp_ith(t).*coordinateChangeJacobian);
% Integral = sum of weighted columns
I_15 = (RelPHI*(Weight.*GK715.Weights_15).');
I_7 = (RelPHI(:,2:2:end)*(Weight(2:2:end).*GK715.Weights_7).');

% Simple stopping criteria for now, and we sub-divide in all despite some
% might be done already or regions being all 0. Must be optimized later
error = sum(abs(part(I_15-I_7)));
if all(error <= max(1e-6,(b-a)*1e-6))
    I = part(I_15);
    % Diagnostics on final interval
    D = [a b];
else
    [I1,n1,D1] = adaptiveGK715(phi_,dphi_,exp_ith,a,(a+b)/2,GK715,n,dg_x,part);
    [I2,n2,D2] = adaptiveGK715(phi_,dphi_,exp_ith,(a+b)/2,b,GK715,n,dg_x,part);
    I = I1+I2;
    % Deepest recursion
    n = max(n1,n2);
    % Diagnostics
    D = [D1 D2];
end

