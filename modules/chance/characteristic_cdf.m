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
        reldphi = distribution.characteristicfunction_relativederivative;
        % Create a function which computes gradient at x
        operator.derivative = @(x)compute_dcdf_using_phi(x,funcs.h,funcs.dh,funcs.g,funcs.dg,phi,dphi,reldphi);
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
end

function cdf = compute_cdf_using_phi(y,phi)
% Perform the inverse Fourier transform to obtain the CDF
% for Probability(z <= y) where z has characterstic function phi(t)
% Not much to exploit so for now we use built-in integrator
integrand = @(t) imag(phi(t) .* exp(-1i * t * y)./t);
cdf =  .5-integral(integrand,0, inf)/pi;


function dcdf = compute_dcdf_using_phi(x,h,dh,g,dg,phi,dphi,reldphi)
% Compute derivative of Probability(h(x)+g(x)'*w <= 0)

% Evaluate the operators
h0  = h(x);
dh0 = dh(x);
g0  = g(x);
dg0 = dg(x);

% Instantiate characteristic functions
phi_z   = @(t) prod(phi(g0(:)*t),1);
phi_    = @(t) phi(g0(:)*t);
dphi_   = @(t) dphi(g0(:)*t);
exp_ith = @(t) exp(1i*t*h0);

% Compute f_z(-h0) using Gil-P
% this will be moved to be computed together with derivative instead but
% for now we do the double work to keep code simple
integrand_pdf = @(t) real(exp_ith(t) .* phi_z(t));
if nnz(dh0)==0
    pdf_val = 0;
else
    pdf_val = (1/pi)*integral(integrand_pdf,0,inf);
end

if nnz(dg0)==0
    dcdf = (-pdf_val.*dh0');
else            
    terms = (-1/pi)*dcdfintegralEvaluator(phi_,dphi_,exp_ith);
    dcdf = (-pdf_val.*dh0') + terms'*dg0;
end

function I = dcdfintegralEvaluator(phi_,dphi_,exp_ith)
% Compute integral. The subdivided regions are not used yet but might be
% used later to speed up repeated calls in some way
[I,regions] = adaptiveGK715(phi_,dphi_,exp_ith,0,1,setupGK715);

function [I,regions] = adaptiveGK715(phi_,dphi_,exp_ith,a,b,GK715)
% Map nodes to [a,b]
midpt = (a+b)/2;
halfh = (b-a)/2;
s = GK715.Nodes*halfh + midpt; 
% and now map to original domain using nonlinear transformation
t = (s./(1-s)).^2;
coordinateChangeJacobian = (2*s./(1-s).^3).*halfh;

% Evaluate all functions in all points
PHI    = phi_(t);
dPHI   = dphi_(t);
PHIZ   = prod(PHI,1);
RelPHI = dPHI./PHI;
% Careful edge-case
RelPHI(isnan(RelPHI))=0;

% Weights, exponential and transformation can be moved into weights
temp = (PHIZ.*exp_ith(t).*coordinateChangeJacobian);
I_15 = imag(RelPHI*(temp.*GK715.Weights_15).');
I_7 = imag(RelPHI(:,2:2:end)*(temp(2:2:end).*GK715.Weights_7(2:2:end)).');
error = abs(I_15-I_7);
if all(error <= 1e-12)
    I = I_15;
    regions = {[a b]};
else
    [I1,regions1] = adaptiveGK715(phi_,dphi_,exp_ith,a,(a+b)/2,GK715);
    [I2,regions2] = adaptiveGK715(phi_,dphi_,exp_ith,(a+b)/2,b,GK715);
    I = I1+I2;
    regions = {regions1{:},regions2{:}};
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
    0, 0.1294849661688697, 0, ...
    0.2797053914892767, 0, 0.3818300505051189, ...
    0, 0.4179591836734694, 0, ...
    0.3818300505051189, 0, 0.2797053914892767, ...
    0, 0.1294849661688697, 0];


%GK715.Nodes = linspace(-1+1e-12,1-1e-12,15);
%GK715.Weights_15 = ones(1,15)/15;
%GK715.Weights_7 = GK715.Weights_15*2;
%GK715.Weights_7(1:2:end-1) = 0;