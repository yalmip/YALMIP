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
        phi = varargin{3};  
        
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
        funcs = varargin{4};
        phi = varargin{5};
        % Create a function which computes gradient at x
        %operator.derivative = @(x)compute_dcdf_using_phi(x,funcs.h,funcs.dh,funcs.g,funcs.dg,phi);
        
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



function dcdf = compute_dcdf_using_phi(x,h,dh,g,dg, phi)
h = h(x);
dh = dh(x);
g = g(x);
dg = dg(x);

% This is wrong of course. fmincon will not converge, and if 
% 'fmincon.derivativecheck','on' is turned on in sdpsettings, it should
% terminate
dcdf = g;



