function varargout = cdf_from_characteristic(varargin)

switch class(varargin{1})
    
    case 'double'
        x = varargin{1};
        phi = varargin{2};
        varargout{1} = compute_cdf(x,phi);
        
    case 'sdpvar'
        varargout{1} = InstantiateElementWise('cdf_from_characteristic',varargin{:});
        
    case 'char'
        
        operator = struct('convexity','none','monotonicity','increasing','definiteness','positive','model','callback');
        operator.bounds = @bounds;
        operator.range = [0 1];
        operator.derivative = @(x)compute_pdf(x,varargin{4:end});      
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
end

function cdf = compute_cdf(x,phi)
% Perform the inverse Fourier transform to obtain the CDF
cdf = zeros(size(x));
for k = 1:length(x)
    integrand = @(t) imag(phi(t) .* exp(-1i * t * x(k))./t);
    cdf(k) = .5-integral(integrand,0, 100)/pi;
end

function pdf = compute_pdf(x,phi)
% Perform the inverse Fourier transform to obtain the PDF
pdf = zeros(size(x));
for k = 1:length(x)
    integrand = @(t) phi(t) .* exp(-1i * t * x(k));
    pdf(k) = real(integral(integrand,0, 100) / (pi));
end
