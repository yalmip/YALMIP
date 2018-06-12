function X=fft(x,N)
%FFT (overloaded)

if min(size(x)) > 1
    error('FFT currently only supported for vectors');
end

% Make sure x is column vector
x = reshape(x,[],1);

% Pad?
if nargin == 1    
    N = length(x);
else    
    x = [x;zeros(N-length(x),1)];
end

% Transform
k = 1:N;
n = 1:N;
e = exp(-j*2*pi*(k-1)'*(n-1)/N);
X=e*x;