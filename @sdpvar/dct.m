function X=dct(x,N)
%DCT (overloaded)

if min(size(x)) > 1
    error('DCT currently only supported for vectors');
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
w=[1/sqrt(N); sqrt(2/N)*ones(N-1,1)];
e = cos(pi*(k-1)'*(2*n-1)/(2*N));
X=w.*(e*x);