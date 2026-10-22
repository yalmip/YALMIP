function X=fft(x,N)
%FFT (overloaded)

if min(size(x)) > 1
    error('FFT currently only supported for vectors');
end

if nargin == 1 || isempty(N)
    N = length(x);
end

% Transform each coefficient column, retaining the input orientation.
X = x;
X.basis = sparse(fft(full(x.basis),N,1));
if size(x,1) == 1
    X.dim = [1 N];
else
    X.dim = [N 1];
end
X.conicinfo = [0 0];
X = clean(X);
