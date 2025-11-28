function I = integral_lazy_rotation(phi, w,theta,part)
% Computes Integral_0^inf [phi_z(t) * exp(i*w*t)] dt
% Uses the pure Contour Deformation (Lazy Rotation)

if nargin < 3 || isempty(theta)
    theta = 30*pi/180;
end

if nargin < 4 
    % Support send real/imag operator
    part = @(x)x;
end

if w ~= 0
    theta = theta*sign(w);
    z = exp(1i*theta);
    integrand_complex = @(r) part((phi(r*z) .* exp(1i * w * (r*z))) .* z);
    I = integral(integrand_complex, 0, Inf, 'RelTol', 1e-12);
else
    % Case w = 0 should have been caught up-stream but ok
    I = integral(phi, 0, Inf, 'RelTol', 1e-12);
end