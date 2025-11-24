function I = integral_lazy_rotation(phi, w)
% Computes Integral_0^inf [phi_z(t) * exp(i*w*t)] dt
% Uses the pure Contour Deformation (Lazy Rotation)

theta = 30 * (pi/180)*sign(w);
line_dir = exp(1i*theta);
integrand_complex = @(r) (phi(r*line_dir) .* exp(1i * w * (r*line_dir))) .* line_dir;
I = integral(integrand_complex, 0, Inf, 'RelTol', 1e-12);

if 0
    phi = @(t)1./(1+t.^2);
    t = 1e-5:0.0001:1000;
    w = 30;
    plot(t,real(phi(t).*exp(w*1i*t)))
    sum(phi(t).*exp(w*1i*t))*mean(diff(t))
    integral_lazy_rotation(phi,w)
    
    theta = 30 * (pi/180)*sign(w);
    line_dir = exp(1i*theta);
    integrand_complex = @(r) (phi(r*line_dir) .* exp(1i * w * (r*line_dir))) .* line_dir;
    r = 1e-5:0.0001:1000;
    plot(r,real(integrand_complex(r)))
end