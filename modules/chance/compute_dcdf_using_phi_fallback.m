function dcdf = compute_dcdf_using_phi_fallback(x,h,dh,g,dg,phi,dphi,reldphi)
% Compute derivative of Probability(h(x)+g(x)'*w <= 0)
h0  = h(x);
dh0 = dh(x);
g0  = g(x);
dg0 = dg(x);

% Define phi_z(t) and weighted relative derivative
phi_z    = @(t) prod(phi(g0(:)*t),1);
reldphi_ = @(t)(reldphi(g0(:)*t));

% Define exponential weight
exp_ith = @(t) exp(1i*t*h0);

% Compute pdf f_z(-h0)
integrand_pdf = @(t) real(exp_ith(t) .* phi_z(t));
pdf_val = (1/pi)*integral(integrand_pdf,0,inf);

% Compute derivative
integrand_dcdf = @(t) imag(exp_ith(t) .* phi_z(t) .* rr(t));
if nnz(dg0)==0
	dcdf = (-pdf_val.*dh0');
else
    terms = (-1/pi)*integral(integrand_dcdf,0,inf,'ArrayValued',true);    
    dcdf = (-pdf_val.*dh0') + terms'*dg0;
end
% dcdf = compute_dcdf_using_phi_finite_difference(x,h,dh,g,dg, phi,dphi);
%[dcdf(:) dcdf_check(:)]

