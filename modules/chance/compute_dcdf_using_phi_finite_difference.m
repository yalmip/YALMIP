function dcdf = compute_dcdf_using_phi_finite_difference(x,h,dh,g,dg, phi,dphi)
% For internal testingedit ch
h0 = h(x);
dh0 = dh(x);
g0 = g(x);
dg0 = dg(x);

% Silly numerical differentiation (which is what a solver does anyway)
phi0 = @(t) prod(phi(g0(:)*t),1);  
cdf0 = compute_cdf_using_phi(-h0,phi0);
eps = 1e-12;
dcdf = [];
for k = 1:length(x)
    x_ = x;x_(k) = x_(k)+eps;
    g_ = g(x_);
    h_ = h(x_);    
    phi_ = @(t) prod(phi(g_(:)*t),1);  
    cdf_ = compute_cdf_using_phi(-h_,phi_);
    dcdf = [dcdf;(cdf_-cdf0)/eps];
end