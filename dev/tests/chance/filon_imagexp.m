function val = filon_imagexp(f, a, b, h0)

% sampling
N = 2000;
x_ = linspace(a,b,N+1);
h_ = (b-a)/N;
fx = arrayfun(f, x_);
fx_real = real(fx);
fx_imag = imag(fx);
sin_h0_x = sin(h0*x_);
cos_h0_x = cos(h0*x_);

% compute integrals of cos and sin separately
partA = filon_sin_linear(@(xx) fx_real, x_, sin_h0_x);
partB = filon_cos_linear(@(xx) fx_imag, x_, cos_h0_x);
val = partA + partB;



function val = filon_sin_linear(f_handle, x, sin_h0_x)
N = length(x)-1;
val = 0;
for k=1:N
    x_ = x(k);
    x1 = x(k+1);
    f0 = f_handle(x_);
    f1 = f_handle(x1);
    s0 = sin_h0_x(k);
    s1 = sin_h0_x(k+1);

    val = val + filon_sin_segment_linear(x_, x1, f0, f1, s0, s1);
end



function seg_val = filon_sin_segment_linear(x0,x1,f0,f1,s0,s1)

h_  = x1 - x0;
a1 = (f1 - f0)/h_;      % slope

y0 = f0 * s0;
y1 = f1 * s1;
seg_val = 0.5*(y0 + y1)*h_;



function val = filon_cos_linear(f_handle, x, cos_h0_x)

N = length(x)-1;
val = 0;
for k=1:N
    x_ = x(k);
    x1 = x(k+1);
    f0 = f_handle(x_);
    f1 = f_handle(x1);
    c0 = cos_h0_x(k);
    c1 = cos_h0_x(k+1);

    val = val + filon_cos_segment_linear(x_,x1,f0,f1,c0,c1);
end


function seg_val = filon_cos_segment_linear(x0,x1,f0,f1,c0,c1)
h_ = x1 - x0;
y0 = f0*c0;
y1 = f1*c1;
seg_val = 0.5*(y0+y1)*h_; 

