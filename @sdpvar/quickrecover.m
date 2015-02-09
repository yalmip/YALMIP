function x = quickrecover(x,v,b);

x.lmi_variables = v;
x.conicinfo = [0 0];
if nargin == 3
    x.basis(2) = b;
end