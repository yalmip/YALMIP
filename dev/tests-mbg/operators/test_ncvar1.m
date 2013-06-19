function test_ncvar1

yalmip('clear')
sdpvar x y z
ncvar a b c
% Test of simple algebra
mbg_assertequal(a*b-a*b,0)
mbg_assertequal(a*x*c*y-a*x*c*y,0)
mbg_assertfalse(isequal(a*c-c*a,0))
mbg_assertequal(b*b-b^2,0)
mbg_assertequal(a*b*b-a*b.^2,0)

mbg_assertequal(b^2-b*b,0)
mbg_assertequal(b^2-b.*b,0)
mbg_assertequal(b^2-b.^2,0)
mbg_assertequal((a*b)^2-a*b*a*b,0)
mbg_assertequal((a*b).^2-a*b*a*b,0)
mbg_assertequal((a*b).^2-(a*b).*(a*b),0)

mbg_assertequal((a*b*x).^2-x^2*(a*b).*(a*b),0)

% TEST COEFFICIENTS
yalmip('clear')
sdpvar x y z
ncvar a b c
p = a*x + (b+c)*x^2 + x*y*z;
[c,v] = coefficients(p,[x;y;z])

yalmip('clear')
sdpvar x y z
ncvar a b c
[coeff,v] = coefficients(a*b+(a*b*b*c)*x,x);
mbg_assertequal(coeff(2)-a*b*b*c,0)

[coeff,v] = coefficients(a*b*x,y);
mbg_assertequal(coeff-a*b*x,0)

% TEST REPLACE
mbg_assertequal(replace(a*b,x,y)-a*b,0)
mbg_assertequal(replace(a*b,b,c)-a*c,0) 
mbg_assertequal(replace(a*b,a,b)-b^2,0) 
mbg_assertequal(replace(a*b,a,a+b+c)-a*b-b*b-c*b,0) 
mbg_assertequal(replace(a*b*c,b,x^2)-x^2*a*c,0) 
mbg_assertequal(replace(a*b*c,b,3)-3*a*c,0) 
mbg_assertequal(replace(a*b*c,[a b c],[1 2 a])-2*a,0) 
mbg_assertequal(replace(a*b*x,x,y)-a*b*y,0) 
mbg_assertequal(replace(a*b*x,x,c)-a*b*c,0) 
mbg_assertequal(replace(a*b*x^2,x,c)-a*b*c^2,0)
mbg_assertequal(replace(a*b*x^2,y,c)-a*b*x^2,0)



