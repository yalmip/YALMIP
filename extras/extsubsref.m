function z=ihatethis(x,i,j);

if nargin < 3
    z = x(i);
else
    z=x(i,j);
end
