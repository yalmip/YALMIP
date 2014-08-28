function [F,vars] = sort_internal(t,X,data);%ii,D,V);

% Hack to figure out all the sorted variables, not just this index.
% Sort is implemented in a slighlt different way (general feature in
% future versions) that allows one element in an operator to modell all
% elements. Reduces the number of calls to the operator code.

ii = data.i;
D = data.D;
V = data.V;

var_start = getvariables(t)-ii+1;
n = length(X);
var_end   = getvariables(t)-ii+1 + n-1;
vars = var_start:var_end;

% Is this a location variable instead of the actuial sort variable. If so,
% shift everything back to get the indicies of the sorted variables.
if data.isthisloc == 1
    vars = vars - n;
end

t   = recover(vars);
loc = recover(vars+n);

[M,m] = derivebounds(X);
X = reshape(X,1,n);

% Standard model
F = (sum(D,1) == 1) + (sum(D,2) == 1);
F = F + (t == sum(V,2));
F = F + (diff(t) >= 0);
for i = 1:n
   di = D(i,:);
   vi = V(i,:);
   F = F + (-(-m)'.*(1-di) <= X-vi <= (M)'.*(1-di));
   F = F + (m'.*di <=  vi <= M'.*di);   
end

% Cuts
F = F + (X == sum(V,1));
F = F + (sum(t) == sum(X));

% Definition of location
F = F + (loc == D*[(1:n)']) + (1 <= loc <= n);

