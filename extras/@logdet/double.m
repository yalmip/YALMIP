function value = double(X);
%double           Overloaded

if isempty(X.cx)
    value = 0;
else
    value = double(X.cx);
end
for i = 1:length(X.P);
    value = value  + X.gain(i)*log(det(double(X.P{i})));
end
