function f = int_sdpvar(f,x,from,to)

for i = 1:length(x)
    [c,v] = coefficients(f,x(i));
    vnew = [];
    for j = 1:length(v)
        di = degree(v(j),x(i));
        vnew = [vnew;to(i)^(di+1)/(di+1)-from(i)^(di+1)/(di+1)];
    end
    f = c'*vnew;
end