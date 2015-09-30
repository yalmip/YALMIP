function f = int_sdpvar(f,x,from,to)

if length(f)>1
    dim = size(f);
    f = f(:);
    F = [];
    for i = 1:length(f)
        F = [F;int_sdpvar(f(i),x,from,to)];
    end
    f = reshape(F,dim);
    return
end

for i = 1:length(x)
    if ~isa(x(i),'sdpvar')
        error('An element in the integration variable is not an sdpvar. Third argument must be an sdpvar');
    end
    [c,v] = coefficients(f,x(i));
    vnew = [];
    for j = 1:length(v)
        di = degree(v(j),x(i));
        vnew = [vnew;to(i)^(di+1)/(di+1)-from(i)^(di+1)/(di+1)];
    end
    f = c'*vnew;
end