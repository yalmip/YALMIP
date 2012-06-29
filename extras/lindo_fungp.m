function [f,err] = fmincon_funn(cbData,nRow,x,njdiff,dXjbase,reserved,inParam)

persistent prob Af A bf b

if nargin == 7    
    prob = inParam;
    ind = find(prob.map==0);
    Af = prob.A(ind,:);
    bf = prob.b(ind);
    A = spalloc(max(prob.map),size(prob.A,1),0);
    for i = 1:max(prob.map)
        ind = find(prob.map==i);
        A(i,ind) = prob.b(ind)';
    end
    return
end

x = x(1:size(prob.A,2));
if nRow == -1    
    z = Af*x;
    f = full(log(bf'*exp(z)));    
else  
    g = A*exp(prob.A*x);
    ind = find(g<1e-2);
    if ~isempty(ind)
        g(ind) = exp(log(1e-2)+(g(ind)-1e-2)/1e-2);
    end
    g = log(g);

    if length(prob.h) > 0
        g = [log(prob.h) + prob.G*x;g];
    end

    f = g(nRow+1);
end
err = 0;