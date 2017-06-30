function [p,feasible,lower] = lpbmitighten(p,lower,upper,lpsolver,xmin,improvethese)

if nargin<6
    improvethese = ones(length(p.lb),1);
end

% Don't use LP to propagate variables which only enter as x + other in one
% single (in)equality. 
ind = sum(p.F_struc(1:p.K.f+p.K.l,:) | p.F_struc(1:p.K.f+p.K.l,:),1);
oneterm = find(ind(2:end) == 1);
for i = 1:length(oneterm)
    a = p.F_struc(1:p.K.f+p.K.l,oneterm(i)+1);
    j = find(a);
    if j <= p.K.f
        b = p.F_struc(j,:);
        if nnz(b(2:end))==2
            improvethese(oneterm(i))=0;
        end
    end
end

% Construct problem with only linear terms
% and add cuts from lower/ upper bounds

p_test = p;
p_test.F_struc = p.lpcuts;
p_test.K.l = size(p.lpcuts,1);
p_test.K.f = 0;
p_test.K.s = 0;
p_test.K.q = 0;
if ~isnan(lower) & ~isinf(lower)
    p_test.F_struc = [-(p.lower-abs(p.lower)*0.01)+p.f p_test.c';p_test.F_struc];
    if p.diagonalized
        n = length(p.c)/2;
        f = p.f;
        c = p.c(1:n);
        d = p.c(n+1:end);
        neg = find(d<0);
        if length(neg)>0
            f = f + sum(d(neg).*xmin(neg).^2 - 2*d(neg).*xmin(neg).^2);
            c(neg) = c(neg) + 2*d(neg).*xmin(neg);
            d(neg) = 0;
            p_test.F_struc = [-(p.lower-abs(p.lower)*0.01)+f c' d';p_test.F_struc];
            p_test.K.l = p_test.K.l + 1;
        end
    end
end
if upper < inf & ~(nnz(p.c)==0 &  nnz(p.Q)==0)
    if p.diagonalized
        n = length(p.c)/2;
        f = p.f;
        c = p.c(1:n);
        d = p.c(n+1:end);
        pos = find(d>0);
        if length(pos)>0
            f = f + sum(d(pos).*xmin(pos).^2 - 2*d(pos).*xmin(pos).^2);
            c(pos) = c(pos) + 2*d(pos).*xmin(pos);
            d(pos) = 0;
            p_test.F_struc = [upper+abs(upper)*0.01-f -c' -d';p_test.F_struc];
            p_test.F_struc = [upper+abs(upper)*0.01-p.f -p_test.c';p_test.F_struc];
            p_test.K.l = p_test.K.l + 2;
        end
    end
    p_test.F_struc = [upper+abs(upper)*0.01-p.f -p_test.c';p_test.F_struc];
    p_test.K.l = p_test.K.l + 1;
end


if p.options.bmibnb.cut.evalvariable
    p_test = addBilinearVariableCuts(p_test);
end
if p.options.bmibnb.cut.evalvariable
    p_test = addEvalVariableCuts(p_test);
end
if p.options.bmibnb.cut.multipliedequality
    p_test = addMultipliedEqualityCuts(p_test);
end

% Try to get rid of numerical noise (this is far from stringent, but it
% works, and help GLPK from crashing in some instances)
%p_test.F_struc = unique(round(p_test.F_struc*1e12)/1e12,'rows');
p_test.K.l = size(p_test.F_struc,1);
p_test.F_struc = [p.F_struc(1:1:p.K.f,:);p_test.F_struc];
p_test.K.f = p.K.f;

feasible = 1;

% Perform reduction on most violating approximations at current solution
if p.options.bmibnb.lpreduce ~= 1
    n = ceil(max(p.options.bmibnb.lpreduce*length(p_test.linears),1));
    res = zeros(length(p.lb),1);
    for i = 1:size(p.bilinears,1)
        res(p.bilinears(i,2)) = res(p.bilinears(i,2)) + abs( p.x0(p.bilinears(i,1))-p.x0(p.bilinears(i,2)).*p.x0(p.bilinears(i,3)));
        res(p.bilinears(i,3)) = res(p.bilinears(i,3)) + abs( p.x0(p.bilinears(i,1))-p.x0(p.bilinears(i,2)).*p.x0(p.bilinears(i,3)));
    end
    res = res(p.linears);
    [ii,jj] = sort(abs(res));
    jj = jj(end-n+1:end);
else
    jj=1:length(p_test.linears);
end

j = 1;
while feasible & j<=length(jj)
    i = p_test.linears(jj(j));
    if abs(p.ub(i)-p.lb(i)>p.options.bmibnb.vartol) & improvethese(i)    
        p_test.c = eyev(length(p_test.c),i);        
        output = feval(lpsolver,removenonlinearity(p_test));
        p.counter.lpsolved = p.counter.lpsolved + 1;
        if output.problem == 0 | output.problem == 2 | output.problem == 12
            if output.problem == 0
                if p_test.lb(i) < output.Primal(i)-1e-5
                    p_test.lb(i) = output.Primal(i);
                    p_test = updateonenonlinearbound(p_test,i);
                end
            end
            p_test.c = -p_test.c;          
            output = feval(lpsolver,removenonlinearity(p_test));   
            p.counter.lpsolved = p.counter.lpsolved + 1;
            if output.problem == 0
                if p_test.ub(i) > output.Primal(i)+1e-5
                    p_test.ub(i) = output.Primal(i);
                    p_test = updateonenonlinearbound(p_test,i);
                end
            end
            if output.problem == 1
                feasible = 0;
            end
        end
        if output.problem == 1
            feasible = 0;
        end
    end
    j = j + 1;
end
p_test = clean_bounds(p_test);
p.lb = p_test.lb;
p.ub = p_test.ub;
% Save these for re-use
p.evalMap = p_test.evalMap;