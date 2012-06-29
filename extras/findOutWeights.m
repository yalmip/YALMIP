function w = findOutWeights(X,weights)

if length(weights)>1
    if length(weights)~=length(X)
        error('The weights vector should have same length as the decision variable');
    end
end

if isa(weights,'double')
    w = weights(:);
else
    for i = 1:length(X)
        x = X(i);
        j = getvariables(x);
        W = getbasematrix(weights,j);
        switch nnz(W)
            case 0
                w(i) = 1;
            case 1
                w(i) = W(find(W));
            otherwise
                error('The weights vector is incorrect (same variable in several locations)');
        end
    end
end