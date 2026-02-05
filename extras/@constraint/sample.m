function S = sample(C,N)

if nargin < 2
    N = 1;
end

if length(C)>1
    error('FIXME: Can only sample indivual constraint for now')
else
    x = sdpvar(C);
    X = sample(x,N);
    S = [];
    for i = 1:N
        S = [S, (X(:,i) >= 0): ['Random sample ' num2str(i)] ];
    end
end