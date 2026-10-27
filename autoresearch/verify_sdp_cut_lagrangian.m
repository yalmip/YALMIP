function verify_sdp_cut_lagrangian
% Exhaustively check the proposed cut bound on small binary models.
rng(1729);
checked = 0;
for trial = 1:300
    n = 5 + mod(trial,4);
    p.binary_variables = (1:n)';
    objectiveScale = 10^randi([-4 4]);
    cutScale = 10^randi([-4 4]);
    p.c = objectiveScale*2*randn(n,1);
    p.f = objectiveScale*randn;
    p.lb = zeros(n,1);
    p.ub = ones(n,1);
    fixed = rand(n,1) < 0.3;
    values = double(rand(n,1) < 0.5);
    p.lb(fixed) = values(fixed);
    p.ub(fixed) = values(fixed);
    p.binarycardinality.down = randi([0 n]);
    p.binarycardinality.up = randi([p.binarycardinality.down n]);

    a = cutScale*randn(n,1);
    b = cutScale*randn;
    sdpCuts.K.l = 1;
    sdpCuts.F_struc = [b a'];
    bound = sdpCutLagrangianBound(p,sdpCuts,1,inf);

    patterns = dec2bin(0:2^n-1)-'0';
    eligible = all(patterns(:,fixed) == values(fixed)',2) & ...
        sum(patterns,2) >= p.binarycardinality.down & ...
        sum(patterns,2) <= p.binarycardinality.up & ...
        b + patterns*a >= -1e-12;
    if any(eligible)
        costs = p.f + patterns(eligible,:)*p.c;
        exact = min(costs);
        assert(bound <= exact + 1e-7*(1+abs(exact)), ...
            'SDP cut bound exceeded a feasible binary objective.');
        checked = checked + 1;
    end
end
fprintf('Verified %d feasible randomized binary-cut instances.\n',checked);
end
