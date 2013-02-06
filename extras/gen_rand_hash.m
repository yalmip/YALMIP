function r = gen_rand_hash(k,n,m);
try
    % Previous approach does not work with parallell toolbox!
    s = rng;
    rng(k);
    r = rand(n,m);
    rng(s);
catch
    % but rng is not available in all versions...
    s = rand('state');
    rand('state',k)
    r = rand(n,m);
    rand('state',s);
end