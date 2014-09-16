function r = gen_rand_hash(k,n,m);
persistent buffer
persistent loc
try
    % Preallocate 100 samples everytime and empty that buffer first
    if isempty(buffer) || length(buffer) < loc + n*m
        % Previous approach does not work with parallell toolbox!
        s = rng;
        rng(k);
        buffer = rand(n*m+100,1);
        r = reshape(buffer(1:n*m),n,m);
        loc = 1;
        rng(s);
    end
     r = reshape(buffer(loc:n*m+loc-1),n,m);
     loc = loc + n*m;
    
catch
    % but rng is not available in all versions...
    s = rand('state');
    rand('state',k)
    r = rand(n,m);
    rand('state',s);
end