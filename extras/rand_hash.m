function r = rand_hash(k,n,m);

% Calling rng is pretty slow, to speed up massive amount of scalar
% definitions of sdpvars, we preallocate random number in batches of at
% least 10000 numbers at a time. We use rng since it is safer in parallell
% architechtures according to earlier user issues
persistent PredefinedRands
persistent PredefinedPointer
try
    if isempty(PredefinedPointer)
        PredefinedPointer = 1;
        s = rng;
        rng(k,'twister');
        PredefinedRands = rand(max(10000,2*n*m),1);
        rng(s);
    end
    if n*m > length(PredefinedRands)-PredefinedPointer+1
        PredefinedPointer = 1;
        s = rng;
        rng(k,'twister');
        PredefinedRands = rand(max(10000,2*n*m),1);
        rng(s);
    end
    r = reshape(PredefinedRands(PredefinedPointer:PredefinedPointer+n*m-1),n,m);
    PredefinedPointer = PredefinedPointer + n*m;
catch
    % Old versions of MATLAB
    s = rand('state');
    rand('state',k)
    r = rand(n,m);
    rand('state',s);
end