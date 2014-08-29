function   [F,obj,options] = create_lrmodel(BlockedA,Blockedb,F_parametric,parobj,options,ParametricVariables)
% Some special code for ther low-rank model in SDPLR
% Experimental code, not official yet
allb = [];
allA = [];
K.s = [];
for i = 1:length(Blockedb)
    allb = [allb;Blockedb{i}];
    Ai = [];
    for j = 1:size(BlockedA{i},2)
        Ai = [Ai BlockedA{i}{j}];
        K.s = [K.s sqrt(size(BlockedA{i}{j},2))];
    end
    %blkdiag bug in 7.0...
    [n1,m1] = size(allA);
    [n2,m2] = size(Ai);
    allA = [allA spalloc(n1,m2,0);spalloc(n2,m1,0) Ai];
end
options.solver = 'sdplr';
z = recover(ParametricVariables)
start = size(BlockedA,2)+1;
Mi = [];
for i = 1:length(allb)
    if isa(allb(i),'sdpvar')
        [Qi,ci,fi,xi,infoi] = quaddecomp(allb(i),z);
    else
        Qi = zeros(length(z));
        ci = zeros(length(z),1);
        fi = allb(i);
    end
    Z = -[fi ci'/2;ci/2 Qi];
    Mi = [Mi;Z(:)'];
end
K.s = [K.s length(z)+1];
zeroRow = zeros(1,size(allA,2));
allA = [allA Mi;zeroRow 1 zeros(1,K.s(end)^2-1)];
b = zeros(size(allA,1),1);b(end) = 1;
y = sdpvar(length(b),1);
CminusAy = -allA'*y;
start = 1;

% Get the cost, expressed in Z
[Qi,ci,fi,xi,infoi] = quaddecomp(parobj,z);
C = [fi ci'/2;ci/2 Qi];
F = ([]);
for i = 1:length(K.s)
    if i<length(K.s)
        F = F + (reshape(CminusAy(start:start+K.s(i)^2-1),K.s(i),K.s(i)) >= 0);
    else
        F = F + (reshape(C(:) + CminusAy(start:start+K.s(i)^2-1),K.s(i),K.s(i))>=0);
    end
    start = start + K.s(i)^2;
end
obj = -b'*y;

options.sdplr.forcerank = ceil(K.s/2);
options.sdplr.forcerank(end) = 1;
options.sdplr.feastol = 1e-7;
