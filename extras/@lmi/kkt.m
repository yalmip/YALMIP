function [KKTConstraints, details] = kkt(F,h,parametricVariables,ops);
%KKT Create KKT system
%
% [KKTConstraints, details] = kkt(Constraints,Objective,parameters,options)

if ~isempty(F)
    if any(is(F,'sos2'))
        error('SOS2 structures not allowed in KKT');
    end
end

[aux1,aux2,aux3,model] = export(F,h,sdpsettings('solver','quadprog','relax',2));
if isempty(model)
    error('KKT system can only be derived for LPs or QPs');
end
model.problemclass.constraint.binary = 0;
model.problemclass.constraint.integer = 0;
model.problemclass.constraint.semicont = 0;
model.problemclass.constraint.sos1 = 0;
model.problemclass.constraint.sos2 = 0;
if ~ismember(problemclass(model), {'LP', 'Convex QP', 'Nonconvex QP'})
    error('KKT system can only be derived for LPs or QPs');
end
if ~isempty(model.binary_variables) | ~isempty(model.integer_variables) | ~isempty(model.semicont_variables)
    comb = [model.used_variables(model.binary_variables) model.used_variables(model.integer_variables) model.used_variables(model.semicont_variables) ];
    if ~isempty(intersect(comb,getvariables(decisionvariables)))
        error('KKT system cannot be derived due to binary or integer decision variables');
    end
end

if nargin < 3
    x = recover(model.used_variables);
    parameters = [];
else
    % Make sure they are sorted    
    parameters = getvariables(parametricVariables);
    x = recover(setdiff(model.used_variables,parameters));
    notparameters = find(ismember(model.used_variables,getvariables(x)));%parameters);
    parameters = setdiff(1:length(model.used_variables),notparameters);
    if ~isempty(parameters)
        y = recover(model.used_variables(parameters));
    end
end

if nargin < 4
    ops = sdpsettings;
end
if isempty(ops)
    ops = sdpsettings;
end

% Ex==f
E = model.F_struc(1:model.K.f,2:end);
f = -model.F_struc(1:model.K.f,1);
% Ax <= b
A = -model.F_struc(model.K.f+1:model.K.f+model.K.l,2:end);
b = model.F_struc(model.K.f+1:model.K.f+model.K.l,1);
c = model.c;
Q = model.Q;

% Are there equalities hidden in the in equalities
moved = zeros(length(b),1);
hash = randn(size(A,2),1);
Ahash = A*hash;

removed = zeros(length(b),1);
for i = 1:length(b)
    if ~removed(i)
        same = setdiff(find(Ahash == Ahash(i)),i);
        k = find(b(i) == b(same));
        if ~isempty(k)
            removed(same(k)) = 1;
        end
        k = find(b(i) < b(same));
        if ~isempty(k)
            removed(same(k)) = 1;
        end
    end
end
if any(removed)
    if ops.verbose
        disp(['Removed ' num2str(nnz(removed)) ' duplicated or redundant inequalities']);
    end
    A(find(removed),:) = [];
    b(find(removed)) = [];
    Ahash = A*hash;
end

for i = 1:length(b)
    if ~moved(i)
        same = setdiff(find(Ahash == -Ahash(i)),i);
        k = find(b(i) == -b(same));
        if ~isempty(k)
            E = [E;A(i,:)];
            f = [f;b(i)];
            moved(i) = 1;
            moved(same(k)) = 1;
        end
    end
end

if any(moved)
    if ops.verbose
        disp(['Transfered ' num2str(nnz(moved)) ' inequalities to equalities']);
    end
    A(find(moved),:) = [];
    b(find(moved)) = [];
    Ahash = A*hash;
end


if ~isempty(b)
    infBounds = find(b>0 & isinf(b));
    if ~isempty(infBounds)
        b(infBounds) = [];
        A(infBounds,:) = [];
    end
end

[dummy,rr] = unique([A b],'rows');
if length(rr)~=size(A,1)
    A = A(rr,:);
    b = b(rr);
end

if ~isempty(parameters)
    b = b-A(:,parameters)*y;
    f = f-E(:,parameters)*y;
    A(:,parameters) = [];
    E(:,parameters) = [];
    c(parameters) = [];
    Q2 = model.Q(notparameters,parameters);
    c = c + 2*Q2*y;
    Q = Q(notparameters,notparameters);
end

used = find(any(A,2));
if isempty(setdiff(1:size(A,1),used))
    parametricDomain = [];
else
    r = setdiff(1:size(A,1),used);
    parametricDomain = b(r)>=0;
    A(r,:)=[];
    b(r)=[];
end

if isempty(E)
    E = [];
    f = [];
end

any_inequalities = size(A,1)>0;
any_equalities = size(E,1)>0;
Lambda = sdpvar(size(A,1),1); % Ax  <=b
mu     = sdpvar(size(E,1),1); % Ex+f==0

if ops.kkt.dualbounds
    if ops.verbose
        disp('Starting derivation of dual bounds (can be turned off using option kkt.dualbounds)');
    end
    [U,pL] = derivedualBounds(2*Q,c,A,b,E,f,ops,parametricDomain);
    %[U,pL] = derivedualBoundsParameterFree(2*Q,c,A,b,E,f,ops,parametricDomain);
else
    U = [];
end
KKTConstraints = [];
s = 2*Q*x + c;
if any_inequalities
    %KKTConstraints = [KKTConstraints, complements(b-A*x, Lambda >= 0):'Compl. slackness and primal-dual inequalities'];
    KKTConstraints = [KKTConstraints, complements(Lambda >= 0,b-A*x):'Compl. slackness and primal-dual inequalities'];
    s = s + A'*Lambda;
end
if any_equalities
    KKTConstraints = [KKTConstraints, (E*x == f):'Primal feasible'];
    s = s + E'*mu;
end
if ~isempty(U)
    finU = find(~isinf(U));
    if ~isempty(finU)
        KKTConstraints = [KKTConstraints, (Lambda(finU) <= U(finU)):'Upper bound on duals'];
    end
end
KKTConstraints = [KKTConstraints, (s == 0):'Stationarity'];

if any_inequalities
    s_ = indicators(KKTConstraints(1));
    if ops.kkt.licqcut
        if ops.verbose
            disp('Generating LICQ cuts');
        end
        [Alicq,blicq] = createLICQCut(A);
        KKTConstraints = [KKTConstraints, Alicq*s_ <= blicq];
    end
    
    
    if ops.kkt.minnormdual;
        MinNorm = [0 <= Lambda <= 10000*s_,s == 0];
        ops.kkt.minnormdual = ops.kkt.minnormdual-1;
        parametricInMinNorm = recover(setdiff(depends(MinNorm),depends(Lambda)));
        [kkt2,info2] = kkt(MinNorm,Lambda'*Lambda,parametricInMinNorm,ops);
        kkt2 = [kkt2, [indicators(kkt2(1)) <= 1-[s_;s_]]];
        kkt2 = [kkt2, info2.dual <= 10000];
        KKTConstraints = [KKTConstraints, kkt2];
        details.info2 = info2;
    end
end

if nnz(Q)>0
    details.info = 'min x''Qx+c''x';
else
    details.info = 'min c''x';
end
if isempty(E)
    details.info = [details.info ' s.t. Ax<=b'];
else
    details.info = [details.info ' s.t. Ax<b, Ex=f'];
end
details.c = c;
details.Q = Q;
details.A = A;
details.b = b;
details.E = E;
details.f = f;
details.dual = Lambda;
details.dualeq = mu;
details.primal = x;
if length(b)>0
    details.inequalities = A*x <= b;
else
    details.inequalities = [];
end
if length(f)>0
    details.equalities = E*x == f;
else
    details.equalities = [];
end
if isempty(U)
    U = repmat(inf,length(b),1);
end
details.dualbounds = U;


function [Alicq,blicq] = createLICQCut(A);
Alicq = ones(1,size(A,1));
blicq = size(A,2);

Atemp = A;
[ii,jj] = sort(sum(A | A,1));
Atemp = Atemp(:,jj);
[ii,jj] = sort(sum(A | A,2));
Atemp = Atemp(jj,:);
for ii = 1:size(Atemp,1)
    rr(ii) = min(find(Atemp(ii,:)));
end
[ii,kk] = sort(rr);
kk = fliplr(kk);
Atemp = Atemp(kk,:);

nel = sum(Atemp|Atemp,2);
top = 1;
while top <= length(nel)
    same = nel == nel(top);
    used_in_first = find(Atemp(top,:));
    oldtop = top;
    while top <= length(nel)-1 && all(ismember(find(Atemp(top+1,:)),used_in_first))
        top = top + 1;
    end
    if length(used_in_first)<size(A,2)
        blicq = [blicq;nnz(used_in_first)];
        e = spalloc(1,size(A,1),0);
        e(jj(kk(oldtop:top)))=1;
        Alicq = [Alicq;e];            
    end
    if  nnz(any(full(Atemp(1:top,:)),1)) < size(A,2)
        blicq = [blicq; nnz(any(full(Atemp(1:top,:)),1))];
        e = spalloc(1,size(A,1),0);
        e(jj(kk(1:top))) = 1;
        Alicq = [Alicq;e];        
    end
    rr=oldtop:top;
    ss = find(nel(rr)==1);
    if ~isempty(ss)
        blicq = [blicq;1];
        e = spalloc(1,size(A,1),0);
        e(jj(kk(rr(ss)))) = 1;
        Alicq = [Alicq;e];       
    end
    top = top + 1;
end




