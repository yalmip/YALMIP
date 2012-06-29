function p = pwaproj(h,w)

if prod(size(h))>1
    error('pwaproj is currently only applicable to scalars');
end

hvars = getvariables(h);
wvars = getvariables(w);
basis = getbase(h);
coefficients = basis(2:end);

Complicating = find(ismember(hvars,yalmip('extvariables')));
if length(Complicating)==length(hvars)
    LinearTerm = 0;
else
    NonComplicating = find(~ismember(hvars,yalmip('extvariables')));
    LinearwIndex = find(ismember(hvars,wvars));
    NonComplicating = setdiff(NonComplicating,LinearwIndex);
    LinearTerm = recover(hvars(NonComplicating))'*coefficients(NonComplicating)';
    if isempty(LinesarwIndex)
        Linearw = 0;
    else
        Linearw = recover(hvars(LinearwIndex))'*coefficients(LinearwIndex)';
    end
end

if isempty(Complicating)
    p = h;
    return
end

basis = getbase(h);
coefficients = basis(2:end);

totalobjective = 0;
totalF = [];
allz = [];
allargs = [];
if ~isempty(Complicating)
    for i = 1:length(Complicating)
        extstruct = yalmip('extstruct',hvars(Complicating(i)));
        [properties{i},F{i},arguments{i}]=model(recover(hvars(Complicating(i))));
        if intersect(depends(arguments{i}),wvars)
            RequiresEpi(i) = 1;
            z = sdpvar(size(arguments{i},1),size(arguments{i},2));
            allz = [allz;z];
            allargs = [allargs;arguments{i}];
            extstruct.arg{1} = z;
            t = feval(properties{i}{1}.name,z);
            f = eval(['@' properties{i}{1}.name]);
            t = f(extstruct.arg{1:end-1});
            [dummy{i},thisF]=model(t);
            totalF = totalF + thisF;
            totalobjective=coefficients(Complicating(i))*t+totalobjective;
        else
           LinearTerm = LinearTerm + recover(hvars(Complicating(i)))'*coefficients(Complicating(i))';  
        end
    end
end
sdpvar t

totalF = [totalF, totalobjective < t];
P = projection([totalF,t<1234],[allz;t]);

P = full(getbase(P));
b = P(:,1);
A = P(:,2:end);
rmv = find(A(:,end)== 0);
A(rmv,:)=[];
b(rmv) = [];
rmv = find(abs(b-1234)<1e-8);
A(rmv,:)=[];
b(rmv) = [];

p = basis(1)+LinearTerm+max((A(:,1:end-1)*allargs-b)./A(:,end));

