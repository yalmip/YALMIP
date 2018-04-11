function [Fnew, obj, M,k,x,u,n,deg,linears,nonlinears,vecConstraints,isinequality,ulong] = momentmodel(F,obj,k,keepnonlinears)

if nargin < 2
    obj = [];
end
if nargin < 3
    k = [];
end
if nargin < 4
    keepnonlinears = 0;
end

vecConstraints = [];
sdpConstraints = [];
isinequality = [];
binaries = [];
xvars = [];
Fnew = ([]);
for i = 1:length(F)
    if is(F(i),'elementwise')
        X = sdpvar(F(i));
        vecConstraints = [vecConstraints;X(:)];
        isinequality = [isinequality ones(1,prod(size(X)))];
        xvars = [xvars depends(X(:))];
    elseif is(F(i),'equality')
        X = sdpvar(F(i));
        if is(X,'symmetric')
            X = X(find(triu(ones(length(X)))));
        end
        vecConstraints = [vecConstraints;-X(:)];
        isinequality = [isinequality zeros(1,prod(size(X)))];
        xvars = [xvars depends(X(:))];
    elseif is(F(i),'sdp')
        sdpConstraints{end+1} = sdpvar(F(i));
        xvars = [xvars depends(F(i))];
    elseif is(F(i),'binary')
        binaries = [binaries getvariables(F(i))];
    else
        Fnew = Fnew+F(i); % Should only be SOCP constraints
    end
end

% Recover the involved variables
x = recover(unique([depends(obj) xvars]));
n = length(x);

% Check degrees of constraints
deg = [];
for i = 1:length(vecConstraints)
    deg(end+1) = degree(vecConstraints(i));
end
for i = 1:length(sdpConstraints)
    deg(end+1) = degree(sdpConstraints{i});
end
if isempty(deg)
    deg = 0;
end

% Create lowest possible relaxation if k=[]
d = ceil((max(degree(obj),max(deg)))/2);
k_min = d;
if isempty(k)
    k = k_min;
else
    if k<k_min
        error('Higher order relaxation needed')
    end
end

% Generate monomials of order k
u{k} = monolist(x,k);
ulong{k} = monolist(x,2*k);

% Largest moment matrix. NOTE SHIFT M{k+1} = M_k.
M{k+1}=u{k}*u{k}';
% Moment matrices easily generated with this trick
% The matrices will NOT be rank-1 since the products
% generate the relaxed variables

% ... and lower degree localization matrices
M{1} = 1;
for i = 1:1:k-1;
    n_i = round(factorial(n+k-i)/(factorial(n)*factorial(k-i)));
    M{k-i+1} = M{k+1}(1:n_i,1:n_i);
end

% Lasserres relaxation (Lasserre, SIAM J. OPTIM, 11(3) 796-817)
Fmoments = (M{k+1}>=0);
for i = 1:length(vecConstraints)   
    if isinequality(i)
        v_k = floor((degree(vecConstraints(i))+1)/2);
        Localizer = vecConstraints(i)*M{k-v_k+1};
        if isa(vecConstraints(i),'double')
            if vecConstraints(i)<0
                error('Problem is trivially infeasible due to negative constant')
            else
                continue
            end
        end
        Fmoments = Fmoments+(Localizer>=0);
    else
        if isa(vecConstraints(i),'double')
            if vecConstraints(i)~=0
                error('Problem is trivially infeasible due to non-zero constant in equality constraints')
            else
                continue
            end
        end        
        Localizer = vecConstraints(i)*monolist(x,2*k-degree(vecConstraints(i)));      
        Fmoments = Fmoments+(Localizer==0);
    end
end
for i = 1:length(sdpConstraints)
    v_k = floor((degree(sdpConstraints{i})+1)/2);
    Fmoments = Fmoments+(kron(M{k-v_k+1},sdpConstraints{i})>=0);
end

% Add them all
Fnew = Fnew + Fmoments;

% Get all binary and reduce problem
binaries = union(binaries,yalmip('binvariables'));
if ~isempty(binaries)
    if isa(obj,'sdpvar')        
        obj = eliminateBinary(obj,binaries);
    end
    for i = 1:length(Fmoments)
        Fnew(i) = eliminateBinary(Fnew(i),binaries);
    end
    for i = 2:1:k+1;
        M{i} = eliminateBinary(M{i},binaries);
    end
end

vars = getvariables(Fnew);
for i = 1:length(M)
    vars = [vars getvariables(M{i})];
end
vars = unique([vars getvariables(obj)]);
[mt,variabletype] = yalmip('monomtable');
nonlinears = vars(find(variabletype(vars)));
newLinear = sdpvar(length(nonlinears),1);

if isa(obj,'sdpvar')
    obj = variablereplace(obj,nonlinears,getvariables(newLinear));
end
for i = 1:length(M)
    if isa(M{i},'sdpvar')
        M{i} = variablereplace(M{i},nonlinears,getvariables(newLinear));
    end
end
linears = getvariables(newLinear);
Fnew = variablereplace(Fnew,nonlinears,getvariables(newLinear));
linears = recover(linears);
nonlinears = recover(nonlinears);

end
