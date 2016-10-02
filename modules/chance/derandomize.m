function Fderandomized = derandomize(F)

% Goes through all probability constraints and checks for cases where we
% can use analytic expressions.

% Find chance constraints
chanceDeclarations = find(is(F,'chance'));
% Find variables with attached distributions
randomDeclarations = find(is(F,'random'));

if isempty(randomDeclarations)
    Fderandomized = F;
    return
end

keep = ones(length(F),1);
%keep(randomDeclarations)=0;
keep(chanceDeclarations)=0;
randomVariables = extractRandomDefinitions(F(randomDeclarations));
[randomVariables,map] = mergeDistributions(randomVariables);
groupedChanceConstraints = groupchanceconstraints(F);

[Fderandomized,eliminatedConstraints,recursive] = deriveChanceModel(groupedChanceConstraints,randomVariables);
Fderandomized = Fderandomized + F(find(keep)) + F(find(keep(~eliminatedConstraints)));
if recursive
    Fderandomized = derandomize(Fderandomized)
end

function [Fderandomized,eliminatedConstraints,recursive] = deriveChanceModel(groupedChanceConstraints,randomVariables);

recursive = 0;
Fderandomized = [];
eliminatedConstraints = zeros(length(groupedChanceConstraints),1);

allwVars = [];
for i = 1:length(randomVariables)
    allwVars = [allwVars;getvariables(randomVariables{i}.variables)];
end

for uncertaintyGroup = 1:length(randomVariables)
    
    wVars = getvariables(randomVariables{uncertaintyGroup}.variables);
    
    for ic = 1:length(groupedChanceConstraints)
        if length(groupedChanceConstraints{ic})>1
            error('Joint chance constraint not supported');
        end
        if ~is(groupedChanceConstraints{ic},'elementwise')
            error('Only elementwise chance constraints supported')
        end
        X = sdpvar(groupedChanceConstraints{ic});
        if length(X)>1
            error('Only single elementwise chance constraints supported')
        end
        
        % Extract quadratic part, X = fX + X, where fx is other stuff
        [fX,X] = functionSeparation(X);
        
        allVars = depends(X);
        if ~isempty(intersect(wVars,allVars))
            xVars = setdiff(allVars,wVars);
            x = recover(xVars);
            w = recover(wVars);
            
            fail = 0;
            [A,cx,b,cw,fail] = quadraticDecomposition(X,x,w);
            
            % Remap to original ordering on variables in distribution
            % Base = full(getbase(randomVariables{uncertaintyGroup}.variables));
            % Base = Base(:,2:end);
            % [ii,jj,kk] = find(Base)
            %  cw = cw*Base;
            %  A = A*Base;
            
            % b(x) + c(x)'*w >= 0
            b = b + fX + cx'*x;
            c = A'*x + cw';
            
            newConstraint = [];            
            if ~fail
                confidencelevel = struct(groupedChanceConstraints{ic}).clauses{1}.confidencelevel;
                if strcmp(func2str(randomVariables{uncertaintyGroup}.distribution.name),'random')
                    switch randomVariables{uncertaintyGroup}.distribution.parameters{1}
                        case 'moment'
                            newConstraint = momentChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,confidencelevel);
                            eliminatedConstraints(ic)=1;
                        case {'normal','normalm'}
                            newConstraint = normalChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,confidencelevel);
                            eliminatedConstraints(ic)=1;
                        case 'normalf'
                            newConstraint = normalfactorizedChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,confidencelevel);
                            eliminatedConstraints(ic)=1;
                        otherwise
                            newConstraint =  sampledchernoffChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,confidencelevel,w);
                            eliminatedConstraints(ic)=1;
                    end
                else
                    newConstraint =  sampledmarkovChanceFilter(b,c,randomVariables{uncertaintyGroup}.distribution,confidencelevel,w);
                    eliminatedConstraints(ic)=1;
                end
            end
            if ~isempty(newConstraint)
                if ~isempty(intersect(depends(newConstraint),allwVars))
                    % New uncertainties popped up,i.e. parameters in a
                    % distribution, are distributions them selves
                    Fderandomized = [Fderandomized, probability(newConstraint)>=confidencelevel];
                    recursive = 1;
                else
                    Fderandomized = [Fderandomized, newConstraint];
                end
            end
        end
    end
end

function [merged,map] = mergeDistributions(randomVariables);

merged = {};
map = [];
used = zeros(1,length(randomVariables));
mergedindex = 1;
for i = 1:length(randomVariables)
    if ~used(i)
        this = randomVariables{i};
        used(i) = 1;
        for j = i+1:length(randomVariables)
            if ~used(j)
                that = randomVariables{j};
                if strcmp(func2str(this.distribution.name),'random') && strcmp(func2str(that.distribution.name),'random')
                    if mergable(this.distribution.parameters{1},that.distribution.parameters{1})
                        % Same distribution
                        this = merge2distributions(this,that);
                        map(j) = mergedindex;
                        used(j) = 1;
                    end
                end
            end
        end
        merged{mergedindex} = this;
        map(i) = mergedindex;
        mergedindex = mergedindex + 1;
    end
end


function C = merge2distributions(A,B)

C = A;
C.variables = [A.variables;B.variables];
if any(strcmp(A.distribution.parameters{1},{'normal','normalf','normalm'})) || any(strcmp(A.distribution.parameters{1},{'normal','normalf','normalm'}))
    
    % A bit messy as there are two additional forms for normal
    % distributions, used to define mv normals, and normals with a
    % factorized covariance
    
    % First, normalize to matrix format
    dimAvar = size(A.distribution.parameters{3});
    dimBvar = size(B.distribution.parameters{3});
    if dimAvar(1) ~= dimAvar(2)
        varA = diag(A.distribution.parameters{3});
    else
        varA = A.distribution.parameters{3};
    end
    if dimBvar(1) ~= dimBvar(2)
        varB = diag(B.distribution.parameters{3});
    else
        varB = B.distribution.parameters{3};
    end
    % Is A in factor form, but not B, and vice versa. If so, put the other
    % one in factor form too
    if strcmp(A.distribution.parameters{1},'normalf') &&  any(strcmp(B.distribution.parameters{1},{'normal','normalm'}))
        varB = chol(varB);
    elseif strcmp(B.distribution.parameters{1},'normalf') &&  any(strcmp(A.distribution.parameters{1},{'normal','normalm'}))
        varA = chol(varA);
    end
    varC = blkdiag(varA,varB);
    
    % Sort out various combinations of normal stuff in different forms
    if strcmp(A.distribution.parameters{1},'normal') && strcmp(B.distribution.parameters{1},'normal')
        varC = diag(varC);
        C.distribution.parameters{1} = 'normal';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = diag(varC);
    elseif strcmp(A.distribution.parameters{1},'normalf') ||  any(strcmp(B.distribution.parameters{1},'normalf'))
        C.distribution.parameters{1} = 'normalf';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = varC;
    else
        C.distribution.parameters{1} = 'normalm';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = varC;
    end
else
    for k = 2:length(A.distribution.parameters)
        C.distribution.parameters{k} = [A.distribution.parameters{k};B.distribution.parameters{k}];
    end
end


function yesno = mergable(a,b)
if isequal(a,b)
    yesno = 1;
elseif any(strcmp(a,{'normal','normalf','normalm'})) && any(strcmp(b,{'normal','normalf','normalm'}))
    yesno = 1;
else
    yesno = 0;
end


function [AAA,ccc,b,c_wTbase,fail] = quadraticDecomposition(X,x,w)
b = [];
A = [];
% Some pre-calc
xw = [x;w];
fail = 0;
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
[Qs,cs,fs,dummy,nonquadratic] = vecquaddecomp(X,xw);
c_wTbase = [];
AAA = [];
ccc = [];
for i = 1:length(X)
    Q = Qs{i};
    c = cs{i};
    f = fs{i};
    if nonquadratic
        error('Constraints can be at most quadratic, with the linear term uncertain');
    end
    Q_ww = Q(wind,wind);
    if nnz(Q_ww)>0
        fail = 1;
        return
    end
    Q_xw = Q(xind,wind);
    Q_xx = Q(xind,xind);
    c_x = c(xind);
    c_w = c(wind);
    
    %b = [b;f + c_w'*w];
    %A = [A;-c_x'-w'*2*Q_xw'];
    % A = [A -c_x-2*Q_xw*w];
    AAA = [AAA;sparse(2*Q_xw)];
    ccc = [ccc;sparse(c_x)];
    b = [b;f+x'*Q_xx*x];
    c_wTbase = [c_wTbase;c_w'];
end


function model = momentChanceFilter(b,c,distribution,confidencelevel)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
X = covariance + theMean*theMean';
S = chol(X);
gamma = sqrtm(confidencelevel);
e = [S*c+b*(inv(S')*theMean);b*sqrtm(1-theMean'*inv(X)*theMean)];
if norm(full(getbase(e'*e - (c'*X*c + 2*b*c'*theMean + b^2))))>1e-10
    error
end
model =  b + c'*theMean >= gamma*norm_callback(e);
%model =  b + c'*theMean*0 >= gamma*sqrtm(c'*covariance*c-2*b*theMean'*c*0+b^2);

function newConstraint = normalChanceFilter(b,c,distribution,confidencelevel)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
gamma = icdf('normal',confidencelevel,0,1);
if min(size(covariance))==1
    covariance = diag(covariance);
end
if isa(covariance,'sdpvar')
    error('Covariance cannot be an SDPVAR in normal distribution. Maybe you meant to use factorized covariance in ''normalf''');
end

e = chol(covariance)*c;
if isa(gamma,'sdpvar')
    newConstraint = b + c'*theMean >= gamma*norm_callback(e);
else
    newConstraint =  b + c'*theMean >= gamma*norm(e);
end

function newConstraint = normalfactorizedChanceFilter(b,c,distribution,confidencelevel);
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
gamma = icdf('normal',confidencelevel,0,1);
e = covariance*c;
if isa(gamma,'sdpvar')
    newConstraint = b + c'*theMean >= gamma*norm_callback(e);
else
    newConstraint = b + c'*theMean >= gamma*norm(e);
end


function newConstraint =  sampledmomentChanceFilter(b,c,distribution,confidencelevel,w);
W = [];for i = 1:1000;W = [W dataSampler(distribution,size(w))];end
d.parameters{2} = mean(W,2);
d.parameters{3} = cov(W');;
newConstraint = momentChanceFilter(b,c,d,confidencelevel);


function newConstraint =  sampledmarkovChanceFilter(b,c,distribution,confidencelevel,w);
W = [];for i = 1:100;W = [W dataSampler(distribution,size(w))];end
alpha = sdpvar(1);
s = sdpvar(1,100);
newConstraint = [0 <= alpha, -b-c'*W+alpha <= s, 0 <= s, sum(s)/length(s) <= alpha*(1-confidencelevel)]

function newConstraint =  sampledchernoffChanceFilter(b,c,distribution,confidencelevel,w);
N = 250;
W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
alpha = sdpvar(1);
e = pexpsum([repmat(alpha,1,N);-b-c'*W])/N;
newConstraint = [0 <= alpha, e <= alpha*(1-confidencelevel)];

