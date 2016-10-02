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




