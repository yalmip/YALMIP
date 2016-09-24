function p = convert_perspective_log(p)

p.kept = 1:length(p.c);
if isempty(p.evalMap)
    return
end

if ~any(p.variabletype == 4)
    return
end

variableIndex = [];
for i=1:length(p.evalMap)
    variableIndex = [variableIndex p.evalMap{i}.variableIndex];
end

removable = [];
for i = 1:length(p.evalMap)
    if strcmp(p.evalMap{i}.fcn,'log')
        argument = p.evalMap{i}.variableIndex;
        if length(argument) == 1
            if p.variabletype(argument)==4
                monoms = p.monomtable(argument,:);
                if nnz(monoms) == 2
                    k = find(monoms);
                    p1 = monoms(k(1));
                    p2 = monoms(k(2));
                    if isequal(sort([p1 p2]) , [-1 1])
                        if p2>p1
                            x = k(2);
                            y = k(1);
                        else
                            x = k(1);
                            y = k(2);
                        end
                        % Ok, so we have log(x/y)
                        % is this multiplied by x somewhere
                        logxy = p.evalMap{i}.computes;
                        enters_in = find(p.monomtable(:,logxy));
                        other = setdiff(enters_in,p.evalMap{i}.computes);
                        if length(nnz(other)) == 1
                            monomsxlog = p.monomtable(other,:);
                            if nnz(monomsxlog) == 2 & (monomsxlog(x) == 1)
                                % Hey, x*log(x/y)!
                                % we change this monomial variable to a
                                % callback variable
                                p.evalMap{i}.fcn = 'negated_perspective_log';
                                p.evalMap{i}.arg{1} = recover([x;y]);
                                p.evalMap{i}.arg{2} = [];
                                p.evalMap{i}.variableIndex = [x y];
                                p.evalMap{i}.computes = other;
                                p.evalMap{i}.properties.bounds = @nplog_bounds;
                                p.evalMap{i}.properties.convexhull = @nplog_convexhull;
                                p.evalMap{i}.properties.derivative = @nplog_derivative;
                                p.evalMap{i}.properties.inverse = [];
                                p.variabletype(other) = 0;
                                p.monomtable(other,:) = 0;
                                p.monomtable(other,other) = 1;
                                p.evalVariables(i) = other;
                                % Figure out if x/y can be removed
                                % This is possible if the x/y term never is
                                % used besides inside the log term
                                if nnz(p.F_struc(:,1+argument)) == 1 & p.c(argument) == 0 & nnz(argument == variableIndex) == 1
                                    removable = [removable argument];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
kept = 1:length(p.c);
kept = setdiff(kept,removable);
aux_used = zeros(1,length(p.c));
aux_used(p.aux_variables) = 1;
aux_used(removable)=[];
p.aux_variables = find(aux_used);
if length(removable) > 0
    kept = 1:length(p.c);
    kept = setdiff(kept,removable);    
    [ii,jj,kk] = find(p.F_struc(:,1+removable));
    p.F_struc(:,1+removable) = [];    
    p.F_struc(ii,:) = [];
    p.K.l = p.K.l - length(removable);
    p.c(removable) = [];
    p.Q(removable,:) = [];
    p.Q(:,removable) = [];
    p.variabletype(removable) = [];
    p.monomtable(:,removable) = [];
    p.monomtable(removable,:) = [];
    for i = 1:length(p.evalVariables)
        p.evalVariables(i) = find(p.evalVariables(i) == kept);
        for j = 1:length(p.evalMap{i}.variableIndex)
           p.evalMap{i}.variableIndex(j) = find(p.evalMap{i}.variableIndex(j) == kept);
        end
        for j = 1:length(p.evalMap{i}.computes)
           p.evalMap{i}.computes(j) = find(p.evalMap{i}.computes(j) == kept);
        end        
    end
    p.lb(removable) = [];
    p.ub(removable) = [];
    p.used_variables(removable) = [];
end

function dp = nplog_derivative(x)
dp = [log(x(1)/x(2)) + 1;-x(1)/x(2)];

function [L,U] = nplog_bounds(xL,xU)
xU(isinf(xU)) = 1e12;
x1 = xL(1)*log(xL(1)/xL(2));
x2 = xU(1)*log(xU(1)/xL(2));
x3 = xL(1)*log(xL(1)/xU(2));
x4 = xU(1)*log(xU(1)/xU(2));
U = max([x1 x2 x3 x4]);
if (exp(-1)*xU(2) > xL(1)) & (exp(-1)*xU(2) < xU(1))
    L = -exp(-1)*xU(2);
else
    L = min([x1 x2 x3 x4]);
end

function [Ax,Ay,b] = nplog_convexhull(xL,xU);

x1 = [xL(1);xL(2)];
x2 = [xU(1);xL(2)];
x3 = [xL(1);xU(2)];
x4 = [xU(1);xU(2)];
x5 = (xL+xU)/2;

f1 = negated_perspective_log(x1);
f2 = negated_perspective_log(x2);
f3 = negated_perspective_log(x3);
f4 = negated_perspective_log(x4);
f5 = negated_perspective_log(x5);

df1 = nplog_derivative(x1);
df2 = nplog_derivative(x2);
df3 = nplog_derivative(x3);
df4 = nplog_derivative(x4);
df5 = nplog_derivative(x5);

[Ax,Ay,b] = convexhullConvex2D(x1,f1,df1,x2,f2,df2,x3,f3,df3,x4,f4,df4,x5,f5,df5);