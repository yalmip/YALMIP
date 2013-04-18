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
                                p.evalMap{i}.fcn = 'perspective_log';
                                p.evalMap{i}.arg{1} = recover([x;y]);
                                p.evalMap{i}.arg{2} = [];
                                p.evalMap{i}.variableIndex = [x y];
                                p.evalMap{i}.computes = other;
                                p.evalMap{i}.properties.bounds = @plog_bounds;
                                p.evalMap{i}.properties.convexhull = @plog_convexhull;
                                p.evalMap{i}.properties.derivative = @plog_derivative;
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

function dp = plog_derivative(x)
dp = [log(x(1)/x(2)) + 1;-x(1)/x(2)];

function [L,U] = plog_bounds(xL,xU)
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

function [Ax,Ay,b] = plog_convexhull(L,U);%xL,xU)
% 
% Ax = [];
% Ay = [];
% b = [];
% return

% L = [1 2];
% U = [3 6];

p1 = [L(1);L(2)];
p2 = [U(1);L(2)];
p3 = [L(1);U(2)];
p4 = [U(1);U(2)];
pm = [(L(1) + U(1))/2;(L(2) + U(2))/2];

z1 = L(1)*log(L(1)/L(2));
z2 = U(1)*log(U(1)/L(2));
z3 = L(1)*log(L(1)/U(2));
z4 = U(1)*log(U(1)/U(2));
zm = pm(1)*log(pm(1)/pm(2));

g1 = [log(L(1)/L(2)) + 1;-L(1)/L(2)];
g2 = [log(U(1)/L(2)) + 1;-U(1)/L(2)];
g3 = [log(L(1)/U(2)) + 1;-L(1)/U(2)];
g4 = [log(U(1)/U(2)) + 1;-U(1)/U(2)];
gm = [log(pm(1)/pm(2)) + 1;-pm(1)/pm(2)];
%sdpvar x y z
%C = [max([z1 z2 z3 z4]) > z > z1 + g1'*([x;y] - p1),z > z2 + g2'*([x;y] - p2), z > z3 + g3'*([x;y] - p3),z > z4 + g4'*([x;y] - p4),[x y]>U, L < [x y]]
Ax = [0 0;g1';g2';g3';g4';gm'];%;eye(2);-eye(2)];
Ay = [1;-1;-1;-1;-1;-1];%;0;0;0;0];
b = [max([z1 z2 z3 z4]);g1'*p1-z1;g2'*p2-z2;g3'*p3-z3;g4'*p4-z4;gm'*pm-zm];%;U(:);-L(:)];

xL = L(1);
yL = L(2);
xU = U(1);
yU = U(2);

p1 = [xL;yL;z1];
p2 = [xU;yL;z2];
p3 = [xL;yU;z3];
p4 = [xU;yU;z4];
U = max([z1 z2 z3 z4]);
H1 = null([p2-p1 p3-p1]')';
H2 = null([p3-p2 p4-p2]')';
A = [H1;H2];
b = [b;H1*p1;H2*p2];
Ax = [Ax;A(:,1:2)];
Ay = [Ay;A(:,3)];
% 
% for i = 1:100
%     for j = 1:100
%         xx(i,j) = xL + (xU-xL)*(i-1)/99;
%         yy(i,j) = yL + (yU-yL)*(j-1)/99;
%         zz(i,j) = xx(i,j)*log(xx(i,j)/yy(i,j));
%     end
% end
% xx=xx(:);
% yy=yy(:);
% zz = zz(:);
% maxx = -inf;
% for i = 1:length(xx)
%     maxx = max([maxx max(Ax*[xx(i);yy(i)] + Ay*zz(i) - b)]);
% end
% maxx






