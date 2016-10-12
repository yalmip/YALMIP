function pcut = addEvalVariableCuts(p)

pcut = p;
if ~isempty(p.evalMap) 
    pcut = emptyNumericalModel;
    for i = 1:length(p.evalMap)
        y = p.evalVariables(i);
        x = p.evalMap{i}.variableIndex;
        xL = p.lb(x);
        xU = p.ub(x);
        
        % Generate a convex hull polytope
        if xL<xU
            if ~isempty(p.evalMap{i}.properties.convexhull)
                % A convex hull generator function is available!
                % Might be able to reuse hull from last run node
                if isfield(p.evalMap{i},'oldhull') && isequal(p.evalMap{i}.oldhull.xL,xL) && isequal(p.evalMap{i}.oldhull.xU,xU)
                    [Ax,Ay,b,K] = getOldHull(p,i);
                else
                    [Ax,Ay,b,K,p] = updateHull(xL,xU,p,i);
                    if isempty(Ax)
                        % Operator bounder does not cover this interval so
                        % use the sample-based instead
                        [Ax,Ay,b,K,p] = convexhullSampled(xL,xU,p,i);
                    end
                end              
            else               
               [Ax,Ay,b,K] = convexhullSampled(xL,xU,p,i);               
            end
            if ~isempty(b)
                if isempty(K)
                    % Compatibility with old code
                    K.f = 0;
                    K.l = length(b);
                end
                F_structemp = zeros(size(b,1),length(p.c)+1);
                F_structemp(:,1+y) = -Ay;
                F_structemp(:,1+x) = -Ax;
                F_structemp(:,1) = b;
                localModel = createNumericalModel(F_structemp,K);
                pcut = mergeNumericalModels(pcut,localModel);             
            end
        end
    end
    
    pcut = mergeNumericalModels(p,pcut);
end

function [Ax,Ay,b,K] = getOldHull(p,i);

Ax = p.evalMap{i}.oldhull.Ax;
Ay = p.evalMap{i}.oldhull.Ay;
b = p.evalMap{i}.oldhull.b;
K = p.evalMap{i}.oldhull.K;

function [Ax,Ay,b,K,p] = updateHull(xL,xU,p,i);

try
    [Ax,Ay,b,K]=feval(p.evalMap{i}.properties.convexhull,xL,xU, p.evalMap{i}.arg{2:end-1});
catch
    [Ax,Ay,b]=feval(p.evalMap{i}.properties.convexhull,xL,xU, p.evalMap{i}.arg{2:end-1});
    if ~isempty(Ax)
        problem = find(any(isinf([Ax Ay b]),2) | any(isnan([Ax Ay b]),2));
        Ax(problem,:) = [];
        Ay(problem,:) = [];
        b(problem) = [];
    end
    K = [];
end
p = saveOldHull(xL,xU,Ax,Ay,b,K,p,i);

function p = saveOldHull(xL,xU,Ax,Ay,b,K,p,i)
p.evalMap{i}.oldhull.xL = xL;
p.evalMap{i}.oldhull.xU = xU;
p.evalMap{i}.oldhull.Ax = Ax;
p.evalMap{i}.oldhull.Ay = Ay;
p.evalMap{i}.oldhull.b = b;
p.evalMap{i}.oldhull.K = K;

function [Ax,Ay,b,K,p] = convexhullSampled(xL,xU,p,i)

if length(xL)>1
    Ax = [];
    Ay = [];
    b = [];
    K = [];
    return
end
% sample function
z = linspace(xL,xU,100);

if isequal(p.evalMap{i}.fcn,'power_internal2')
    % Special code for automatically converting sigmonial
    % terms to be solvable with bmibnb
    fz = feval(p.evalMap{i}.fcn,z,p.evalMap{i}.arg{2});
    
else
    arg = p.evalMap{i}.arg;
    arg{1} = z;
    fz = real(feval(p.evalMap{i}.fcn,arg{1:end-1}));
    % end
    [minval,minpos] = min(fz);
    [maxval,maxpos] = max(fz);
    xtestmin = linspace(z(max([1 minpos-5])),z(min([100 minpos+5])),100);
    xtestmax = linspace(z(max([1 maxpos-5])),z(min([100 maxpos+5])),100);
    arg{1} = xtestmin;
    fz1 = real(feval(p.evalMap{i}.fcn,arg{1:end-1}));
    arg{1} = xtestmax;
    fz2 = real(feval(p.evalMap{i}.fcn,arg{1:end-1}));
    z = [z(:);xtestmin(:);xtestmax(:)];
    fz = [fz(:);fz1(:);fz2(:)];
    [z,sorter] = sort(z);
    fz = fz(sorter);
    [z,ii,jj]=unique(z);
    fz = fz(ii);
end

[Ax,Ay,b] = convexhullFromSampled(z,fz,xL,xU);
K = [];

p = saveOldHull(xL,xU,Ax,Ay,b,K,p,i);