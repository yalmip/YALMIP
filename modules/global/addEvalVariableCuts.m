function pcut = addEvalVariableCuts(p)

temp = p.F_struc(1:p.K.f,:);
p.F_struc(1:p.K.f,:) = [];
pcut = p;
for i = 1:length(p.evalMap)
    y = p.evalVariables(i);
    x = p.evalMap{i}.variableIndex;
    xL = p.lb(x);
    xU = p.ub(x);

    % Generate a convex hull polytope
    if xL<xU
        if ~isempty(p.evalMap{i}.properties.convexhull)
            % A convex hull generator function is available!
            K = [];
            try
                if isfield(p.evalMap{i},'oldhull') & isequal(p.evalMap{i}.oldhull.xL,xL) & isequal(p.evalMap{i}.oldhull.xU,xU)
                    Ax = p.evalMap{i}.oldhull.Ax;
                    Ay = p.evalMap{i}.oldhull.Ay;
                    b = p.evalMap{i}.oldhull.b;
                    K = p.evalMap{i}.oldhull.K;
                else
                    [Ax,Ay,b,K]=feval(p.evalMap{i}.properties.convexhull,xL,xU, p.evalMap{i}.arg{2:end-1});
                    pcut.evalMap{i}.oldhull.xL = xL;
                    pcut.evalMap{i}.oldhull.xU = xU;
                    pcut.evalMap{i}.oldhull.Ax = Ax;
                    pcut.evalMap{i}.oldhull.Ay = Ay;
                    pcut.evalMap{i}.oldhull.b = b;
                    pcut.evalMap{i}.oldhull.K = K;
                end
            catch
                if isfield(p.evalMap{i},'oldhull') & isequal(p.evalMap{i}.oldhull.xL,xL) & isequal(p.evalMap{i}.oldhull.xU,xU)
                    Ax = p.evalMap{i}.oldhull.Ax;
                    Ay = p.evalMap{i}.oldhull.Ay;
                    b = p.evalMap{i}.oldhull.b;
                    K = p.evalMap{i}.oldhull.K;
                else
                    [Ax,Ay,b]=feval(p.evalMap{i}.properties.convexhull,xL,xU, p.evalMap{i}.arg{2:end-1});
                    pcut.evalMap{i}.oldhull.xL = xL;
                    pcut.evalMap{i}.oldhull.xU = xU;
                    pcut.evalMap{i}.oldhull.Ax = Ax;
                    pcut.evalMap{i}.oldhull.Ay = Ay;
                    pcut.evalMap{i}.oldhull.b = b;
                    pcut.evalMap{i}.oldhull.K = K;
                end
            end
        else
            if length(xL)>1
                disp(['The ' p.evalMap{i}.fcn ' operator does not have a convex hull operator'])
                disp('This is required for multi-input single output operators');
                disp('Sampling approximation does not work in this case.');
                error('Missing convex hull operator');
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
                fz = feval(p.evalMap{i}.fcn,arg{1:end-1});
                % end
                [minval,minpos] = min(fz);
                [maxval,maxpos] = max(fz);
                xtestmin = linspace(z(max([1 minpos-5])),z(min([100 minpos+5])),100);
                xtestmax = linspace(z(max([1 maxpos-5])),z(min([100 maxpos+5])),100);
                arg{1} = xtestmin;
                fz1 = feval(p.evalMap{i}.fcn,arg{1:end-1});
                arg{1} = xtestmax;
                fz2 = feval(p.evalMap{i}.fcn,arg{1:end-1});
                z = [z(:);xtestmin(:);xtestmax(:)];
                fz = [fz(:);fz1(:);fz2(:)];
                [z,sorter] = sort(z);
                fz = fz(sorter);
                [z,ii,jj]=unique(z);
                fz = fz(ii);
            end

            %                    fz = feval(p.evalMap{i}.fcn,z);
            %   end
            % create 4 bounding planes
            % f(z) < k1*(x-XL) + f(xL)
            % f(z) > k2*(x-XL) + f(xL)
            % f(z) < k3*(x-XU) + f(xU)
            % f(z) > k4*(x-XU) + f(xU)
            k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
            k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
            k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
            k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
            Ax = [-k1;k2;-k3;k4];
            Ay = [1;-1;1;-1];
            b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
            K = [];
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
            pcut.F_struc = [F_structemp; pcut.F_struc];
            pcut.K.l = pcut.K.l + K.l;%length(b);
            pcut.K.f = pcut.K.f + K.f;%length(b);
        end
    end
end
pcut.F_struc = [temp;pcut.F_struc];


