function model = mpt_parbb(Matrices,options)

% For simple development, the code is currently implemented in high-level
% YALMIP and MPT code. Hence, a substantial part of the computation time is
% stupid over-head.

[Matrices.lb,Matrices.ub] = mpt_detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);

U = sdpvar(Matrices.nu,1);
x = sdpvar(Matrices.nx,1);
F = set(Matrices.G*U < Matrices.W + Matrices.E*x);
F = F + set(Matrices.lb < [U;x] < Matrices.ub);
F = F + set(binary(U(Matrices.binary_var_index)));
F = F + set(Matrices.Aeq*U + Matrices.Beq*x == Matrices.beq);
h = Matrices.H*U + Matrices.F*x;

Universe = polytope(set(Matrices.lb(end-Matrices.nx+1:end) <= x <= Matrices.ub(end-Matrices.nx+1:end)));

model = parametric_bb(F,h,options,x,Universe);

function sol = parametric_bb(F,obj,ops,x,Universe)

% F   : All constraints
% obj : Objective
% x   : parametric variables
% y   : all binary variables

if isempty(ops)
    ops = sdpsettings;
end
ops.mp.algorithm = 1;
ops.cachesolvers = 0;
ops.mp.presolve=1;
ops.solver = '';

% Expand nonlinear operators only once
F = expandmodel(F,obj);
ops.expand = 0;

% Gather all binary variables
y = unique([depends(F) depends(obj)]);
n = length(y)-length(x);
y = intersect(y,[yalmip('binvariables') depends(F(find(is(F,'binary'))))]);
y = recover(y);

% Make sure binary relaxations satisfy 0-1 constraints
F = F + set(0 <= y <= 1);

% recursive, starting in maximum universe
sol = sub_bb(F,obj,ops,x,y,Universe);

% Nice, however, we have introduced variables along the process, so the
% parametric solutions contain variables we don't care about
for i = 1:length(sol)
    for j = 1:length(sol{i}.Fi)
        sol{i}.Fi{j} = sol{i}.Fi{j}(1:n,:);
        sol{i}.Gi{j} = sol{i}.Gi{j}(1:n,:);
    end
end

function sol = sub_bb(F,obj,ops,x,y,Universe)

sol = {};

% Find a feasible point in this region. Note that it may be the case that a
% point is feasible, but the feasible space is flat. This will cause the
% mplp solver to return an empty solution, and we have to pick a new
% binary solution.

localsol = {[]};
intsol.problem = 0;

if 1%while intsol.problem == 0
    localsol = {[]};
    while isempty(localsol{1}) & (intsol.problem == 0)
        ops.verbose = ops.verbose-1;
        intsol = solvesdp(F,obj,sdpsettings(ops,'solver','glpk'));
        ops.verbose = ops.verbose+1;
        if intsol.problem == 0
            y_feasible = round(double(y));
            ops.relax = 1;
            localsol = solvemp(F+set(y == y_feasible),obj,ops,x);
            ops.relax = 0;
            if isempty(localsol{1})
                F = F + not_equal(y,y_feasible);
            end
            F = F + not_equal(y,y_feasible);

        end
    end
    if ~isempty(localsol{1})

        % YALMIP syntax...
        if isa(localsol,'cell')
            localsol = localsol{1};
        end

        % Now we want to find solutions with other binary combinations, in
        % order to find the best one. Cut away the current bionary using
        % overloaded not equal
        F = F + not_equal(y,y_feasible);

        % Could be that the binary was feasible, but the feasible space in the
        % other variables is empty/lower-dimensional
        if ~isempty(localsol)
            % Dig into this solution. Try to find another feasible binary
            % combination, with a better cost, in each of the regions
            for i = 1:length(localsol.Pn)
                G = F;
                % Better cost
                G = G + set(obj <= localsol.Bi{i}*x + localsol.Ci{i});
                % In this region
                [H,K] = double(localsol.Pn(i));
                G = G + set(H*x <= K);
                % Recurse
                diggsol{i} = sub_bb(G,obj,ops,x,y,localsol.Pn(i));
            end

            % Create all neighbour regions, and compute solutions in them too
            flipped = regiondiff(union(Universe),union(localsol.Pn));
            flipsol={};
            for i = 1:length(flipped)
                [H,K] = double(flipped(i));
                flipsol{i} = sub_bb(F+ set(H*x <= K),obj,ops,x,y,flipped(i));
            end

            % Just place all solutions in one big cell. We should do some
            % intersect and compare already here, but I am lazy now.
            sol = appendlists(sol,localsol,diggsol,flipsol);
        end
    end
end

function sol = appendlists(sol,localsol,diggsol,flipsol)

sol{end+1} = localsol;
for i = 1:length(diggsol)
    if ~isempty(diggsol{i})
        if isa(diggsol{i},'cell')
            for j = 1:length(diggsol{i})
                sol{end+1} = diggsol{i}{j};
            end
        else
            sol{end+1} = diggsol{i};
        end
    end
end
for i = 1:length(flipsol)
    if ~isempty(flipsol{i})
        if isa(flipsol{i},'cell')
            for j = 1:length(flipsol{i})
                sol{end+1} = flipsol{i}{j};
            end
        else
            sol{end+1} = flipsol{i};
        end
    end
end


function F = not_equal(X,Y)
zv = find((Y == 0));
ov = find((Y == 1));
lhs = 0;
if ~isempty(zv)
    lhs = lhs + sum(extsubsref(X,zv));
end
if ~isempty(ov)
    lhs = lhs + sum(1-extsubsref(X,ov));
end
F = set(lhs >=1);
