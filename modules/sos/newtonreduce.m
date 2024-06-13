function [exponent_m,changed,no_lp_solved,keep] = newtonreduce(exponent_m,exponent_p,ops,interfacedata)
% NEWTONREDUCE Removes monomials outside half Newton polytope
%
% V = NEWTONREDUCE(V,P)
%
% Input
%  P : Scalar SDPVAR object
%  V : Vector with SDPVAR objects
%
% Output
%  V : Vector with SDPVAR objects
%
% Example:
%
% sdpvar x y
% p = 1+x^4*y^2+x^2*y^4;
% v = monolist([x y],degree(p)/2);
% sdisplay(v)
% v = newtonreduce(v,p);
% sdisplay(v)
%
% See also NEWTONMONOMS, CONSISTENT, CONGRUENCEBLOCKS

% Higher level call with SDPVAR polynomials
% Convert to exponent form
sdpvarout = 0;
if isa(exponent_m,'sdpvar')
    z = depends(exponent_p);
    z = recover(unique([depends(exponent_p) depends(exponent_m)]));
    [exponent_p,p_base] = getexponentbase(exponent_p,z);
    [m,m_base] = getexponentbase(exponent_m,z);
    exponent_m = cell(1);exponent_m{1} = m;
    sdpvarout = 1;
end

% Put exponents in cell (compability with monomialreduction.m
doubleout = 0;
if isa(exponent_m,'double')
    temp = exponent_m;
    exponent_m = cell(1);
    exponent_m{1} = temp;
    doubleout = 1;
end

if nargin<3
    ops = sdpsettings;
end

ops.solver = 'cdd,glpk,*';  % CDD is generally robust on these problems
ops.verbose = 0;
ops.saveduals = 0;

% Create an LP structure
if nargin < 4
    x=sdpvar(1,1);
    [aux1,aux2,aux3,interfacedata] = export((x>=0),x,ops);
    interfacedata.getsolvertime = 0;
end

bfixed = [zeros(size(exponent_p,1),1);1];
interfacedata.K.l = length(bfixed);
interfacedata.lb = [];
interfacedata.ub = [];

no_lp_solved = 0;
for j = 1:length(exponent_m)
    try
        % Basic idea : Try to find separating hyperplane between Newton polytope and candidate,
        % for all candidates.
        % Pro : Numerical stability, polynomial (does NOT calculate H-representation of Newton polytope)
        % Con : Slightly slower than explicit convex hull calculation in most cases.

        
        keep = ismember(exponent_m{j}*2,exponent_p,'rows');
        dontkeep = 0*keep;
        
        if 0
            % Requires that cdd is installed.
            H=cddmex('hull',struct('V',full(exponent_p)));
            for i = 1:size(exponent_m{j},1)
                if ~all(H.A*exponent_m{j}(i,:)'*2 <= H.B)
                    dontkeep(i) = 1;
                end
            end
        else
            ii = 0;
            
            for i = 1:length(keep)
                if ~keep(i) & ~dontkeep(i)
                    q = exponent_m{j}(i,:)'*2;
                    interfacedata.c = ([-q' 1])';
                    interfacedata.Q = spalloc(length(interfacedata.c),length(interfacedata.c),0);
                    interfacedata.variabletype = zeros(1,length(interfacedata.c));
                    interfacedata.F_struc = ([bfixed -[exponent_p -ones(size(exponent_p,1),1);q' -1]]);
                    solution = feval(interfacedata.solver.call,interfacedata);
                    no_lp_solved = no_lp_solved  + 1;
                    ii = ii + 1;
                    if (solution.problem == 0 & ([-q' 1]*solution.Primal < 0)) | solution.problem == 2
                        a = solution.Primal(1:end-1);
                        b = solution.Primal(end);
                        u = find(a'*2*exponent_m{j}' - b > sqrt(eps));
                        dontkeep(u) = 1;
                    end
                end
                if keep(i)
                    q = exponent_m{j}(i,:)'*2;
                end
            end
        end
        keep = find(~dontkeep);
    catch
        keep = newtonpolytope(exponent_m{j}*2,exponent_p);
    end
    changed = length(keep)<size(exponent_m{j},1);
    exponent_m{j} = exponent_m{j}(keep,:);
end

if doubleout
    exponent_m = exponent_m{1};
end

if sdpvarout
    exponent_m = recovermonoms(exponent_m{1},z);
end

