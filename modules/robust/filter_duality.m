function [F,feasible] = filter_duality(F_xw,Zmodel,x,w,ops)

% Creates robustified version of the uncertain set of linear inequalities
% s.t A(w)*x <= b(w) for all F(w) >= 0 where F(w) is a conic set, here
% given in YALMIP numerical format.
%
% Based on Robust Optimization - Methodology and Applications. A. Ben-Tal
% and A. Nemerovskii. Mathematical Programming (Series B), 92:453-480, 2002
% Note, there are some sign errors in the paper.
%
% The method introduces a large amount of new variables, equalities and
% inequalities. By turning on robust.reducedual, the equalities are
% eliminated, thus reducing the number of variables slightly

feasible = 1;

if length(F_xw) == 0
    F = [];
    return
end

X = sdpvar(F_xw);
b = [];
A = [];
% Some pre-calc
xw = [x;w];
xind = find(ismembcYALMIP(getvariables(xw),getvariables(x)));
wind = find(ismembcYALMIP(getvariables(xw),getvariables(w)));
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
    Q_xw = Q(xind,wind);
    Q_xx = Q(xind,xind);
    c_x = c(xind);
    c_w = c(wind);

    %b = [b;f + c_w'*w];
    %A = [A;-c_x'-w'*2*Q_xw'];
    % A = [A -c_x-2*Q_xw*w];
    AAA = [AAA;sparse(-2*Q_xw)];
    ccc = [ccc;-sparse(c_x)];
    b = [b;f];
    c_wTbase = [c_wTbase;c_w'];       
end
b = b + c_wTbase*w;
% Ac = A';
A = reshape(ccc + AAA*w,size(c_x,1),[]);
if isa(A,'double')
    A = sparse(A);
end
A = A';
% Try to find variables that only have simple bound constraints. These
% variables can explicitly be optimized and thus speed up the construction,
% and allow a model with fewer variables.
[Zmodel2,lower,upper] = find_simple_variable_bounds(Zmodel);

% Partition the uncertain variables
simple_w  = find( ~isinf(lower) & ~isinf(upper));
general_w = find( isinf(lower) |  isinf(upper));
simple_w = recover(simple_w);
general_w = recover(general_w);

% Linear uncertain constraint is (Bbetai*x + cdi) >= 0 for all w, or
% (bi' + (Bi*w)')*x + (ci'*w + di).
cd    = b;
Bbeta = -A;

F = ([]);
top = 1;

% To speed up the construction, compute the ci vectors for all constraints
% in one call ci_basis = [c1 c2 ...]
ci_basis = basis(cd',w);

if ops.verbose
    disp(' - Using duality to eliminate uncertainty');
end

nv = yalmip('nvars');
simple_rows = [];
for i = 1:length(b)
   
    Bbetai = Bbeta(i,:);

    if (nnz(ci_basis(:,i))==0) & isa(Bbetai,'double')
        % This constraint row does not depend on uncertainty
        %row = Bbetai*x + cdi;
        row = X(i);
        if isa(row,'sdpvar')
          %  F = F + (row >= 0);
            simple_rows =  [simple_rows;i];
        else
            if row<0
                feasible = 0;
                return
            end
        end
            
    else
        cdi = cd(i);
        if isempty(general_w)

            ci = ci_basis(:,i);

            di = basis(cdi,0);
            if isa(Bbetai,'double')
                Bi = zeros(1,length(w));
            else
                Bi = basis(Bbetai,w)';
            end
            bi = basis(Bbeta(i,:),0)';
            % Scale to -1,1 uncertainty
            T = diag((upper-lower))/2;
            e = (upper+lower)/2;
            if nnz(Bi) == 0
                if nnz(bi)==0
                    % Basically constant + w > 0
                    if  (di+e'*ci) - norm(T*ci,1) < 0
                        error('Problem is trivially infeasible');
                        feasible = 0;
                        return
                    end
                else
                    F = F + (bi'*x + (di+e'*ci) - norm(T*ci,1) >= 0);
                end
            else
                non_zeroBirow = find(sum(abs(Bi'),2));
                zeroBirow = find(sum(abs(Bi'),2) == 0);
                if length(non_zeroBirow)>1
                    t = sdpvar(length(non_zeroBirow),1);
                    F = F + ((bi'+e'*Bi')*x + (di+e'*ci) - sum(t) >= 0) + (-t <= T(non_zeroBirow,:)*(ci+Bi'*x) <= t);
                else
                    F = F + ((bi'+e'*Bi')*x + (di+e'*ci) - T(non_zeroBirow,:)*(ci+Bi'*x) >= 0) ;
                    F = F + ((bi'+e'*Bi')*x + (di+e'*ci) + T(non_zeroBirow,:)*(ci+Bi'*x) >= 0) ;
                end
            end

        else

            lhs1 = 0;
            lhs2 = 0;
            top = 1;

            Flocal = [];
            
            if Zmodel.K.f > 0
                zeta = sdpvar(Zmodel.K.f,1);            
                lhs1 = lhs1 + Zmodel.F_struc(top:top + Zmodel.K.f-1,2:end)'*zeta;
                lhs2 = lhs2 - Zmodel.F_struc(top:top + Zmodel.K.f-1,1)'*zeta;
                top = top + Zmodel.K.f;
            end

            if Zmodel.K.l > 0
                zeta = sdpvar(Zmodel.K.l,1);
                Flocal = Flocal + (zeta >= 0);
                lhs1 = lhs1 + Zmodel.F_struc(top:top + Zmodel.K.l-1,2:end)'*zeta;
                lhs2 = lhs2 - Zmodel.F_struc(top:top + Zmodel.K.l-1,1)'*zeta;
                top = top + Zmodel.K.l;
            end

            if Zmodel.K.q(1) > 0
                for j = 1:length(Zmodel.K.q)
                    zeta = sdpvar(Zmodel.K.q(j),1);
                    if length(zeta)>2
                        Flocal = Flocal + (cone(zeta));
                    else
                        Flocal = Flocal + (zeta(2) <= zeta(1)) + (-zeta(2) <= zeta(1));
                    end
                    lhs1 = lhs1 + Zmodel.F_struc(top:top + Zmodel.K.q(j)-1,2:end)'*zeta(:);
                    lhs2 = lhs2 - Zmodel.F_struc(top:top + Zmodel.K.q(j)-1,1)'*zeta(:);
                    top = top + Zmodel.K.q(j);
                end
            end
            
            if Zmodel.K.s(1) > 0
                for j = 1:length(Zmodel.K.s)
                    zeta = sdpvar(Zmodel.K.s(j));
                    Flocal = Flocal + (zeta >= 0);
                    lhs1 = lhs1 + Zmodel.F_struc(top:top + Zmodel.K.s(j)^2-1,2:end)'*zeta(:);
                    lhs2 = lhs2 - Zmodel.F_struc(top:top + Zmodel.K.s(j)^2-1,1)'*zeta(:);
                    top = top + Zmodel.K.s(j)^2;
                end
            end

            %  if isempty(simple_w)
            ci = basis(cd(i),w);
            di = basis(cd(i),0);
            Bi = basis(Bbeta(i,:),w)';
            bi = basis(Bbeta(i,:),0)';
            if ops.robust.reducedual
                Ablhs = getbase(lhs1);
                blhs = Ablhs(:,1);
                Alhs = Ablhs(:,2:end);
                % b+A*zeta ==  Bi'*x + ci
                % A*zeta == -b +  Bi'*x + ci
                % zets ==...
                Anull = null(full(Alhs));
                zeta2 = (Alhs\(-blhs +  Bi'*x + ci))+Anull*sdpvar(size(Anull,2),1);
                lhs2 = replace(lhs2,recover(depends(lhs1)),zeta2);
                Flocal = replace(Flocal,recover(depends(lhs1)),zeta2);
            else
                Flocal = [Flocal,lhs1 == Bi'*x + ci]; 
               % F = F + (lhs1 == Bi'*x + ci);
            end
            if isa(Bi,'double') & ops.robust.reducesemiexplicit
                ops2=ops;ops2.verbose = 0;
                sol = solvesdp([Flocal],-lhs2,ops);
                if sol.problem == 0
                    F = F + (double(lhs2) >= - (bi'*x + di));
                else
                    F = F + Flocal;
                    F = F + (lhs2 >= - (bi'*x + di));
                end
            else
                F = F + Flocal;
                F = F + (lhs2 >= - (bi'*x + di));
            end
        end
    end
end
if ~isempty(simple_rows)
    F = [F, X(simple_rows)>=0];
end

function b = basis(p,w)

if isequal(w,0)
    b = getbasematrix(p,0);
else
    n = length(w);
    if  isequal(getbase(w),[spalloc(n,1,0) speye(n)])
        if 0
            b = [];
            lmi_variables = getvariables(w);
            for i = 1:length(w)
                b = [b ; getbasematrix(p,lmi_variables(i))];
            end
        else
            lmi_variables = getvariables(w);
            b = spalloc(n,length(p),0);
            [~,loc] = ismember(getvariables(w),getvariables(p));
            p_basis = getbase(p);p_basis = p_basis(:,2:end);
            used = find(loc);
            for i = 1:length(used)
                b(used(i),:) = p_basis(:,loc(used(i)))';
            end
        end
    %    if norm(b-b2)>1e-10
    %        error('sdfsdfsd')
    %    end
        
        
    else
        b = [];
        for i = 1:length(w)
            b = [b ; getbasematrix(p,getvariables(w(i)))];
        end
    end
end
b = full(b);












