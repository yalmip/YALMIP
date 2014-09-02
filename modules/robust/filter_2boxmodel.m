function [F,feasible] = filter_2boxmodel(F_xw,F_struc,x,w,allw)

% Creates robustified version of the uncertain set of linear inequalities
% s.t A(w)*x <= b(w) for all w^2< r^2

% As a first step, we figure out the radius
r = sqrt((F_struc(1,1)^2-F_struc(end,1)^2)/4);

feasible = 1;

if length(F_xw) == 0
    F = [];
    return
end

X = sdpvar(F_xw);
% Some pre-calc
x = recover(unique([getvariables(x),setdiff(getvariables(allw),getvariables(w))]));

xw = [x;w];
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
[Qs,cs,fs,dummy,nonquadratic] = vecquaddecomp(X,xw);
all_f = [];
all_c_w = [];
all_c_x = [];
all_Q_xw = [];
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

    all_f = [all_f;f];
    all_c_w = [all_c_w;c_w'];
    all_c_x = [all_c_x;sparse(c_x)];
    all_Q_xw = [all_Q_xw Q_xw'];
end
% Linear uncertain constraint is (Bbetai*x + cdi) >= 0 for all w, or
% (bi' + (Bi*w)')*x + (ci'*w + di).

F = ([]);

% To speed up the construction, compute the ci vectors for all constraints
% in one call ci_basis = [c1 c2 ...]
ci_basis = all_c_w';
for i = 1:length(all_f)
    Bi = 2*all_Q_xw(:,length(x)*(i-1)+1:length(x)*i)';
    bi = all_c_x(length(x)*(i-1)+1:length(x)*i);
    if (nnz(ci_basis(:,i))==0) & nnz(Bi)==0
        % This constraint row is constant
        F = F + (bi(:)'*x + all_f(i) >= 0);
    else
        ci = ci_basis(:,i);
        di = all_f(i);
        F = F + (bi'*x + di - r*norm(Bi'*x+ci) > 0);
    end
end
