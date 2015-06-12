function [F,feasible] = filter_normball(F_xw,Zmodel,x,w,allw,norm_p,ops,VariableType)

% Creates robustified version of the uncertain set of linear inequalities
% s.t A(w)*x <= b(w) for all norm(w,p)<r

feasible = 1;
if length(F_xw) == 0
    F = [];
    return
end

X = sdpvar(F_xw);

% Some pre-calc. Consider uncertainties not in this class as decision
% variables (i.e. other than the ones we are robustifying)
vars = [getvariables(x) depends(X) getvariables(allw)];
vars = setdiff(vars,getvariables(w));
x = recover(vars);

% Trying to cope with cases as x+norm(w,1)+norm(w,2)<1, s.t norm(w,2)<1
% In test_robust_5
%extended_in_w = intersect(VariableType.w_variables,yalmip('extvariables'));
%extended_in_F = intersect(vars,yalmip('extvariables'));
if ~(strcmp(ops.robust.auxreduce,'none') | strcmp(ops.robust.auxreduce,'affine'))
    if ~isempty(intersect(vars,intersect(VariableType.w_variables,yalmip('extvariables'))))
        F = F_xw;
        return
    end
end


% Create a bilinear decomposition of the constraints
try
xw = [x;w];
xind = find(ismembcYALMIP(getvariables(xw),getvariables(x)));
wind = find(ismembcYALMIP(getvariables(xw),getvariables(w)));
[Qs,cs,fs,dummy,nonquadratic] = vecquaddecomp(X,xw);
all_f = [];
all_c_w = [];
all_c_x = [];
all_Q_xw = [];
HigherOrder = [];
top = 1;
for i = 1:length(X)
    Q = Qs{i};
    c = cs{i};
    f = fs{i};
    if nonquadratic
        HigherOrder = [HigherOrder;i];
        %error('Constraints can be at most quadratic, with the linear term uncertain');
    else        
        Q_ww = Q(wind,wind);
        Q_xw = Q(xind,wind);
        Q_xx{top} = Q(xind,xind);
        c_x = c(xind);
        c_w = c(wind);
        
        all_f = [all_f;f];
        all_c_w = [all_c_w;c_w'];
        all_c_x = [all_c_x;sparse(c_x)];
        all_Q_xw = [all_Q_xw Q_xw'];
        top = top + 1;
    end
end
% Linear uncertain constraint is (Bbetai*x + cdi) >= 0 for all w, or
% (bi' + (Bi*w)')*x + (ci'*w + di).
catch
    % Complicated epression, so go for polynomial decomposition below
    Q = [];
end

% Currently 3 special cases implemented
if ops.verbose
    disp([' - Eliminating uncertainty using explicit maximization of ' norm_p])
end

if length(Q)> 0
    switch norm_p
        case '1-norm'
            [F,feasible] = filter_norm_1(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X,Q_xx,VariableType);
            
        case '2-norm'
            [F,feasible] = filter_norm_2(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X,Q_xx,VariableType);
            
        case 'inf-norm'
            [F,feasible] = filter_norm_inf(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X,Q_xx,VariableType);
            
        otherwise
            error('The detected norm-ball has not been implemented yet')
    end
elseif isequal(norm_p,'inf-norm')
    % Quadratic decomposition failed. Let's try a general basis if simple
    % inf-norm
    F = [];
    feasible = 1;
    HigherOrder = [];
   % [c] = coefficients([1+sum(allw);X],allw);
    for i = 1:length(X)
       % try            
            [c,v] = coefficients(X(i),allw,[1;allw]);
            % Great, affine in uncertainty            
            if isa(c(2:end),'double') && nnz(c(2:end))==0
                % Trivial case, no uncertainty in model
                F = [F, X(i) >= 0];
            else
                % c(1) +  c(2:end)'*w >= 0
                % c(1) +  c(2:end)'*(center + r*wnormalized) >= 0
                % c(1) +  c(2:end)'*center + r*c(2:end)'*wnormalized >= 0
                center = (Zmodel.lb + Zmodel.ub)/2;
                r = (Zmodel.ub - Zmodel.lb)/2;
                t = sdpvar(length(allw),1);
                F = [F, c(1) + c(2:end)'*center - sum(t) >= 0, -t <= r.*c(2:end) <= t];
            end
      %  catch
      %      HigherOrder = [HigherOrder;i];
      %  end
    end
else
    F = [];
    feasible = 1;
end

if ~isempty(HigherOrder)
    F = [F,X(HigherOrder) >= 0];
end