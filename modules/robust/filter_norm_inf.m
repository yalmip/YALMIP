function [F,feasible] = filter_norm_inf(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X,Q_xx,VariableType)

if any(Zmodel.lb > Zmodel.ub)
    error('The uncertainty model is weird. Some of the lower bounds on the uncertainty are larger than the upper bounds. In other words, the feasible set for the uncertainty is empty.');    
end

feasible = 1;
lower = Zmodel.lb;
upper = Zmodel.ub;

F = ([]);

% To speed up the construction, compute the ci vectors for all constraints
% in one call ci_basis = [c1 c2 ...]
ci_basis = all_c_w';

collectA = [];
collectb = [];
collectE = [];
collectd = [];
collectC = [];
           
x = flush(x);

allbi = [];
alldi = [];
allnormT=[];
for i = 1:length(all_f)
    Bi = 2*all_Q_xw(:,length(x)*(i-1)+1:length(x)*i)';
    bi = all_c_x(length(x)*(i-1)+1:length(x)*i);
    
    if (nnz(ci_basis(:,i))==0) & nnz(Bi)==0
        % This constraint row is constant
        F = F + (X(i)>=0);
    else
        ci = ci_basis(:,i);

        di = all_f(i);
        % Scale to -1,1 uncertainty
        T = diag(sparse((upper-lower)))/2;
        e = (upper+lower)/2;
        if nnz(Bi) == 0
            if nnz(bi)==0 & nnz(Q_xx{i})==0
                % Basically constant + w > 0
                if  (di+e'*ci) - norm(T*ci,1) < 0
                    error('Problem is trivially infeasible');
                    feasible = 0;
                    return
                end
            elseif nnz(Q_xx{i})==0
                allbi = [allbi;bi'];
                alldi = [alldi;(di+e'*ci)-norm(T*ci,1) ];    
               % F = F + (bi'*x + (di+e'*ci) - norm(T*ci,1) > 0);
            else                           
                F = F + (x'*Q_xx{i}*x+bi'*x + (di+e'*ci) - norm(T*ci,1) >= 0);                
            end
        else
            % (bi' + (Bi*w)')*x + (ci'*w + di)
            % (bi' + (Bi*(e+Tw))')*x + (ci'*(e+Tw) + di). |w|<1
            % (bi'+e'*Bi')*x + (di+e'*ci) +(ci'T+x'*Bi*T')*w
            non_zeroBirow = find(sum(abs(Bi'),2));
            zeroBirow = find(sum(abs(Bi'),2) == 0);
            if length(non_zeroBirow)>1 & nnz(Q_xx{i})==0
                % This is what we are doing...
                %t = sdpvar(length(non_zeroBirow),1);
                %F = F + ((bi'+e'*Bi')*x + (di+e'*ci) - norm(T(zeroBirow,:)*ci,1)-sum(t) >= 0) + (-t < T(non_zeroBirow,:)*(ci+Bi'*x) < t);
                % However, to save over-head, we save information and
                % post-pone the generation of a constraint 
                % A*x+b+C*t>0,-t<d+E*x<t
                collectA = [collectA;sparse((bi'+e'*Bi'))];
                collectb = [collectb;(di+e'*ci) - norm(T(zeroBirow,:)*ci,1)];
                collectC = blkdiag(collectC,-sparse(ones(1,length(non_zeroBirow))));
                collectE = [collectE; T(non_zeroBirow,:)*Bi'];
                collectd = [collectd; T(non_zeroBirow,:)*ci];
                %F = F + ((bi'+e'*Bi')*x + (di+e'*ci) - norm(T(zeroBirow,:)*ci,1)-sum(t) >= 0) + (-t < T(non_zeroBirow,:)*(ci+Bi'*x) < t);
            else
                % There is only one expression involving product between x
                % and w. We explicitly construct the absolut value
                % constraint projection
                extremecase = x'*Q_xx{i}*x+(bi'+e'*Bi')*x + (di+e'*ci) - norm(T(zeroBirow,:)*ci,1)-T(non_zeroBirow,:)*(ci+Bi'*x);
                if isnumeric(extremecase)
                    if extremecase < 0
                        error('A worst-case computation reveals a trivially infeasible constraint');
                    end
                else
                    F = F + (extremecase >= 0) ;
                end
                extremecase = x'*Q_xx{i}*x+(bi'+e'*Bi')*x + (di+e'*ci) - norm(T(zeroBirow,:)*ci,1)+T(non_zeroBirow,:)*(ci+Bi'*x);
                if isnumeric(extremecase)
                    if extremecase < 0
                        error('A worst-case computation reveals a trivially infeasible constraint.');
                    end
                else
                    F = F + (extremecase >= 0) ;
                end
            end
        end
    end
end
if ~isempty(allbi)
    F = F+[allbi*x+alldi>=0];
end
if ~isempty(collectA)
    z = collectE*x + collectd;
    [U,L] = derivebounds(z);
    if all(L>=0)
        % All variables are non-negative, hence no reason to introduce an
        % epigraph variable to model abs(Ex+d)
          F = F + (collectA*x + collectb+collectC*z >= 0);
    elseif all(U<=0)
          F = F + (collectA*x + collectb+collectC*(-z) >= 0);
    else
        t = sdpvar(size(collectC,2),1);
        %  F = F + (collectA*x + collectb+collectC*t >= 0);
        %  F = F + (-t <= collectE*x + collectd <= t);
        F = F + ([-collectA -collectC;collectE -speye(length(t));-collectE -speye(length(t))]*[x;t] + [-collectb;collectd;-collectd] <= 0);
    end
end